/*
Bullet Continuous Collision Detection and Physics Library
Copyright (c) 2003-2006 Erwin Coumans  http://continuousphysics.com/Bullet/

This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from the use of this software.
Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it freely,
subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/


#include "btSequentialImpulseConstraintSolverMt.h"

#include "LinearMath/btIDebugDraw.h"
#include "LinearMath/btMinMax.h"
#include "LinearMath/btStackAlloc.h"
#include "LinearMath/btQuickprof.h"
#include "LinearMath/btAlignedObjectArray.h"

#include "BulletCollision/NarrowPhaseCollision/btPersistentManifold.h"
#include "BulletCollision/CollisionDispatch/btUnionFind.h"

#include "BulletDynamics/ConstraintSolver/btTypedConstraint.h"
#include "BulletDynamics/Dynamics/btRigidBody.h"

#include <new>
#include <string.h> //for memset


bool btSequentialImpulseConstraintSolverMt::sAllowNestedParallelForLoops = false;  // some task schedulers don't like nested loops
int btSequentialImpulseConstraintSolverMt::sMinimumContactManifoldsForBatching = 1500;
int btSequentialImpulseConstraintSolverMt::sMinBatchSize = 40;
int btSequentialImpulseConstraintSolverMt::sMaxBatchSize = 50;

const int kUnassignedBatch = -1;
const int kNoMerge = -1;

struct btBatchedConstraintInfo
{
    int bodyIds[2];
};


struct btBatchInfo
{
    int phaseId;
    int numConstraints;
    int mergeIndex;

    btBatchInfo(int _phaseId = -1) : numConstraints(0), mergeIndex(-1), phaseId(_phaseId) {}
};


bool BatchedConstraints::validate(btConstraintArray* constraints, const btAlignedObjectArray<btSolverBody>& bodies) const
{
    //
    // validate: for debugging only. Verify coloring of bodies, that no body is touched by more than one batch in any given phase
    //
    int errors = 0;
    btAlignedObjectArray<int> bodyBatchId;
    for (int iPhase = 0; iPhase < m_phases.size(); ++iPhase)
    {
        bodyBatchId.resizeNoInitialize(0);
        bodyBatchId.resize( bodies.size(), kUnassignedBatch );
        const Range& phase = m_phases[iPhase];
        for (int iBatch = phase.begin; iBatch < phase.end; ++iBatch)
        {
            const Range& batch = m_batches[iBatch];
            for (int iiCons = batch.begin; iiCons < batch.end; ++iiCons)
            {
                int iCons = m_constraintIndices[iiCons];
                const btSolverConstraint& cons = constraints->at(iCons);
                const btSolverBody& bodyA = bodies[cons.m_solverBodyIdA];
                const btSolverBody& bodyB = bodies[cons.m_solverBodyIdB];
                if (! bodyA.internalGetInvMass().isZero())
                {
                    int thisBodyBatchId = bodyBatchId[cons.m_solverBodyIdA];
                    if (thisBodyBatchId == kUnassignedBatch)
                    {
                        bodyBatchId[cons.m_solverBodyIdA] = iBatch;
                    }
                    else if (thisBodyBatchId != iBatch)
                    {
                        btAssert( !"dynamic body is used in 2 different batches in the same phase" );
                        errors++;
                    }
                }
                if (! bodyB.internalGetInvMass().isZero())
                {
                    int thisBodyBatchId = bodyBatchId[cons.m_solverBodyIdB];
                    if (thisBodyBatchId == kUnassignedBatch)
                    {
                        bodyBatchId[cons.m_solverBodyIdB] = iBatch;
                    }
                    else if (thisBodyBatchId != iBatch)
                    {
                        btAssert( !"dynamic body is used in 2 different batches in the same phase" );
                        errors++;
                    }
                }
            }
        }
    }
    return errors == 0;
}


static void initBatchedBodyDynamicFlags(btAlignedObjectArray<bool>* outBodyDynamicFlags, const btAlignedObjectArray<btSolverBody>& bodies)
{
    BT_PROFILE("initBatchedBodyDynamicFlags");
    btAlignedObjectArray<bool>& bodyDynamicFlags = *outBodyDynamicFlags;
    bodyDynamicFlags.resizeNoInitialize(bodies.size());
    for (int i = 0; i < bodies.size(); ++i)
    {
        const btSolverBody& body = bodies[ i ];
        bodyDynamicFlags[i] = ( body.internalGetInvMass().x() > btScalar( 0 ) );
    }
}


static void initBatchedConstraintInfo(btAlignedObjectArray<btBatchedConstraintInfo>* outConInfos, btConstraintArray* constraints)
{
    BT_PROFILE("initBatchedConstraintInfo");
    btAlignedObjectArray<btBatchedConstraintInfo>& conInfos = *outConInfos;
    int numConstraints = constraints->size();
    conInfos.resizeNoInitialize(numConstraints);

    for (int i = 0; i < numConstraints; ++i)
    {
        btBatchedConstraintInfo& conInfo = conInfos[i];
        const btSolverConstraint& con = constraints->at( i );
        conInfo.bodyIds[0] = con.m_solverBodyIdA;
        conInfo.bodyIds[1] = con.m_solverBodyIdB;
    }
}


static void createBatchesForPhaseGreedy(int curPhaseId,
    btAlignedObjectArray<int>& constraintBatchIds,
    const btAlignedObjectArray<btBatchedConstraintInfo>& conInfos,
    btAlignedObjectArray<int>& curConstraints,
    const btAlignedObjectArray<bool>& bodyDynamicFlags,
    btAlignedObjectArray<int>& bodyBatchIds,
    btAlignedObjectArray<btBatchInfo>& batches,
    int minBatchSize,
    int maxBatchSize
    )
{
    BT_PROFILE("createBatchesForPhaseGreedy");
    const bool allowMerging = true;

    // mark all bodies unassigned
    const int numBodies = bodyDynamicFlags.size();
    bodyBatchIds.resize( 0 );
    bodyBatchIds.resize( numBodies, kUnassignedBatch );

    // for each unassigned constraint,
    for ( int iiCon = curConstraints.size() - 1; iiCon >= 0; --iiCon )
    {
        int iCon = curConstraints[ iiCon ];
        const btBatchedConstraintInfo& conInfo = conInfos[ iCon ];
        btAssert( constraintBatchIds[iCon] == kUnassignedBatch );
        int batchIds[ 2 ] = { kUnassignedBatch, kUnassignedBatch };
        bool isDynamic[ 2 ];
        for ( int iiBody = 0; iiBody < 2; ++iiBody )
        {
            int iBody = conInfo.bodyIds[ iiBody ];
            isDynamic[ iiBody ] = bodyDynamicFlags[iBody];
            if ( isDynamic[ iiBody ] )
            {
                int iBatch = bodyBatchIds[iBody];
                if ( allowMerging )
                {
                    // check for merge remapping
                    while ( iBatch != kUnassignedBatch && batches[ iBatch ].mergeIndex != kNoMerge )
                    {
                        iBatch = batches[ iBatch ].mergeIndex;
                    }
                }
                batchIds[ iiBody ] = iBatch;
            }
        }
        int assignBatchId = kUnassignedBatch;
        int numBodiesToBeAdded = 0;
        // if both unassigned,
        if ( batchIds[ 0 ] == kUnassignedBatch && batchIds[ 1 ] == kUnassignedBatch )
        {
            // create a new batch
            assignBatchId = batches.size();
            btBatchInfo batch;
            batch.phaseId = curPhaseId;
            batches.push_back( batch );
            btAssert( assignBatchId < batches.size() );
        }
        else if ( batchIds[ 0 ] == kUnassignedBatch )
        {
            // one is unassigned, the other not
            assignBatchId = batchIds[ 1 ];
            if ( isDynamic[ 0 ] )
            {
                numBodiesToBeAdded = 1;
            }
        }
        else if ( batchIds[ 1 ] == kUnassignedBatch )
        {
            // one is unassigned, the other not
            assignBatchId = batchIds[ 0 ];
            if ( isDynamic[ 1 ] )
            {
                numBodiesToBeAdded = 1;
            }
        }
        else if ( batchIds[ 0 ] == batchIds[ 1 ] )
        {
            // both bodies in same batch already
            // if it won't make us exceed max batch size,
            if ( batches[ batchIds[ 0 ] ].numConstraints < maxBatchSize )
            {
                assignBatchId = batchIds[ 0 ];
            }
        }
        else if ( allowMerging )
        {
            // bodies already assigned to different batches
            // we could either merge batches, or postpone until next phase
            int batch0Size = batches[ batchIds[ 0 ] ].numConstraints;
            int batch1Size = batches[ batchIds[ 1 ] ].numConstraints;
            btAssert( batches[ batchIds[ 0 ] ].phaseId == batches[ batchIds[ 1 ] ].phaseId );
            if ( batch0Size < minBatchSize && batch1Size < minBatchSize && ( batch0Size + batch1Size ) < maxBatchSize )
            {
                // merge higher index batch into lower index batch
                assignBatchId = btMin( batchIds[ 0 ], batchIds[ 1 ] );
                int otherBatch = btMax( batchIds[ 0 ], batchIds[ 1 ] );
                batches[ assignBatchId ].numConstraints += batches[ otherBatch ].numConstraints;
                btAssert( batches[ assignBatchId ].mergeIndex == kNoMerge );
                batches[ otherBatch ].numConstraints = 0;
                batches[ otherBatch ].mergeIndex = assignBatchId;
            }
        }
        // if batch has reached minBatchSize and this constraint would add a body to the batch,
        if ( assignBatchId != kUnassignedBatch && numBodiesToBeAdded > 0 && batches[ assignBatchId ].numConstraints >= minBatchSize )
        {
            // don't assign
            assignBatchId = kUnassignedBatch;
        }
        if ( assignBatchId != kUnassignedBatch )
        {
            constraintBatchIds[iCon] = assignBatchId;
            batches[ assignBatchId ].numConstraints++;
            // swap and pop
            curConstraints.swap( iiCon, curConstraints.size() - 1 );
            curConstraints.pop_back();
            // make sure bodies are assigned
            for ( int iiBody = 0; iiBody < 2; ++iiBody )
            {
                if ( isDynamic[ iiBody ] && batchIds[ iiBody ] != assignBatchId )
                {
                    int iBody = conInfo.bodyIds[ iiBody ];
                    bodyBatchIds[ iBody ] = assignBatchId;
                }
            }
        }
    }
}


static void mergeSmallBatches(btBatchInfo* batches, int iBeginBatch, int iEndBatch, int minBatchSize, int maxBatchSize)
{
    BT_PROFILE("mergeSmallBatches");
    for ( int iBatch = iEndBatch - 1; iBatch >= iBeginBatch; --iBatch )
    {
        btBatchInfo& batch = batches[ iBatch ];
        if ( batch.mergeIndex == kNoMerge && batch.numConstraints < minBatchSize )
        {
            for ( int iDestBatch = iBatch - 1; iDestBatch >= iBeginBatch; --iDestBatch )
            {
                btBatchInfo& destBatch = batches[ iDestBatch ];
                if ( destBatch.mergeIndex == kNoMerge && ( destBatch.numConstraints + batch.numConstraints ) < maxBatchSize )
                {
                    destBatch.numConstraints += batch.numConstraints;
                    batch.numConstraints = 0;
                    batch.mergeIndex = iDestBatch;
                    break;
                }
            }
        }
    }
    // flatten mergeIndexes
    // e.g. in case where A was merged into B and then B was merged into C, we need A to point to C instead of B
    // Note: loop goes forward through batches because batches always merge from higher indexes to lower,
    //     so by going from low to high it reduces the amount of trail-following
    for ( int iBatch = iBeginBatch; iBatch < iEndBatch; ++iBatch )
    {
        btBatchInfo& batch = batches[ iBatch ];
        if ( batch.mergeIndex != kNoMerge )
        {
            int iMergeDest = batches[ batch.mergeIndex ].mergeIndex;
            // follow trail of merges to the end
            while ( iMergeDest != kNoMerge )
            {
                int iNext = batches[ iMergeDest ].mergeIndex;
                if ( iNext == kNoMerge )
                {
                    batch.mergeIndex = iMergeDest;
                    break;
                }
                iMergeDest = iNext;
            }
        }
    }
}


static void updateConstraintBatchIdsForMerges(int* constraintBatchIds, int numConstraints, const btBatchInfo* batches, int numBatches)
{
    BT_PROFILE("updateConstraintBatchIdsForMerges");
    // update batchIds to account for merges
    for (int i = 0; i < numConstraints; ++i)
    {
        int iBatch = constraintBatchIds[i];
        btAssert(iBatch < numBatches);
        // if this constraint references a batch that was merged into another batch
        if (batches[iBatch].mergeIndex != kNoMerge)
        {
            // update batchId
            constraintBatchIds[i] = batches[iBatch].mergeIndex;
        }
    }
}


struct UpdateConstraintBatchIdsForMergesLoop : public btIParallelForBody
{
    int* m_constraintBatchIds;
    const btBatchInfo* m_batches;
    int m_numBatches;

    UpdateConstraintBatchIdsForMergesLoop( int* constraintBatchIds, const btBatchInfo* batches, int numBatches )
    {
        m_constraintBatchIds = constraintBatchIds;
        m_batches = batches;
        m_numBatches = numBatches;
    }
    void forLoop( int iBegin, int iEnd ) const BT_OVERRIDE
    {
        BT_PROFILE( "UpdateConstraintBatchIdsForMergesLoop" );
        updateConstraintBatchIdsForMerges( m_constraintBatchIds + iBegin, iEnd - iBegin, m_batches, m_numBatches );
    }
};


static void updateConstraintBatchIdsForMergesMt(int* constraintBatchIds, int numConstraints, const btBatchInfo* batches, int numBatches)
{
    BT_PROFILE( "updateConstraintBatchIdsForMergesMt" );
    UpdateConstraintBatchIdsForMergesLoop loop(constraintBatchIds, batches, numBatches);
    int grainSize = 800;
    btParallelFor(0, numConstraints, grainSize, loop);
}


inline bool BatchCompare(const BatchedConstraints::Range& a, const BatchedConstraints::Range& b)
{
    int lenA = a.end - a.begin;
    int lenB = b.end - b.begin;
    return lenA > lenB;
}


static void writeOutBatches(BatchedConstraints* bc,
    const int* constraintBatchIds,
    int numConstraints,
    const btBatchInfo* batches,
    int* batchWork,
    int maxNumBatchesPerPhase,
    int numPhases
)
{
    BT_PROFILE("writeOutBatches");
    typedef BatchedConstraints::Range Range;
    bc->m_constraintIndices.reserve( numConstraints );
    bc->m_batches.resizeNoInitialize( 0 );
    bc->m_phases.resizeNoInitialize( 0 );

    int maxNumBatches = numPhases * maxNumBatchesPerPhase;
    {
        int* constraintIdPerBatch = batchWork;  // for each batch, keep an index into the next available slot in the m_constraintIndices array
        int iConstraint = 0;
        int curPhaseBegin = 0;
        int curPhaseId = 0;
        for (int iPhase = 0; iPhase < numPhases; ++iPhase)
        {
            int iBegin = iPhase * maxNumBatchesPerPhase;
            int iEnd = iBegin + maxNumBatchesPerPhase;
            for ( int i = iBegin; i < iEnd; ++i )
            {
                const btBatchInfo& batch = batches[ i ];
                int curBatchBegin = iConstraint;
                constraintIdPerBatch[ i ] = curBatchBegin;  // record the start of each batch in m_constraintIndices array
                int numConstraints = batch.numConstraints;
                iConstraint += numConstraints;
                if ( numConstraints > 0 )
                {
                    // if new phase is starting
                    if ( batch.phaseId != curPhaseId )
                    {
                        // output the previous phase
                        bc->m_phases.push_back( Range( curPhaseBegin, bc->m_batches.size() ) );
                        curPhaseBegin = bc->m_batches.size();
                        curPhaseId = batch.phaseId;
                    }
                    bc->m_batches.push_back( Range( curBatchBegin, iConstraint ) );
                }
            }
        }
        if (bc->m_batches.size() > curPhaseBegin)
        {
            // output last phase
            bc->m_phases.push_back( Range( curPhaseBegin, bc->m_batches.size() ) );
        }

        btAssert(iConstraint == numConstraints);
        bc->m_constraintIndices.resizeNoInitialize( numConstraints );
        for ( int iCon = 0; iCon < numConstraints; ++iCon )
        {
            int iBatch = constraintBatchIds[iCon];
            //int iMergeBatch = batches[iBatch].mergeIndex;
            //if (iMergeBatch != kNoMerge)
            //{
            //    iBatch = iMergeBatch;
            //}
            int iDestCon = constraintIdPerBatch[iBatch];
            constraintIdPerBatch[iBatch] = iDestCon + 1;
            bc->m_constraintIndices[iDestCon] = iCon;
        }
    }
    // for each phase
    for (int iPhase = 0; iPhase < bc->m_phases.size(); ++iPhase)
    {
        // sort the batches from largest to smallest (can be helpful to some task schedulers)
        const Range& curBatches = bc->m_phases[iPhase];
        bc->m_batches.quickSortInternal(BatchCompare, curBatches.begin, curBatches.end-1);
    }
    bc->m_phaseOrder.resize(bc->m_phases.size());
    for (int i = 0; i < bc->m_phases.size(); ++i)
    {
        bc->m_phaseOrder[i] = i;
    }
}


static int getSumNumConstraints(const btAlignedObjectArray<btBatchInfo>& batches)
{
    int sum = 0;
    for (int i = 0; i < batches.size(); ++i)
    {
        sum += batches[i].numConstraints;
    }
    return sum;
}


//
// setupBatchesGreedyWithMerging
//
// This method of generating batches produces phases that are uneven in
// the number of batches and constraints in each. The first phase generated
// has the most constraints and batches, and then subsequent phases tend to have
// fewer and fewer batches/constraints with the last few phases having just
// a single batch each with just a handful of constraints.
// This is rather inefficent for thread scheduling. The first couple of phases
// can run in parallel nicely, but then the parallelism declines and the last few
// phases are serial.
// Also the generation of batches was not particularly fast and due to the sequential
// nature of how it works, not easily parallelisable.
//
void setupBatchesGreedyWithMerging(
    BatchedConstraints* batchedConstraints,
    btConstraintArray* constraints,
    const btAlignedObjectArray<btSolverBody>& bodies,
    int minBatchSize,
    int maxBatchSize
    )
{
    BT_PROFILE("setupBatchesGreedyWithMerging");
    int numConstraints = constraints->size();

    const bool allowMerging = false;

    // TODO: get rid of all these dynamic allocations
    btAlignedObjectArray<bool> bodyDynamicFlags;
    initBatchedBodyDynamicFlags(&bodyDynamicFlags, bodies);

    btAlignedObjectArray<btBatchedConstraintInfo> conInfos;
    initBatchedConstraintInfo(&conInfos, constraints);

    btAlignedObjectArray<btBatchInfo> batches;
    batches.reserve( numConstraints/4 );
    btAlignedObjectArray<int> curConstraints;
    curConstraints.resizeNoInitialize(numConstraints);
    for (int i = 0; i < numConstraints; ++i)
    {
        curConstraints[i] = i;
    }

    btAlignedObjectArray<int> bodyBatchIds;
    bodyBatchIds.resize(bodies.size(), kUnassignedBatch);

    btAlignedObjectArray<int> constraintBatchIds;
    constraintBatchIds.resize(numConstraints, kUnassignedBatch);

    int curPhaseId = 0;
    while (curConstraints.size() > 0)
    {
        // begin new phase
        int curPhaseBegin = batches.size();
        createBatchesForPhaseGreedy(curPhaseId,
            constraintBatchIds,
            conInfos,
            curConstraints,
            bodyDynamicFlags,
            bodyBatchIds,
            batches,
            minBatchSize,
            maxBatchSize
        );

        btAssert(getSumNumConstraints(batches) == (numConstraints - curConstraints.size()));
        // merge small batches into others
        mergeSmallBatches(&batches[0], curPhaseBegin, batches.size(), minBatchSize, maxBatchSize);

        btAssert(getSumNumConstraints(batches) == (numConstraints - curConstraints.size()));
        curPhaseId++;
    }
    // all constraints have been assigned a batchId
    updateConstraintBatchIdsForMergesMt(&constraintBatchIds[0], numConstraints, &batches[0], batches.size());

    btAlignedObjectArray<int> batchWork;
    batchWork.resizeNoInitialize(batches.size());

    //writeOutBatches(batchedConstraints, &constraintBatchIds[0], numConstraints, &batches[0], &batchWork[0], batches.size());
    btAssert(batchedConstraints->validate(constraints, bodies));
}


void initTableRecursive(btAlignedObjectArray<int>& table, int numGroups, int x, int y, int step, int phase)
{
    btAssert(step>=1);
    if (step == 1)
    {
        table[x + y*numGroups] = phase;
        table[y + x*numGroups] = phase;
    }
    else
    {
        int newStep = step/2;
        initTableRecursive(table, numGroups, x,         y,         newStep, phase);
        initTableRecursive(table, numGroups, x+newStep, y,         newStep, phase+newStep);
        initTableRecursive(table, numGroups, x,         y+newStep, newStep, phase+newStep);
        initTableRecursive(table, numGroups, x+newStep, y+newStep, newStep, phase);
    }
}


static bool validatePhaseTable(const btAlignedObjectArray<int>& table, int numGroups)
{
    // validate
    bool valid = true;
    btAlignedObjectArray<bool> phaseFound;
    for ( int x = 0; x < numGroups; ++x )
    {
        phaseFound.resize(0);
        phaseFound.resize( numGroups, false );
        for ( int y = 0; y < numGroups; ++y )
        {
            int iPhase = table[x + y*numGroups];
            if (iPhase != table[y + x*numGroups])
            {
                btAssert(!"should be symmetric");
                valid = false;
            }
            if (phaseFound[iPhase])
            {
                btAssert(!"each phase should only appear once per column/once per row");
                valid = false;
            }
            phaseFound[iPhase] = true;
        }
    }
    return valid;
}


static void initBatchPhaseTable(btAlignedObjectArray<int>* phaseTable, int numGroups)
{
    btAlignedObjectArray<int>& table = *phaseTable;
    table.resizeNoInitialize(numGroups*numGroups);
    initTableRecursive(table, numGroups, 0, 0, numGroups, 0);
    btAssert(validatePhaseTable(table, numGroups));
}


//
// setupMtSimple
//
// For this method of batching constraints, each body is arbitrarily assigned to a 'group'
// of which there are a fixed number (numGroups). For each constraint, we find the group for
// each body of the constraint, and the 2 group ids are used to completely determine which
// batch and phase the constraint is assigned to.
// A more detailed description can be found in:
//   bullet3/docs/GPU_rigidbody_using_OpenCL.pdf (page 16)
//
// Pros:
//   - very fast to generate batches -- can use a table lookup
//   - could generate batches in parallel if needed
//
// Cons:
//   - need roughly twice as many phases as there are batches per phase (and number of phases must be a power of 2)
//   - size of batches is uneven
//
// The problem with having many phases is that the worker threads need to be synchronized at the end of
// each phase which is a relatively slow operation. And having uneven batches means many of the worker threads
// are going to be spending much of their time waiting for the one thread with the biggest batch to finish working.
//
static void setupMtSimple(
    BatchedConstraints* batchedConstraints,
    btConstraintArray* constraints,
    const btAlignedObjectArray<btSolverBody>& bodies,
    int minBatchSize,
    int maxBatchSize
    )
{
    BT_PROFILE("setupMtSimple");
    const int numGroups = 8;  // must be a power of 2
    const int groupMask = numGroups - 1;
    btAlignedObjectArray<int> groupPhaseTable;
    initBatchPhaseTable(&groupPhaseTable, numGroups);

    btAlignedObjectArray<btBatchInfo> batches;
    batches.reserve(numGroups + (numGroups-1)*numGroups/2);
    btAlignedObjectArray<int> groupBatchTable;
    groupBatchTable.resizeNoInitialize(numGroups*numGroups);
    for ( int phase = 0; phase < numGroups; ++phase )
    {
        for ( int y = 0; y < numGroups; ++y )
        {
            for ( int x = y; x < numGroups; ++x )
            {
                int iPhase = groupPhaseTable[ x + y*numGroups ];
                if ( iPhase == phase )
                {
                    groupBatchTable[x + y*numGroups] = batches.size();
                    groupBatchTable[y + x*numGroups] = batches.size();
                    btBatchInfo& batch = batches.expandNonInitializing();
                    batch.mergeIndex = -1;
                    batch.numConstraints = 0;
                    batch.phaseId = phase;
                }
            }
        }
    }

    int numConstraints = constraints->size();

    // TODO: get rid of all these dynamic allocations
    btAlignedObjectArray<btBatchedConstraintInfo> conInfos;
    initBatchedConstraintInfo(&conInfos, constraints);

    btAlignedObjectArray<int> constraintBatchIds;
    constraintBatchIds.resize(numConstraints, kUnassignedBatch);
    for (int iCon = 0; iCon < numConstraints; ++iCon)
    {
        const btBatchedConstraintInfo& conInfo = conInfos[iCon];
        int groupA = conInfo.bodyIds[0] & groupMask;
        int groupB = conInfo.bodyIds[1] & groupMask;
        int iBatch = groupBatchTable[groupA + groupB*numGroups];
        constraintBatchIds[iCon] = iBatch;
        batches[iBatch].numConstraints++;
    }

    btAlignedObjectArray<int> batchWork;
    batchWork.resizeNoInitialize(batches.size());

    //writeOutBatches(batchedConstraints, &constraintBatchIds[0], numConstraints, &batches[0], &batchWork[0], batches.size());
    btAssert(batchedConstraints->validate(constraints, bodies));
}


static void mergeSmallPhases(btAlignedObjectArray<int>* inoutGroupPhaseTable, int* inoutNumPhases, int targetNumPhases, const btAlignedObjectArray<btBatchedConstraintInfo>& conInfos)
{
    BT_PROFILE("mergeSmallPhases");
    btAlignedObjectArray<int>& groupPhaseTable = *inoutGroupPhaseTable;
    const int numGroups = *inoutNumPhases;
    const int groupMask = numGroups - 1;
    int numPhases = numGroups;

    // count the number of constraints in each phase
    btAlignedObjectArray<int> numConstraintsPerPhase;
    numConstraintsPerPhase.resize( numPhases, 0 );
    int numConstraints = conInfos.size();
    for ( int iCon = 0; iCon < numConstraints; ++iCon )
    {
        const btBatchedConstraintInfo& conInfo = conInfos[ iCon ];
        int body0 = conInfo.bodyIds[ 0 ];
        int body1 = conInfo.bodyIds[ 1 ];
        int groupA = body0 & groupMask;
        int groupB = body1 & groupMask;
        int iPhase = groupPhaseTable[ groupA + groupB*numGroups ];
        numConstraintsPerPhase[ iPhase ]++;
    }

    // merge smallest phases together until we have more evenly distributed sizes
    btAlignedObjectArray<int> phaseMergeIds;
    phaseMergeIds.resize( numPhases );
    for ( int i = 0; i < phaseMergeIds.size(); ++i )
    {
        phaseMergeIds[ i ] = i;
    }

    int maxPhaseSize = 0;
    for ( int i = 0; i < numPhases; ++i )
    {
        maxPhaseSize = btMax( maxPhaseSize, numConstraintsPerPhase[ i ] );
    }
    maxPhaseSize *= 2;

    int numNonZeroPhases = numPhases;
    while ( numNonZeroPhases > targetNumPhases )
    {
        // find the smallest phase
        int bestPhaseSize = numConstraints;
        int iSmallestPhase = -1;
        for ( int i = 0; i < numPhases; ++i )
        {
            int phaseSize = numConstraintsPerPhase[ i ];
            if ( phaseSize > 0 && phaseSize < bestPhaseSize )
            {
                bestPhaseSize = phaseSize;
                iSmallestPhase = i;
            }
        }
        // find the second-smallest phase
        bestPhaseSize = numConstraints;
        int iNextSmallestPhase = -1;
        for ( int i = 0; i < numPhases; ++i )
        {
            int phaseSize = numConstraintsPerPhase[ i ];
            if ( i != iSmallestPhase && phaseSize > 0 && phaseSize < bestPhaseSize )
            {
                bestPhaseSize = phaseSize;
                iNextSmallestPhase = i;
            }
        }
        if ( iNextSmallestPhase == -1 || numConstraintsPerPhase[ iSmallestPhase ] + numConstraintsPerPhase[ iNextSmallestPhase ] > maxPhaseSize )
        {
            break;
        }
        // merge smallest into next smallest phase
        numConstraintsPerPhase[ iNextSmallestPhase ] += numConstraintsPerPhase[ iSmallestPhase ];
        numConstraintsPerPhase[ iSmallestPhase ] = 0;
        btAssert( phaseMergeIds[ iSmallestPhase ] == iSmallestPhase );
        phaseMergeIds[ iSmallestPhase ] = iNextSmallestPhase;
        numNonZeroPhases--;
    }
    // flatten merge ids
    for ( int i = 0; i < phaseMergeIds.size(); ++i )
    {
        int iMergeDest = i;
        while ( phaseMergeIds[ iMergeDest ] != iMergeDest )
        {
            iMergeDest = phaseMergeIds[ iMergeDest ];
        }
        phaseMergeIds[ i ] = iMergeDest;
    }
    // fill in holes so that all valid phases are contiguously numbered
    int iPhaseSrc = phaseMergeIds.size() - 1;
    int iPhaseDst = 0;
    while ( true )
    {
        // find highest filled slot
        while ( iPhaseSrc > 0 && numConstraintsPerPhase[ iPhaseSrc ] == 0 )
        {
            iPhaseSrc--;
        }
        // find lowest empty slot
        while ( iPhaseDst < iPhaseSrc && numConstraintsPerPhase[ iPhaseDst ] != 0 )
        {
            iPhaseDst++;
        }
        // if src and dest crossed, we're done
        if ( iPhaseDst >= iPhaseSrc )
        {
            break;
        }
        // move src to dst
        btAssert( numConstraintsPerPhase[ iPhaseDst ] == 0 );
        numConstraintsPerPhase[ iPhaseDst ] = numConstraintsPerPhase[ iPhaseSrc ];
        numConstraintsPerPhase[ iPhaseSrc ] = 0;
        for ( int i = 0; i < phaseMergeIds.size(); ++i )
        {
            if ( phaseMergeIds[ i ] == iPhaseSrc )
            {
                phaseMergeIds[ i ] = iPhaseDst;
            }
        }
    }
    *inoutNumPhases = numNonZeroPhases;
    // remap groupPhaseTable according to merged phases
    for ( int i = 0; i < groupPhaseTable.size(); ++i )
    {
        int oldPhase = groupPhaseTable[ i ];
        int newPhase = phaseMergeIds[ oldPhase ];
        groupPhaseTable[ i ] = newPhase;
    }
}


struct CreateBatchesParams
{
    int* constraintBatchIds;
    btBatchInfo* batches;
    const char* constraintPhaseIds;
    const char* phaseMappingTable;
    const int* groupPhaseTable;
    const btBatchedConstraintInfo* conInfos;
    const bool* bodyDynamicFlags;
    int numBodies;
    int numConstraints;
    int numGroups;
    int minNumBatchesPerPhase;
    int maxNumBatchesPerPhase;
    int minBatchSize;
    int maxBatchSize;
    CreateBatchesWork* workArray;

    CreateBatchesParams()
    {
        memset(this, 0, sizeof(*this));
    }
};


static void createBatchesForPhase(int iPhase, const CreateBatchesParams& params)
{
    BT_PROFILE("createBatchesForPhase");
    int* constraintBatchIds = params.constraintBatchIds;
    const btBatchedConstraintInfo* conInfos = params.conInfos;
    const bool* bodyDynamicFlags = params.bodyDynamicFlags;
    btBatchInfo* batches = params.batches;
    const int numGroups = params.numGroups;
    const int groupMask = numGroups - 1;
    const int minNumBatchesPerPhase = params.minNumBatchesPerPhase;
    const int maxNumBatchesPerPhase = params.maxNumBatchesPerPhase;
    const int minBatchSize = params.minBatchSize;
    const int maxBatchSize = params.maxBatchSize;
    const int numBodies = params.numBodies;
    const int numConstraints = params.numConstraints;
    btAssert( numBodies > 0 );
    btAssert( numConstraints > 0 );

    int iThread = btGetCurrentThreadIndex();
    CreateBatchesWork& work = params.workArray[iThread];

    btAlignedObjectArray<int>& curConstraints = work.m_curConstraints;
    curConstraints.reserve(numConstraints);

    btAlignedObjectArray<int>& curFlexConstraints = work.m_curFlexConstraints;  // flex constraints involve one non-dynamic body, so if the dynamic body is unassigned, it can be added to any batch
    curFlexConstraints.reserve(numConstraints);

    btAlignedObjectArray<int>& bodyBatchIds = work.m_bodyBatchIds;

    btUnionFind& unionFind = work.m_unionFind;

    unionFind.reset( numBodies );
    curConstraints.resize( 0 );
    curFlexConstraints.resize( 0 );
    if (params.constraintPhaseIds)
    {
        const char* constraintPhaseIds = params.constraintPhaseIds;
        const char* phaseMappingTable = params.phaseMappingTable;
        for ( int iCon = 0; iCon < numConstraints; ++iCon )
        {
            char unmappedPhase = constraintPhaseIds[iCon];
            btAssert( unmappedPhase < params.numGroups );
            char mappedPhase = phaseMappingTable[ unmappedPhase ];
            if ( iPhase == mappedPhase )
            {
                const btBatchedConstraintInfo& conInfo = conInfos[ iCon ];
                int body0 = conInfo.bodyIds[ 0 ];
                int body1 = conInfo.bodyIds[ 1 ];
                if ( bodyDynamicFlags[ body0 ] && bodyDynamicFlags[ body1 ] )
                {
                    unionFind.unite( body0, body1 );
                    curConstraints.push_back( iCon );
                }
                else
                {
                    curFlexConstraints.push_back( iCon );
                }
            }
        }
    }
    else if (params.groupPhaseTable)
    {
        const int* groupPhaseTable = params.groupPhaseTable;
        for ( int iCon = 0; iCon < numConstraints; ++iCon )
        {
            const btBatchedConstraintInfo& conInfo = conInfos[ iCon ];
            int body0 = conInfo.bodyIds[ 0 ];
            int body1 = conInfo.bodyIds[ 1 ];
            int groupA = body0 & groupMask;
            int groupB = body1 & groupMask;
            if ( iPhase == groupPhaseTable[ groupA + groupB*numGroups ] )
            {
                if ( bodyDynamicFlags[ body0 ] && bodyDynamicFlags[ body1 ] )
                {
                    unionFind.unite( body0, body1 );
                    curConstraints.push_back( iCon );
                }
                else
                {
                    curFlexConstraints.push_back( iCon );
                }
            }
        }
    }
    bodyBatchIds.resize( 0 );
    bodyBatchIds.resize( numBodies, kUnassignedBatch );
    int iBatchBegin = iPhase * maxNumBatchesPerPhase;
    int curBatch = iBatchBegin;
    for (int iBatch = iBatchBegin; iBatch < iBatchBegin + maxNumBatchesPerPhase; ++iBatch)
    {
        btBatchInfo& batch = batches[ iBatch ];
        batch = btBatchInfo(-1);
    }
    for ( int iiCon = 0; iiCon < curConstraints.size(); ++iiCon )
    {
        int iCon = curConstraints[ iiCon ];
        const btBatchedConstraintInfo& conInfo = conInfos[ iCon ];
        int body = conInfo.bodyIds[ 0 ];
        btAssert( bodyDynamicFlags[ body ] );
        int island = unionFind.find( body );
        int iBatch = bodyBatchIds[ island ];
        if ( iBatch == kUnassignedBatch )
        {
            iBatch = curBatch++;
            bodyBatchIds[ island ] = iBatch;
            batches[ iBatch ].phaseId = iPhase;
        }
        constraintBatchIds[ iCon ] = iBatch;
        btBatchInfo& batch = batches[ iBatch ];
        batch.numConstraints++;
    }
    // assign any flex constraints that belong to already created batch
    for ( int iiCon = curFlexConstraints.size() - 1; iiCon >= 0; --iiCon )
    {
        int iCon = curFlexConstraints[ iiCon ];
        const btBatchedConstraintInfo& conInfo = conInfos[ iCon ];
        int body0 = conInfo.bodyIds[ 0 ];
        int body1 = conInfo.bodyIds[ 1 ];
        int body = bodyDynamicFlags[ body0 ] ? body0 : body1;
        btAssert( bodyDynamicFlags[ body ] );
        int island = unionFind.find( body );
        int iBatch = bodyBatchIds[ island ];
        // if dynamic body has already been assigned a batch,
        if ( iBatch != kUnassignedBatch )
        {
            // add it to the batch
            constraintBatchIds[ iCon ] = iBatch;
            btBatchInfo& batch = batches[ iBatch ];
            batch.numConstraints++;
            // swap and pop
            curFlexConstraints.swap( iiCon, curFlexConstraints.size() - 1 );
            curFlexConstraints.pop_back();
        }
    }
    mergeSmallBatches( batches, iBatchBegin, curBatch, minBatchSize, maxBatchSize );
    // if we have leftover flex constraints,
    if ( curFlexConstraints.size() > 0 )
    {
        // decide how to distribute them
        int numBatchesInPhase = 0;
        int numConstraintsInPhase = 0;
        for ( int iBatch = iBatchBegin; iBatch < curBatch; ++iBatch )
        {
            const btBatchInfo& batch = batches[ iBatch ];
            if ( batch.numConstraints > 0 )
            {
                numBatchesInPhase++;
                numConstraintsInPhase += batch.numConstraints;
            }
        }
        numConstraintsInPhase += curFlexConstraints.size();
        if ( numBatchesInPhase < minNumBatchesPerPhase )
        {
            int iNewBatch = kUnassignedBatch;
            // generate more batches
            while ( curFlexConstraints.size() > 0 )
            {
                int iCon = curFlexConstraints[ curFlexConstraints.size() - 1 ];
                curFlexConstraints.pop_back();

                const btBatchedConstraintInfo& conInfo = conInfos[ iCon ];
                int body0 = conInfo.bodyIds[ 0 ];
                int body1 = conInfo.bodyIds[ 1 ];
                int body = bodyDynamicFlags[ body0 ] ? body0 : body1;
                btAssert( bodyDynamicFlags[ body ] );
                int iBatch = bodyBatchIds[ body ];
                if ( iBatch == kUnassignedBatch )
                {
                    // create a new batch
                    if ( iNewBatch == kUnassignedBatch )
                    {
                        iNewBatch = curBatch++;
                    }
                    iBatch = iNewBatch;
                    bodyBatchIds[ body ] = iBatch;
                    batches[ iBatch ].phaseId = iPhase;
                }
                btAssert( iBatch != kUnassignedBatch );
                constraintBatchIds[ iCon ] = iBatch;
                btBatchInfo& batch = batches[ iBatch ];
                batch.numConstraints++;
            }
        }
        if ( numBatchesInPhase > 0 )
        {
            // try to even out the sizes of existing batches
            int avgNumConstraintsInBatch = numConstraintsInPhase / numBatchesInPhase;
            while ( curFlexConstraints.size() > 0 )
            {
                // find smallest batch
                int iSmallestBatch = kUnassignedBatch;
                int smallestBatchSize = numConstraints;
                for ( int iBatch = iBatchBegin; iBatch < curBatch; ++iBatch )
                {
                    int batchSize = batches[ iBatch ].numConstraints;
                    if ( batchSize > 0 && batchSize < smallestBatchSize )
                    {
                        smallestBatchSize = batchSize;
                        iSmallestBatch = iBatch;
                    }
                }
                btAssert( iSmallestBatch != kUnassignedBatch );
                for ( int iiCon = curFlexConstraints.size() - 1; iiCon >= 0; --iiCon )
                {
                    int iCon = curFlexConstraints[ iiCon ];
                    curFlexConstraints.pop_back();
                    const btBatchedConstraintInfo& conInfo = conInfos[ iCon ];
                    int body0 = conInfo.bodyIds[ 0 ];
                    int body1 = conInfo.bodyIds[ 1 ];
                    int body = bodyDynamicFlags[ body0 ] ? body0 : body1;
                    btAssert( bodyDynamicFlags[ body ] );
                    // there could be multiple flex constraints on a given body, so we always need to check bodyBatchIds array
                    int iBatch = bodyBatchIds[ body ];
                    // if dynamic body has not already been assigned a batch,
                    if ( iBatch == kUnassignedBatch )
                    {
                        // add it to the smallest batch
                        bodyBatchIds[ body ] = iSmallestBatch;
                        iBatch = iSmallestBatch;
                    }
                    // add constraint to batch
                    constraintBatchIds[ iCon ] = iBatch;
                    btBatchInfo& batch = batches[ iBatch ];
                    batch.numConstraints++;
                    if ( iBatch == iSmallestBatch )
                    {
                        smallestBatchSize = batch.numConstraints;
                        // if batch is no longer the smallest,
                        if ( smallestBatchSize >= avgNumConstraintsInBatch )
                        {
                            break;
                        }
                    }
                }
            }
        }
    }
    btAssert( curFlexConstraints.size() == 0 );
}


struct CreateBatchesForPhaseLoop : public btIParallelForBody
{
    const CreateBatchesParams* m_params;

    CreateBatchesForPhaseLoop( const CreateBatchesParams& params )
    {
        m_params = &params;
    }
    void forLoop( int iBegin, int iEnd ) const BT_OVERRIDE
    {
        BT_PROFILE( "CreateBatchesForPhaseLoop" );
        for ( int iPhase = iBegin; iPhase < iEnd; ++iPhase )
        {
            createBatchesForPhase(iPhase, *m_params);
        }
    }
};


//
// setupMt
//
// This method is a hybrid of the setupMtSimple method and the greedyWithMerges method.
//
// It works as follows:
//  - First it distributes all constraints into 16 different phases using the table-lookup
//      method from the setupMtSimple method.
//  - Then it merges the smallest phases together repeatedly, until there are around 12 roughly
//      equal size phases.
//  - Then each of the phases can be processed in parallel:
//      - constraints are gathered into a list of regular constraints and a list of "flexible" constraints
//          (flexible constraints only involve a single dynamic body, and therefore if the body is unassigned
//          the constraint can be added to any batch)
//      - the regular constraints are grouped into batches in greedy fashion
//      - small batches are merged together to form larger ones
//      - flexible constraints are then added to the smallest batches to make the batches as uniformly sized as possible
//
// Pros:
//   - usually generates 8 evenly sized phases each with a good number of evenly sized batches
//   - each of the ~8 phases is processed in parallel
//
// Cons:
//   - complicated
//   - tuned for a single use case (~1000 bodies in a stack, ~20K contact points)
//
static void setupMt(
    BatchedConstraints* batchedConstraints,
    btConstraintArray* constraints,
    const btAlignedObjectArray<btSolverBody>& bodies,
    int minBatchSize,
    int maxBatchSize,
    CreateBatchesWork* workArray
    )
{
    BT_PROFILE("setupMt");
    const int numGroups = 16;
    const int groupMask = numGroups - 1;
    int numPhases = numGroups;
    btAlignedObjectArray<int> groupPhaseTable;
    initBatchPhaseTable(&groupPhaseTable, numGroups);

    int numConstraints = constraints->size();

    int maxNumBatchesPerPhase = bodies.size();
    int minNumBatchesPerPhase = 16;

    btAlignedObjectArray<btBatchInfo> batches;
    batches.resize(maxNumBatchesPerPhase * numPhases);

    btAlignedObjectArray<btBatchedConstraintInfo> conInfos;
    initBatchedConstraintInfo(&conInfos, constraints);

    mergeSmallPhases( &groupPhaseTable, &numPhases, 12, conInfos );

    btAlignedObjectArray<bool> bodyDynamicFlags;
    initBatchedBodyDynamicFlags(&bodyDynamicFlags, bodies);

    btAlignedObjectArray<int> constraintBatchIds;
    constraintBatchIds.resize(numConstraints, kUnassignedBatch);

    CreateBatchesParams params;
    params.constraintBatchIds = &constraintBatchIds[0];
    params.batches = &batches[0];
    params.bodyDynamicFlags = &bodyDynamicFlags[0];
    params.conInfos = &conInfos[0];
    params.groupPhaseTable = &groupPhaseTable[0];
    params.numGroups = numGroups;
    params.numBodies = bodies.size();
    params.numConstraints = numConstraints;
    params.minNumBatchesPerPhase = minNumBatchesPerPhase;
    params.maxNumBatchesPerPhase = maxNumBatchesPerPhase;
    params.minBatchSize = minBatchSize;
    params.maxBatchSize = maxBatchSize;
    params.workArray = workArray;

    if (true)
    {
        // parallel batch creation (deterministic)
        CreateBatchesForPhaseLoop loop(params);
        btParallelFor( 0, numPhases, 1, loop );
    }
    else
    {
        // sequential batch creation
        for ( int iPhase = 0; iPhase < numPhases; ++iPhase )
        {
            createBatchesForPhase( iPhase, params );
        }
    }
    // all constraints have been assigned a batchId
    updateConstraintBatchIdsForMergesMt(&constraintBatchIds[0], numConstraints, &batches[0], batches.size());

    btAlignedObjectArray<int> batchWork;
    batchWork.resizeNoInitialize(batches.size());

    writeOutBatches(batchedConstraints, &constraintBatchIds[0], numConstraints, &batches[0], &batchWork[0], maxNumBatchesPerPhase, numPhases);
    btAssert(batchedConstraints->validate(constraints, bodies));
}


struct FaceMap
{
    // defines a planar uv-mapping over a face
    float m_u;
    float m_v;
    btVector3 m_uVec;
    btVector3 m_vVec;
    btVector3 m_normal; // face normal
};


static const FaceMap gCubeFaceMaps[ 8 ] =
{
    {0.5, 0.5, btVector3(0, 1, 0), btVector3(0,  0,  1), btVector3( 1,  0,  0)},  // face0 +X
    {0.5, 1.5, btVector3(1, 0, 0), btVector3(0,  0, -1), btVector3( 0,  1,  0)},  // face0 +Y
    {0.5, 2.5, btVector3(1, 0, 0), btVector3(0,  1,  0), btVector3( 0,  0,  1)},  // face0 +Z
    {0.5, 3.5, btVector3(0, 1, 0), btVector3(0,  0, -1), btVector3(-1,  0,  0)},  // face0 -X
    {0.5, 4.5, btVector3(1, 0, 0), btVector3(0,  0,  1), btVector3( 0, -1,  0)},  // face0 -Y
    {0.5, 5.5, btVector3(1, 0, 0), btVector3(0, -1,  0), btVector3( 0,  0, -1)},  // face0 -Z
};


static btVector3 normalizedDirectionToCubeUv(const btVector3& dir)
{
    btAssert(btFabs(dir.length()-1.0f) < 0.001f);
    // select which face of the cube we should map to
    int iFace = dir.closestAxis();
    if (dir[iFace] < 0.0f)
    {
        iFace += 3;
    }
    btAssert(iFace < 6);
    const FaceMap& face = gCubeFaceMaps[iFace];
    btVector3 n = face.m_normal;
    // project point on unit sphere onto the cube face
    float tt = 1.0f / n.dot(dir);
    btVector3 ptOnFace = dir * tt;
    float u = face.m_u + 0.5*btClamped(ptOnFace.dot(face.m_uVec), btScalar(-1), btScalar(1));
    float v = face.m_v + 0.5*btClamped(ptOnFace.dot(face.m_vVec), btScalar(-1), btScalar(1));
    btAssert(u >= 0 && u <= 1.0);
    btAssert(v >= 0 && v < 6.0);
    return btVector3(u, v, 0.0f);
}


static btVector3 cubeUvToNormalizedDirection(float u, float v)
{
    btAssert(u >= 0 && u <= 1.0);
    btAssert(v >= 0 && v < 6.0);
    int iFace = int(floor(v));
    v -= float(iFace);
    btAssert(iFace >= 0 && iFace < 6);
    const FaceMap& face = gCubeFaceMaps[iFace];
    btVector3 dir = face.m_normal + face.m_uVec*(u-0.5)*2.0 + face.m_vVec*(v-0.5)*2.0;
    dir.normalize();
    return dir;
}


struct CubeMap
{
    static const int N = 32;
    char data[6*N][N];

    void buildCubeMapLookupFromAxes(const btVector3* axes, int numAxes)
    {
        float cubeStep = 1.0 / N;
        for ( int iy = 0; iy < N*6; ++iy )
        {
            float v = cubeStep*( 0.5f + iy );
            for ( int ix = 0; ix < N; ++ix )
            {
                float u = cubeStep*( 0.5f + ix );
                btVector3 dir = cubeUvToNormalizedDirection( u, v );
                btVector3 uv = normalizedDirectionToCubeUv( dir );
                btAssert( btFabs( uv.x() - u ) < cubeStep*1.5 );
                btAssert( btFabs( uv.y() - v ) < cubeStep*1.5 );
                // find closest axis
                int axis = 0;
                float bestDot = 0.0f;
                for ( int i = 0; i < numAxes; ++i )
                {
                    float dot = btFabs( dir.dot( axes[ i ] ) );
                    if ( dot > bestDot )
                    {
                        bestDot = dot;
                        axis = i;
                    }
                }
                data[ iy ][ ix ] = axis;
            }
        }
    }
};


struct AssignConstraintsToPhasesParams
{
    const CubeMap* cubeMap;
    btBatchedConstraintInfo* conInfos;
    const btSolverConstraint* constraints;
    const btVector3* bodyPositions;
    int numAxes;
    char* constraintPhaseIds;
    int* constraintBatchIds;

    AssignConstraintsToPhasesParams(char* _constraintPhaseIds,
        int* _constraintBatchIds,
        btBatchedConstraintInfo* _conInfos,
        const btSolverConstraint* _constraints,
        const CubeMap* _cubeMap,
        const btVector3* _bodyPositions,
        int _numAxes
    )
    {
        constraintPhaseIds = _constraintPhaseIds;
        constraintBatchIds = _constraintBatchIds;
        conInfos = _conInfos;
        constraints = _constraints;
        numAxes = _numAxes;
        cubeMap = _cubeMap;
        bodyPositions = _bodyPositions;
    }
};


static void assignConstraintsToPhases(const AssignConstraintsToPhasesParams& params, int begin, int end)
{
    BT_PROFILE("assignConstraintsToPhases");
    int cubeMapSize = CubeMap::N;
    int numPhases = params.numAxes*2;
    for (int iCon = begin; iCon < end; ++iCon)
    {
        btBatchedConstraintInfo& conInfo = params.conInfos[iCon];
        const btSolverConstraint& con = params.constraints[iCon];
        int body0 = con.m_solverBodyIdA;
        int body1 = con.m_solverBodyIdB;
        conInfo.bodyIds[0] = body0;
        conInfo.bodyIds[1] = body1;
        btVector3 v = params.bodyPositions[body1] - params.bodyPositions[body0];
        btVector3 vDir = v.normalized();

        btVector3 uv = normalizedDirectionToCubeUv( vDir );
        int ix = int( uv.x() * cubeMapSize );
        int iy = int( uv.y() * cubeMapSize );
        int axis = params.cubeMap->data[ iy ][ ix ];
        btAssert( axis >= 0 && axis < params.numAxes );

        int iPhase = axis*2 + ((body0)&1);
        btAssert(iPhase >= 0 && iPhase < numPhases);
        params.constraintPhaseIds[iCon] = iPhase;
        params.constraintBatchIds[iCon] = -1;
    }
}


class AssignConstraintsToPhasesLoop : public btIParallelForBody
{
    const AssignConstraintsToPhasesParams* m_params;
public:
    AssignConstraintsToPhasesLoop(const AssignConstraintsToPhasesParams& params)
    {
        m_params = &params;
    }

    void forLoop( int iBegin, int iEnd ) const BT_OVERRIDE
    {
        BT_PROFILE( "AssignConstraintsToPhasesLoop" );
        assignConstraintsToPhases(*m_params, iBegin, iEnd);
    }
};


static int makePhaseRemappingTable(char* phaseMappingTable, int* numConstraintsPerPhase, int numPhases)
{
    BT_PROFILE("makePhaseRemappingTable");
    // join together small phases to try to create fewer more evenly-sized phases
    for (int i = 0; i < numPhases; ++i)
    {
        phaseMappingTable[i] = i;
    }
    int maxConstraintsPerPhase = 0;
    for (int i = 0; i < numPhases; ++i)
    {
        maxConstraintsPerPhase = btMax(maxConstraintsPerPhase, numConstraintsPerPhase[i]);
    }
    int minAllowedConstraintsPerPhase = maxConstraintsPerPhase/4;
    for (int iSrc = 0; iSrc < numPhases; ++iSrc)
    {
        if (numConstraintsPerPhase[iSrc] > 0 && numConstraintsPerPhase[iSrc] < minAllowedConstraintsPerPhase)
        {
            // find dest phase to merge into
            // find smallest valid phase
            int bestCount = INT_MAX;
            int iDest = 0;
            for ( int i = 0; i < numPhases; ++i )
            {
                int n = numConstraintsPerPhase[i];
                if (n >= minAllowedConstraintsPerPhase && n < bestCount)
                {
                    bestCount = n;
                    iDest = i;
                }
            }
            // merge phases
            numConstraintsPerPhase[iDest] += numConstraintsPerPhase[iSrc];
            numConstraintsPerPhase[iSrc] = 0;
            phaseMappingTable[iSrc] = iDest;
        }
    }
    for (int iDest = 0; iDest < numPhases; ++iDest)
    {
        // if slot is empty,
        if (numConstraintsPerPhase[iDest] == 0)
        {
            // find next non-empty slot,
            int iSrc = iDest;
            for ( int i = iDest + 1; i < numPhases; ++i )
            {
                if ( numConstraintsPerPhase[ i ] > 0 )
                {
                    iSrc = i;
                    break;
                }
            }
            // if found
            if (iSrc != iDest)
            {
                for ( int i = 0; i < numPhases; ++i )
                {
                    if (phaseMappingTable[i] == iSrc)
                    {
                        phaseMappingTable[i] = iDest;
                    }
                }
                btSwap(numConstraintsPerPhase[ iSrc ], numConstraintsPerPhase[ iDest ]);
            }
        }
    }
    int numActualPhases = 0;
    for ( int i = 0; i < numPhases; ++i )
    {
        if (numConstraintsPerPhase[ i ] > 0)
        {
            numActualPhases++;
        }
    }
    return numActualPhases;
}


template <int N>
class PreallocatedMemoryHelper
{
    struct Chunk
    {
        void** ptr;
        size_t size;
    };
    Chunk m_chunks[N];
    int m_numChunks;
public:
    PreallocatedMemoryHelper() {m_numChunks=0;}
    void addChunk( void** ptr, size_t sz )
    {
        btAssert( m_numChunks < N );
        if ( m_numChunks < N )
        {
            Chunk& chunk = m_chunks[ m_numChunks ];
            chunk.ptr = ptr;
            chunk.size = sz;
            m_numChunks++;
        }
    }
    size_t getSizeToAllocate() const
    {
        size_t totalSize = 0;
        for (int i = 0; i < m_numChunks; ++i)
        {
            totalSize += m_chunks[i].size;
        }
        return totalSize;
    }
    void setChunkPointers(void* mem) const
    {
        size_t totalSize = 0;
        for (int i = 0; i < m_numChunks; ++i)
        {
            const Chunk& chunk = m_chunks[ i ];
            char* chunkPtr = static_cast<char*>(mem) + totalSize;
            *chunk.ptr = chunkPtr;
            totalSize += chunk.size;
        }
    }
};


//
// setupMt2 -- we want a way to ensure that a single batch does not get too large
//
// KEY IDEA: if we assign constraints to phases based on the DIRECTION vector between
//    the 2 bodies of the constraint, then that would tend to limit how large a
//    batch within a phase could get since it could only string together a clump of
//    bodies along a particular axis until it hits the ends of the island of bodies.
//    Imagine tightly packed boxes on an axis-aligned 3d grid. There would be only
//    3 phases in the initial binning: constraints in the up/down direction, constraints
//    north/south and constraints east/west.
//
//    We'd need to quantize the direction vector somehow, and bin the constraints.
//    It may be useful to compute the bounding box of the body positions.
//
//    We pre-generate a small cubemap, which we can use to map from the direction
//    vector between bodies.
//    Imagine a cube centered at the origin.
//    For each of the 6 faces, we have a direction from the origin to the center of
//    each face.
//    For each of the 12 edges, we can define a direction from the origin to the
//    center of each edge.
//    For each of the 8 corners, we can define a direction from the origin to the corner.
//
//    That gives 26 directions in all, but only half as many axes. Each axis passes
//    through 2 faces, or 2 edges, or 2 corners. So there are 13 axes defined.
//    The choice of axes is fairly arbitrary, other choices are possible.
//    For each axis, we pre-allocate 2 phases. Once we map a constraint to one of
//    the 13 axes, we assign it to one of the 2 phases in a way that leads to an
//    even distribution of constraints.
//
//    Once all of the constraints are assigned to phases, we even out the size of
//    phases by merging together the 2 smallest phases in the group until the smallest
//    phase is at least one-fourth as large as the largest phase.
//
//    Once all constraints have been assigned to a phase we start generating batches for
//    each phase. There is no guaranttee that we will get a decent number of batches
//    for every phase. We could end up with a phase that has all its constraints in a
//    single batch. However, this doesn't seem to be an issue with testing so far.
//
static void setupMt2(
    BatchedConstraints* batchedConstraints,
    btAlignedObjectArray<char>* scratchMemory,
    btConstraintArray* constraints,
    const btAlignedObjectArray<btSolverBody>& bodies,
    int minBatchSize,
    int maxBatchSize,
    CreateBatchesWork* workArray
    )
{
    BT_PROFILE("setupMt2");
    btVector3 axes[] =
    {
        // cube face axes
        btVector3(1,0,0),
        btVector3(0,1,0),
        btVector3(0,0,1),
        // cube edge axes
        btVector3(1,0,1),
        btVector3(-1,0,1),
        btVector3(0,1,1),
        btVector3(0,-1,1),
        btVector3(1,1,0),
        btVector3(-1,1,0),
        // cube corner axes
        btVector3(1,1,1),
        btVector3(-1,1,1),
        btVector3(1,-1,1),
        btVector3(-1,-1,1),
    };
    const int numAxes = sizeof(axes)/sizeof(axes[0]);
    const bool useCubeMap = true;
    static CubeMap cubeMap;
    const int cubeMapSize = CubeMap::N;
    static bool builtCubeMap = false;
    if (!builtCubeMap)
    {
        for ( int i = 0; i < numAxes; ++i )
        {
            axes[ i ].normalize();
        }

        cubeMap.buildCubeMapLookupFromAxes(axes, numAxes);
        builtCubeMap = true;
    }

    const int numPhases = numAxes*2;
    int numConstraints = constraints->size();

    int maxNumBatchesPerPhase = bodies.size();
    int minNumBatchesPerPhase = 16;
    int maxNumBatches = maxNumBatchesPerPhase * numPhases;

    btVector3* bodyPositions = NULL;
    bool* bodyDynamicFlags = NULL;
    btBatchInfo* batches = NULL;
    int* batchWork = NULL;
    btBatchedConstraintInfo* conInfos = NULL;
    char* constraintPhaseIds = NULL;
    int* constraintBatchIds = NULL;
    {
        PreallocatedMemoryHelper<7> memHelper;
        memHelper.addChunk( (void**) &bodyPositions, sizeof( btVector3 ) * bodies.size() );
        memHelper.addChunk( (void**) &bodyDynamicFlags, sizeof( bool ) * bodies.size() );
        memHelper.addChunk( (void**) &batches, sizeof( btBatchInfo )* maxNumBatches );
        memHelper.addChunk( (void**) &batchWork, sizeof( int )* maxNumBatches );
        memHelper.addChunk( (void**) &conInfos, sizeof( btBatchedConstraintInfo ) * numConstraints );
        memHelper.addChunk( (void**) &constraintPhaseIds, sizeof( char ) * numConstraints );
        memHelper.addChunk( (void**) &constraintBatchIds, sizeof( int ) * numConstraints );
        size_t scratchSize = memHelper.getSizeToAllocate();
        scratchMemory->resizeNoInitialize( scratchSize );
        char* memPtr = &scratchMemory->at(0);
        memHelper.setChunkPointers( memPtr );
    }
    //bodyDynamicFlags.resizeNoInitialize(bodies.size());

    for (int i = 0; i < bodies.size(); ++i)
    {
        const btSolverBody& body = bodies[i];
        bodyPositions[i] = body.getWorldTransform().getOrigin();
        bodyDynamicFlags[i] = ( body.internalGetInvMass().x() > btScalar( 0 ) );
    }

    //for (int i = 0; i < numConstraints; ++i)
    //{
    //    constraintBatchIds[i] = -1;
    //}
    //for (int i = 0; i < maxNumBatches; ++i)
    //{
    //    btBatchInfo& batch = batches[i];
    //    batch = btBatchInfo();
    //}

    int numConstraintsPerPhase[numPhases];
    for (int i = 0; i < numPhases; ++i)
    {
        numConstraintsPerPhase[i] = 0;
    }

    bool assignInParallel = true;
    {
        AssignConstraintsToPhasesParams params( &constraintPhaseIds[ 0 ],
            constraintBatchIds,
            &conInfos[ 0 ],
            &constraints->at(0),
            &cubeMap,
            &bodyPositions[ 0 ],
            numAxes
        );
        if ( assignInParallel )
        {
            AssignConstraintsToPhasesLoop loop( params );
            btParallelFor( 0, numConstraints, 512, loop );
        }
        else
        {
            assignConstraintsToPhases( params, 0, numConstraints );
        }
    }

    for (int iCon = 0; iCon < numConstraints; ++iCon)
    {
        int iPhase = constraintPhaseIds[iCon];
        numConstraintsPerPhase[iPhase]++;
    }
    char phaseMappingTable[numPhases];
    int numActualPhases = makePhaseRemappingTable(phaseMappingTable, numConstraintsPerPhase, numPhases);

    CreateBatchesParams params;
    params.constraintBatchIds = &constraintBatchIds[0];
    params.batches = &batches[0];
    params.bodyDynamicFlags = &bodyDynamicFlags[0];
    params.conInfos = &conInfos[0];
    params.groupPhaseTable = NULL;
    params.constraintPhaseIds = &constraintPhaseIds[0];
    params.phaseMappingTable = phaseMappingTable;
    params.numBodies = bodies.size();
    params.numConstraints = numConstraints;
    params.numGroups = numPhases;
    params.minNumBatchesPerPhase = minNumBatchesPerPhase;
    params.maxNumBatchesPerPhase = maxNumBatchesPerPhase;
    params.minBatchSize = minBatchSize;
    params.maxBatchSize = maxBatchSize;
    params.workArray = workArray;

    if (true)
    {
        // parallel batch creation (deterministic)
        CreateBatchesForPhaseLoop loop(params);
        btParallelFor( 0, numActualPhases, 1, loop );
    }
    else
    {
        // sequential batch creation
        for ( int iPhase = 0; iPhase < numActualPhases; ++iPhase )
        {
            createBatchesForPhase( iPhase, params );
        }
    }
    // all constraints have been assigned a batchId
    updateConstraintBatchIdsForMergesMt(constraintBatchIds, numConstraints, batches, maxNumBatches);

    writeOutBatches(batchedConstraints, constraintBatchIds, numConstraints, batches, batchWork, maxNumBatchesPerPhase, numActualPhases);
    btAssert(batchedConstraints->validate(constraints, bodies));
}


static void setupSingleBatch(
    BatchedConstraints* bc,
    int numConstraints
)
{
    BT_PROFILE("setupSingleBatch");
    typedef BatchedConstraints::Range Range;

    bc->m_constraintIndices.resize( numConstraints );
    for ( int i = 0; i < numConstraints; ++i )
    {
        bc->m_constraintIndices[ i ] = i;
    }

    bc->m_batches.resizeNoInitialize( 0 );
    bc->m_phases.resizeNoInitialize( 0 );
    bc->m_phaseOrder.resizeNoInitialize( 0 );

    if (numConstraints > 0)
    {
        bc->m_batches.push_back( Range( 0, numConstraints ) );
        bc->m_phases.push_back( Range( 0, 1 ) );
        bc->m_phaseOrder.push_back(0);
    }
}


void BatchedConstraints::setup(
    btConstraintArray* constraints,
    const btAlignedObjectArray<btSolverBody>& bodies,
    int minBatchSize,
    int maxBatchSize,
    btAlignedObjectArray<char>* scratchMemory,
    CreateBatchesWork* workArray
    )
{
    if (constraints->size() >= 200)
    {
        setupMt2( this, scratchMemory, constraints, bodies, minBatchSize, maxBatchSize, workArray );
        //setupBatchesGreedyWithMerging(this, constraints, bodies, minBatchSize, maxBatchSize);
    }
    else
    {
        setupSingleBatch( this, constraints->size() );
    }
}


btSequentialImpulseConstraintSolverMt::btSequentialImpulseConstraintSolverMt()
{
    m_numFrictionDirections = 1;
    m_useBatching = false;
    m_useObsoleteJointConstraints = false;
}


btSequentialImpulseConstraintSolverMt::~btSequentialImpulseConstraintSolverMt()
{
}


void btSequentialImpulseConstraintSolverMt::setupBatchedContactConstraints()
{
    BT_PROFILE("setupBatchedContactConstraints");
    m_batchedContactConstraints.setup( &m_tmpSolverContactConstraintPool, m_tmpSolverBodyPool, sMinBatchSize, sMaxBatchSize, &m_scratchMemory, m_createBatchesWorkArray );
}


void btSequentialImpulseConstraintSolverMt::setupBatchedJointConstraints()
{
    BT_PROFILE("setupBatchedJointConstraints");
    m_batchedJointConstraints.setup( &m_tmpSolverNonContactConstraintPool, m_tmpSolverBodyPool, sMinBatchSize, sMaxBatchSize, &m_scratchMemory, m_createBatchesWorkArray );
}


void btSequentialImpulseConstraintSolverMt::internalSetupContactConstraints(int iContactConstraint, const btContactSolverInfo& infoGlobal)
{
    btSolverConstraint& contactConstraint = m_tmpSolverContactConstraintPool[iContactConstraint];

    btVector3 rel_pos1;
    btVector3 rel_pos2;
    btScalar relaxation;

    int solverBodyIdA = contactConstraint.m_solverBodyIdA;
    int solverBodyIdB = contactConstraint.m_solverBodyIdB;

    btSolverBody* solverBodyA = &m_tmpSolverBodyPool[ solverBodyIdA ];
    btSolverBody* solverBodyB = &m_tmpSolverBodyPool[ solverBodyIdB ];

    btRigidBody* colObj0 = solverBodyA->m_originalBody;
    btRigidBody* colObj1 = solverBodyB->m_originalBody;

    btManifoldPoint& cp = *static_cast<btManifoldPoint*>( contactConstraint.m_originalContactPoint );

    const btVector3& pos1 = cp.getPositionWorldOnA();
    const btVector3& pos2 = cp.getPositionWorldOnB();

    rel_pos1 = pos1 - solverBodyA->getWorldTransform().getOrigin();
    rel_pos2 = pos2 - solverBodyB->getWorldTransform().getOrigin();

    btVector3 vel1;
    btVector3 vel2;

    solverBodyA->getVelocityInLocalPointNoDelta( rel_pos1, vel1 );
    solverBodyB->getVelocityInLocalPointNoDelta( rel_pos2, vel2 );

    btVector3 vel = vel1 - vel2;
    btScalar rel_vel = cp.m_normalWorldOnB.dot( vel );

    setupContactConstraint( contactConstraint, solverBodyIdA, solverBodyIdB, cp, infoGlobal, relaxation, rel_pos1, rel_pos2 );

    // setup rolling friction constraints
    int rollingFrictionIndex = m_rollingFrictionIndexTable[iContactConstraint];
    if (rollingFrictionIndex >= 0)
    {
        btSolverConstraint& spinningFrictionConstraint = m_tmpSolverContactRollingFrictionConstraintPool[ rollingFrictionIndex ];
        btAssert( spinningFrictionConstraint.m_frictionIndex == iContactConstraint );
        setupTorsionalFrictionConstraint( spinningFrictionConstraint,
            cp.m_normalWorldOnB,
            solverBodyIdA,
            solverBodyIdB,
            cp,
            cp.m_combinedSpinningFriction,
            rel_pos1,
            rel_pos2,
            colObj0,
            colObj1,
            relaxation,
            0.0f,
            0.0f
        );
        btVector3 axis[2];
        btPlaneSpace1( cp.m_normalWorldOnB, axis[0], axis[1] );
        axis[0].normalize();
        axis[1].normalize();

        applyAnisotropicFriction( colObj0, axis[0], btCollisionObject::CF_ANISOTROPIC_ROLLING_FRICTION );
        applyAnisotropicFriction( colObj1, axis[0], btCollisionObject::CF_ANISOTROPIC_ROLLING_FRICTION );
        applyAnisotropicFriction( colObj0, axis[1], btCollisionObject::CF_ANISOTROPIC_ROLLING_FRICTION );
        applyAnisotropicFriction( colObj1, axis[1], btCollisionObject::CF_ANISOTROPIC_ROLLING_FRICTION );
        // put the largest axis first
        if (axis[1].length2() > axis[0].length2())
        {
            btSwap(axis[0], axis[1]);
        }
        const btScalar kRollingFrictionThreshold = 0.001f;
        for (int i = 0; i < 2; ++i)
        {
            int iRollingFric = rollingFrictionIndex + 1 + i;
            btSolverConstraint& rollingFrictionConstraint = m_tmpSolverContactRollingFrictionConstraintPool[ iRollingFric ];
            btAssert(rollingFrictionConstraint.m_frictionIndex == iContactConstraint);
            btVector3 dir = axis[i];
            if ( dir.length() > kRollingFrictionThreshold )
            {
                setupTorsionalFrictionConstraint( rollingFrictionConstraint,
                    dir,
                    solverBodyIdA,
                    solverBodyIdB,
                    cp,
                    cp.m_combinedRollingFriction,
                    rel_pos1,
                    rel_pos2,
                    colObj0,
                    colObj1,
                    relaxation,
                    0.0f,
                    0.0f
                );
            }
            else
            {
                rollingFrictionConstraint.m_frictionIndex = -1;  // disable constraint
            }
        }
    }

    // setup friction constraints
    //	setupFrictionConstraint(solverConstraint, normalAxis, solverBodyIdA, solverBodyIdB, cp, rel_pos1, rel_pos2, colObj0, colObj1, relaxation, infoGlobal, desiredVelocity, cfmSlip);
    {
        ///Bullet has several options to set the friction directions
        ///By default, each contact has only a single friction direction that is recomputed automatically very frame
        ///based on the relative linear velocity.
        ///If the relative velocity it zero, it will automatically compute a friction direction.

        ///You can also enable two friction directions, using the SOLVER_USE_2_FRICTION_DIRECTIONS.
        ///In that case, the second friction direction will be orthogonal to both contact normal and first friction direction.
        ///
        ///If you choose SOLVER_DISABLE_VELOCITY_DEPENDENT_FRICTION_DIRECTION, then the friction will be independent from the relative projected velocity.
        ///
        ///The user can manually override the friction directions for certain contacts using a contact callback,
        ///and set the cp.m_lateralFrictionInitialized to true
        ///In that case, you can set the target relative motion in each friction direction (cp.m_contactMotion1 and cp.m_contactMotion2)
        ///this will give a conveyor belt effect
        ///
	    btSolverConstraint* frictionConstraint1 = &m_tmpSolverContactFrictionConstraintPool[contactConstraint.m_frictionIndex];
        btAssert(frictionConstraint1->m_frictionIndex == iContactConstraint);

        btSolverConstraint* frictionConstraint2 = NULL;
        if ( infoGlobal.m_solverMode & SOLVER_USE_2_FRICTION_DIRECTIONS )
        {
            frictionConstraint2 = &m_tmpSolverContactFrictionConstraintPool[contactConstraint.m_frictionIndex + 1];
            btAssert( frictionConstraint2->m_frictionIndex == iContactConstraint );
        }

        if ( !( infoGlobal.m_solverMode & SOLVER_ENABLE_FRICTION_DIRECTION_CACHING ) || !( cp.m_contactPointFlags&BT_CONTACT_FLAG_LATERAL_FRICTION_INITIALIZED ) )
        {
            cp.m_lateralFrictionDir1 = vel - cp.m_normalWorldOnB * rel_vel;
            btScalar lat_rel_vel = cp.m_lateralFrictionDir1.length2();
            if ( !( infoGlobal.m_solverMode & SOLVER_DISABLE_VELOCITY_DEPENDENT_FRICTION_DIRECTION ) && lat_rel_vel > SIMD_EPSILON )
            {
                cp.m_lateralFrictionDir1 *= 1.f / btSqrt( lat_rel_vel );
                applyAnisotropicFriction( colObj0, cp.m_lateralFrictionDir1, btCollisionObject::CF_ANISOTROPIC_FRICTION );
                applyAnisotropicFriction( colObj1, cp.m_lateralFrictionDir1, btCollisionObject::CF_ANISOTROPIC_FRICTION );
                setupFrictionConstraint( *frictionConstraint1, cp.m_lateralFrictionDir1, solverBodyIdA, solverBodyIdB, cp, rel_pos1, rel_pos2, colObj0, colObj1, relaxation, infoGlobal );

                if ( frictionConstraint2 )
                {
                    cp.m_lateralFrictionDir2 = cp.m_lateralFrictionDir1.cross( cp.m_normalWorldOnB );
                    cp.m_lateralFrictionDir2.normalize();//??
                    applyAnisotropicFriction( colObj0, cp.m_lateralFrictionDir2, btCollisionObject::CF_ANISOTROPIC_FRICTION );
                    applyAnisotropicFriction( colObj1, cp.m_lateralFrictionDir2, btCollisionObject::CF_ANISOTROPIC_FRICTION );
                    setupFrictionConstraint( *frictionConstraint2, cp.m_lateralFrictionDir2, solverBodyIdA, solverBodyIdB, cp, rel_pos1, rel_pos2, colObj0, colObj1, relaxation, infoGlobal );
                }
            }
            else
            {
                btPlaneSpace1( cp.m_normalWorldOnB, cp.m_lateralFrictionDir1, cp.m_lateralFrictionDir2 );

                applyAnisotropicFriction( colObj0, cp.m_lateralFrictionDir1, btCollisionObject::CF_ANISOTROPIC_FRICTION );
                applyAnisotropicFriction( colObj1, cp.m_lateralFrictionDir1, btCollisionObject::CF_ANISOTROPIC_FRICTION );
                setupFrictionConstraint( *frictionConstraint1, cp.m_lateralFrictionDir1, solverBodyIdA, solverBodyIdB, cp, rel_pos1, rel_pos2, colObj0, colObj1, relaxation, infoGlobal );

                if ( frictionConstraint2 )
                {
                    applyAnisotropicFriction( colObj0, cp.m_lateralFrictionDir2, btCollisionObject::CF_ANISOTROPIC_FRICTION );
                    applyAnisotropicFriction( colObj1, cp.m_lateralFrictionDir2, btCollisionObject::CF_ANISOTROPIC_FRICTION );
                    setupFrictionConstraint( *frictionConstraint2, cp.m_lateralFrictionDir2, solverBodyIdA, solverBodyIdB, cp, rel_pos1, rel_pos2, colObj0, colObj1, relaxation, infoGlobal );
                }

                if ( ( infoGlobal.m_solverMode & SOLVER_USE_2_FRICTION_DIRECTIONS ) && ( infoGlobal.m_solverMode & SOLVER_DISABLE_VELOCITY_DEPENDENT_FRICTION_DIRECTION ) )
                {
                    cp.m_contactPointFlags |= BT_CONTACT_FLAG_LATERAL_FRICTION_INITIALIZED;
                }
            }
        }
        else
        {
            setupFrictionConstraint( *frictionConstraint1, cp.m_lateralFrictionDir1, solverBodyIdA, solverBodyIdB, cp, rel_pos1, rel_pos2, colObj0, colObj1, relaxation, infoGlobal, cp.m_contactMotion1, cp.m_frictionCFM );
            if ( frictionConstraint2 )
            {
                setupFrictionConstraint( *frictionConstraint2, cp.m_lateralFrictionDir2, solverBodyIdA, solverBodyIdB, cp, rel_pos1, rel_pos2, colObj0, colObj1, relaxation, infoGlobal, cp.m_contactMotion2, cp.m_frictionCFM );
            }
        }
    }

    setFrictionConstraintImpulse( contactConstraint, solverBodyIdA, solverBodyIdB, cp, infoGlobal );
}


struct SetupContactConstraintsLoop : public btIParallelForBody
{
    btSequentialImpulseConstraintSolverMt* m_solver;
    const BatchedConstraints* m_bc;
    const btContactSolverInfo* m_infoGlobal;

    SetupContactConstraintsLoop( btSequentialImpulseConstraintSolverMt* solver, const BatchedConstraints* bc, const btContactSolverInfo& infoGlobal )
    {
        m_solver = solver;
        m_bc = bc;
        m_infoGlobal = &infoGlobal;
    }
    void forLoop( int iBegin, int iEnd ) const BT_OVERRIDE
    {
        BT_PROFILE( "SetupContactConstraintsLoop" );
        for ( int iBatch = iBegin; iBatch < iEnd; ++iBatch )
        {
            const BatchedConstraints::Range& batch = m_bc->m_batches[ iBatch ];
            for (int i = batch.begin; i < batch.end; ++i)
            {
                int iContact = m_bc->m_constraintIndices[i];
                m_solver->internalSetupContactConstraints( iContact, *m_infoGlobal );
            }
        }
    }
};


void btSequentialImpulseConstraintSolverMt::setupAllContactConstraints(const btContactSolverInfo& infoGlobal)
{
    BT_PROFILE( "setupAllContactConstraints" );
    if ( m_useBatching )
    {
        const BatchedConstraints& batchedCons = m_batchedContactConstraints;
        SetupContactConstraintsLoop loop( this, &batchedCons, infoGlobal );
        for ( int iiPhase = 0; iiPhase < batchedCons.m_phases.size(); ++iiPhase )
        {
            int iPhase = batchedCons.m_phaseOrder[ iiPhase ];
            const BatchedConstraints::Range& phase = batchedCons.m_phases[ iPhase ];
            int grainSize = 1;
            btParallelFor( phase.begin, phase.end, grainSize, loop );
        }
    }
    else
    {
        for ( int i = 0; i < m_tmpSolverContactConstraintPool.size(); ++i )
        {
            internalSetupContactConstraints( i, infoGlobal );
        }
    }
}


int	btSequentialImpulseConstraintSolverMt::getOrInitSolverBodyThreadsafe(btCollisionObject& body,btScalar timeStep)
{
    //
    // getOrInitSolverBody is threadsafe only for a single thread per solver (with potentially multiple solvers)
    //
    // getOrInitSolverBodyThreadsafe -- attempts to be fully threadsafe (however may affect determinism)
    //
    int solverBodyId = -1;
    if ( !body.isStaticOrKinematicObject() )
    {
        // dynamic body
        // Dynamic bodies can only be in one island, so it's safe to write to the companionId
        solverBodyId = body.getCompanionId();
        if ( solverBodyId < 0 )
        {
            m_bodySolverArrayMutex.lock();
            // now that we have the lock, check again
            solverBodyId = body.getCompanionId();
            if ( solverBodyId < 0 )
            {
                if ( btRigidBody* rb = btRigidBody::upcast( &body ) )
                {
                    solverBodyId = m_tmpSolverBodyPool.size();
                    btSolverBody& solverBody = m_tmpSolverBodyPool.expand();
                    initSolverBody( &solverBody, &body, timeStep );
                    body.setCompanionId( solverBodyId );
                }
            }
            m_bodySolverArrayMutex.unlock();
        }
    }
    else if (body.isKinematicObject())
    {
        //
        // NOTE: must test for kinematic before static because some kinematic objects also
        //   identify as "static"
        //
        // Kinematic bodies can be in multiple islands at once, so it is a
        // race condition to write to them, so we use an alternate method
        // to record the solverBodyId
        int uniqueId = body.getWorldArrayIndex();
        const int INVALID_SOLVER_BODY_ID = -1;
        if (m_kinematicBodyUniqueIdToSolverBodyTable.size() <= uniqueId )
        {
            m_kinematicBodyUniqueIdToSolverBodyTableMutex.lock();
            // now that we have the lock, check again
            if ( m_kinematicBodyUniqueIdToSolverBodyTable.size() <= uniqueId )
            {
                m_kinematicBodyUniqueIdToSolverBodyTable.resize( uniqueId + 1, INVALID_SOLVER_BODY_ID );
            }
            m_kinematicBodyUniqueIdToSolverBodyTableMutex.unlock();
        }
        solverBodyId = m_kinematicBodyUniqueIdToSolverBodyTable[ uniqueId ];
        // if no table entry yet,
        if ( INVALID_SOLVER_BODY_ID == solverBodyId )
        {
            // need to acquire both locks
            m_kinematicBodyUniqueIdToSolverBodyTableMutex.lock();
            m_bodySolverArrayMutex.lock();
            // now that we have the lock, check again
            solverBodyId = m_kinematicBodyUniqueIdToSolverBodyTable[ uniqueId ];
            if ( INVALID_SOLVER_BODY_ID == solverBodyId )
            {
                // create a table entry for this body
                btRigidBody* rb = btRigidBody::upcast( &body );
                solverBodyId = m_tmpSolverBodyPool.size();
                btSolverBody& solverBody = m_tmpSolverBodyPool.expand();
                initSolverBody( &solverBody, &body, timeStep );
                m_kinematicBodyUniqueIdToSolverBodyTable[ uniqueId ] = solverBodyId;
            }
            m_bodySolverArrayMutex.unlock();
            m_kinematicBodyUniqueIdToSolverBodyTableMutex.unlock();
        }
    }
    else
    {
        // all fixed bodies (inf mass) get mapped to a single solver id
        if ( m_fixedBodyId < 0 )
        {
            m_bodySolverArrayMutex.lock();
            // now that we have the lock, check again
            if ( m_fixedBodyId < 0 )
            {
                m_fixedBodyId = m_tmpSolverBodyPool.size();
                btSolverBody& fixedBody = m_tmpSolverBodyPool.expand();
                initSolverBody( &fixedBody, 0, timeStep );
            }
            m_bodySolverArrayMutex.unlock();
        }
        solverBodyId = m_fixedBodyId;
    }
    btAssert( solverBodyId < m_tmpSolverBodyPool.size() );
	return solverBodyId;
}


void btSequentialImpulseConstraintSolverMt::internalCollectContactManifoldCachedInfo(btContactManifoldCachedInfo* cachedInfoArray, btPersistentManifold** manifoldPtr, int numManifolds, const btContactSolverInfo& infoGlobal)
{
    BT_PROFILE("internalCollectContactManifoldCachedInfo");
    for (int i = 0; i < numManifolds; ++i)
    {
        btContactManifoldCachedInfo* cachedInfo = &cachedInfoArray[i];
        btPersistentManifold* manifold = manifoldPtr[i];
        btCollisionObject* colObj0 = (btCollisionObject*) manifold->getBody0();
        btCollisionObject* colObj1 = (btCollisionObject*) manifold->getBody1();

        int solverBodyIdA = getOrInitSolverBodyThreadsafe( *colObj0, infoGlobal.m_timeStep );
        int solverBodyIdB = getOrInitSolverBodyThreadsafe( *colObj1, infoGlobal.m_timeStep );

        cachedInfo->solverBodyIds[ 0 ] = solverBodyIdA;
        cachedInfo->solverBodyIds[ 1 ] = solverBodyIdB;
        cachedInfo->numTouchingContacts = 0;

        btSolverBody* solverBodyA = &m_tmpSolverBodyPool[ solverBodyIdA ];
        btSolverBody* solverBodyB = &m_tmpSolverBodyPool[ solverBodyIdB ];

        ///avoid collision response between two static objects
        if ( solverBodyA->m_invMass.fuzzyZero() && solverBodyB->m_invMass.fuzzyZero() )
            break;

        int iContact = 0;
        for ( int j = 0; j < manifold->getNumContacts(); j++ )
        {
            btManifoldPoint& cp = manifold->getContactPoint( j );

            if ( cp.getDistance() <= manifold->getContactProcessingThreshold() )
            {
                cachedInfo->contactPoints[ iContact ] = &cp;
                cachedInfo->contactHasRollingFriction[ iContact ] = ( cp.m_combinedRollingFriction > 0.f );
                iContact++;
            }
        }
        cachedInfo->numTouchingContacts = iContact;
    }
}


struct CollectContactManifoldCachedInfoLoop : public btIParallelForBody
{
    btSequentialImpulseConstraintSolverMt* m_solver;
    btSequentialImpulseConstraintSolverMt::btContactManifoldCachedInfo* m_cachedInfoArray;
    btPersistentManifold** m_manifoldPtr;
    const btContactSolverInfo* m_infoGlobal;

    CollectContactManifoldCachedInfoLoop( btSequentialImpulseConstraintSolverMt* solver, btSequentialImpulseConstraintSolverMt::btContactManifoldCachedInfo* cachedInfoArray, btPersistentManifold** manifoldPtr, const btContactSolverInfo& infoGlobal )
    {
        m_solver = solver;
        m_cachedInfoArray = cachedInfoArray;
        m_manifoldPtr = manifoldPtr;
        m_infoGlobal = &infoGlobal;
    }
    void forLoop( int iBegin, int iEnd ) const BT_OVERRIDE
    {
        m_solver->internalCollectContactManifoldCachedInfo( m_cachedInfoArray + iBegin, m_manifoldPtr + iBegin, iEnd - iBegin, *m_infoGlobal );
    }
};


void btSequentialImpulseConstraintSolverMt::internalAllocContactConstraints(const btContactManifoldCachedInfo* cachedInfoArray, int numManifolds)
{
    BT_PROFILE("internalAllocContactConstraints");
    // possibly parallel part
    for ( int iManifold = 0; iManifold < numManifolds; ++iManifold )
    {
        const btContactManifoldCachedInfo& cachedInfo = cachedInfoArray[ iManifold ];
        int contactIndex = cachedInfo.contactIndex;
        int frictionIndex = contactIndex * m_numFrictionDirections;
        int rollingFrictionIndex = cachedInfo.rollingFrictionIndex;
        for ( int i = 0; i < cachedInfo.numTouchingContacts; i++ )
        {
            btSolverConstraint& contactConstraint = m_tmpSolverContactConstraintPool[contactIndex];
            contactConstraint.m_solverBodyIdA = cachedInfo.solverBodyIds[ 0 ];
            contactConstraint.m_solverBodyIdB = cachedInfo.solverBodyIds[ 1 ];
            contactConstraint.m_originalContactPoint = cachedInfo.contactPoints[ i ];

            // allocate the friction constraints
            contactConstraint.m_frictionIndex = frictionIndex;
            for ( int i = 0; i < m_numFrictionDirections; ++i )
            {
                btSolverConstraint& frictionConstraint = m_tmpSolverContactFrictionConstraintPool[frictionIndex];
                frictionConstraint.m_frictionIndex = contactIndex;
                frictionIndex++;
            }

            // allocate rolling friction constraints
            if ( cachedInfo.contactHasRollingFriction[ i ] )
            {
                m_rollingFrictionIndexTable[ contactIndex ] = rollingFrictionIndex;
                // allocate 3 (although we may use only 2 sometimes)
                for ( int i = 0; i < 3; i++ )
                {
                    m_tmpSolverContactRollingFrictionConstraintPool[ rollingFrictionIndex ].m_frictionIndex = contactIndex;
                    rollingFrictionIndex++;
                }
            }
            else
            {
                // indicate there is no rolling friction for this contact point
                m_rollingFrictionIndexTable[ contactIndex ] = -1;
            }
            contactIndex++;
        }
    }
}


struct AllocContactConstraintsLoop : public btIParallelForBody
{
    btSequentialImpulseConstraintSolverMt* m_solver;
    const btSequentialImpulseConstraintSolverMt::btContactManifoldCachedInfo* m_cachedInfoArray;

    AllocContactConstraintsLoop( btSequentialImpulseConstraintSolverMt* solver, btSequentialImpulseConstraintSolverMt::btContactManifoldCachedInfo* cachedInfoArray )
    {
        m_solver = solver;
        m_cachedInfoArray = cachedInfoArray;
    }
    void forLoop( int iBegin, int iEnd ) const BT_OVERRIDE
    {
        m_solver->internalAllocContactConstraints( m_cachedInfoArray + iBegin, iEnd - iBegin );
    }
};


void btSequentialImpulseConstraintSolverMt::allocAllContactConstraints(btPersistentManifold** manifoldPtr, int numManifolds, const btContactSolverInfo& infoGlobal)
{
    BT_PROFILE( "allocAllContactConstraints" );
    btAlignedObjectArray<btContactManifoldCachedInfo> cachedInfoArray; // = m_manifoldCachedInfoArray;
    cachedInfoArray.resizeNoInitialize( numManifolds );
    if (true)
    {
        // sequential
        internalCollectContactManifoldCachedInfo(&cachedInfoArray[ 0 ], manifoldPtr, numManifolds, infoGlobal);
    }
    else
    {
        // may alter ordering of bodies which affects determinism
        CollectContactManifoldCachedInfoLoop loop( this, &cachedInfoArray[ 0 ], manifoldPtr, infoGlobal );
        int grainSize = 200;
        btParallelFor( 0, numManifolds, grainSize, loop );
    }

    {
        // serial part
        int numContacts = 0;
        int numRollingFrictionConstraints = 0;
        for ( int iManifold = 0; iManifold < numManifolds; ++iManifold )
        {
            btContactManifoldCachedInfo& cachedInfo = cachedInfoArray[ iManifold ];
            cachedInfo.contactIndex = numContacts;
            cachedInfo.rollingFrictionIndex = numRollingFrictionConstraints;
            numContacts += cachedInfo.numTouchingContacts;
            for (int i = 0; i < cachedInfo.numTouchingContacts; ++i)
            {
                if (cachedInfo.contactHasRollingFriction[i])
                {
                    numRollingFrictionConstraints += 3;
                }
            }
        }
        m_tmpSolverContactConstraintPool.resizeNoInitialize(numContacts);
        m_rollingFrictionIndexTable.resizeNoInitialize(numContacts);
        m_tmpSolverContactFrictionConstraintPool.resizeNoInitialize(numContacts*m_numFrictionDirections);
        m_tmpSolverContactRollingFrictionConstraintPool.resizeNoInitialize(numRollingFrictionConstraints);
    }
    {
        AllocContactConstraintsLoop loop(this, &cachedInfoArray[0]);
        int grainSize = 200;
        btParallelFor( 0, numManifolds, grainSize, loop );
    }
}


void btSequentialImpulseConstraintSolverMt::convertContacts(btPersistentManifold** manifoldPtr, int numManifolds, const btContactSolverInfo& infoGlobal)
{
    if (!m_useBatching)
    {
        btSequentialImpulseConstraintSolver::convertContacts(manifoldPtr, numManifolds, infoGlobal);
        return;
    }
    BT_PROFILE( "convertContacts" );
    if (numManifolds > 0)
    {
        if ( m_fixedBodyId < 0 )
        {
            m_fixedBodyId = m_tmpSolverBodyPool.size();
            btSolverBody& fixedBody = m_tmpSolverBodyPool.expand();
            initSolverBody( &fixedBody, 0, infoGlobal.m_timeStep );
        }
        allocAllContactConstraints( manifoldPtr, numManifolds, infoGlobal );
        if ( m_useBatching )
        {
            setupBatchedContactConstraints();
        }
        setupAllContactConstraints( infoGlobal );
    }
}


btScalar btSequentialImpulseConstraintSolverMt::solveGroupCacheFriendlySetup(
     btCollisionObject** bodies,
     int numBodies,
     btPersistentManifold** manifoldPtr,
     int numManifolds,
     btTypedConstraint** constraints,
     int numConstraints,
     const btContactSolverInfo& infoGlobal,
     btIDebugDraw* debugDrawer
     )
{
    m_numFrictionDirections = (infoGlobal.m_solverMode & SOLVER_USE_2_FRICTION_DIRECTIONS) ? 2 : 1;
    m_useBatching = false;
    if ( numManifolds >= sMinimumContactManifoldsForBatching &&
        (sAllowNestedParallelForLoops || !btThreadsAreRunning())
        )
    {
        m_useBatching = true;
    }
    btSequentialImpulseConstraintSolver::solveGroupCacheFriendlySetup( bodies,
                                                                       numBodies,
                                                                       manifoldPtr,
                                                                       numManifolds,
                                                                       constraints,
                                                                       numConstraints,
                                                                       infoGlobal,
                                                                       debugDrawer
                                                                       );
    if ( m_useBatching )
    {
        setupBatchedJointConstraints();
    }
    return 0.0f;
}


btScalar btSequentialImpulseConstraintSolverMt::resolveMultipleContactSplitPenetrationImpulseConstraints( const btAlignedObjectArray<int>& consIndices, int batchBegin, int batchEnd )
{
    btScalar leastSquaresResidual = 0.f;
    for ( int iiCons = batchBegin; iiCons < batchEnd; ++iiCons )
    {
        int iCons = consIndices[ iiCons ];
        const btSolverConstraint& solveManifold = m_tmpSolverContactConstraintPool[ iCons ];
        btSolverBody& bodyA = m_tmpSolverBodyPool[ solveManifold.m_solverBodyIdA ];
        btSolverBody& bodyB = m_tmpSolverBodyPool[ solveManifold.m_solverBodyIdB ];
        btScalar residual = resolveSplitPenetrationImpulse( bodyA, bodyB, solveManifold );
        leastSquaresResidual += residual*residual;
    }
    return leastSquaresResidual;
}


struct ContactSplitPenetrationImpulseSolverLoop : public btIParallelSumBody
{
    btSequentialImpulseConstraintSolverMt* m_solver;
    const BatchedConstraints* m_bc;

    ContactSplitPenetrationImpulseSolverLoop( btSequentialImpulseConstraintSolverMt* solver, const BatchedConstraints* bc )
    {
        m_solver = solver;
        m_bc = bc;
    }
    btScalar sumLoop( int iBegin, int iEnd ) const BT_OVERRIDE
    {
        BT_PROFILE( "ContactSplitPenetrationImpulseSolverLoop" );
        btScalar sum = 0;
        for ( int iBatch = iBegin; iBatch < iEnd; ++iBatch )
        {
            const BatchedConstraints::Range& batch = m_bc->m_batches[ iBatch ];
            sum += m_solver->resolveMultipleContactSplitPenetrationImpulseConstraints( m_bc->m_constraintIndices, batch.begin, batch.end );
        }
        return sum;
    }
};


void btSequentialImpulseConstraintSolverMt::solveGroupCacheFriendlySplitImpulseIterations(btCollisionObject** bodies,int numBodies,btPersistentManifold** manifoldPtr, int numManifolds,btTypedConstraint** constraints,int numConstraints,const btContactSolverInfo& infoGlobal,btIDebugDraw* debugDrawer)
{
	BT_PROFILE("solveGroupCacheFriendlySplitImpulseIterations");
	if (infoGlobal.m_splitImpulse)
	{
        for ( int iteration = 0; iteration < infoGlobal.m_numIterations; iteration++ )
        {
            btScalar leastSquaresResidual = 0.f;
            if (m_useBatching)
            {
                const BatchedConstraints& batchedCons = m_batchedContactConstraints;
                ContactSplitPenetrationImpulseSolverLoop loop( this, &batchedCons );
                btScalar leastSquaresResidual = 0.f;
                for ( int iiPhase = 0; iiPhase < batchedCons.m_phases.size(); ++iiPhase )
                {
                    int iPhase = batchedCons.m_phaseOrder[ iiPhase ];
                    const BatchedConstraints::Range& phase = batchedCons.m_phases[ iPhase ];
                    int grainSize = 8;
                    leastSquaresResidual += btParallelSum( phase.begin, phase.end, grainSize, loop );
                }
            }
            else
            {
                // non-batched
                leastSquaresResidual = resolveMultipleContactSplitPenetrationImpulseConstraints(m_orderTmpConstraintPool, 0, m_tmpSolverContactConstraintPool.size());
            }
            if ( leastSquaresResidual <= infoGlobal.m_leastSquaresResidualThreshold || iteration >= ( infoGlobal.m_numIterations - 1 ) )
            {
#ifdef VERBOSE_RESIDUAL_PRINTF
                printf( "residual = %f at iteration #%d\n", leastSquaresResidual, iteration );
#endif
                break;
            }
        }
	}
}


btScalar btSequentialImpulseConstraintSolverMt::solveSingleIteration(int iteration, btCollisionObject** bodies, int numBodies, btPersistentManifold** manifoldPtr, int numManifolds, btTypedConstraint** constraints,int numConstraints,const btContactSolverInfo& infoGlobal,btIDebugDraw* debugDrawer)
{
    if ( !m_useBatching )
    {
        return btSequentialImpulseConstraintSolver::solveSingleIteration( iteration, bodies, numBodies, manifoldPtr, numManifolds, constraints, numConstraints, infoGlobal, debugDrawer );
    }
    BT_PROFILE( "solveSingleIterationMt" );
    btScalar leastSquaresResidual = 0.f;

	if (infoGlobal.m_solverMode & SOLVER_RANDMIZE_ORDER)
	{
		if (1)			// uncomment this for a bit less random ((iteration & 7) == 0)
		{
            randomizeConstraintOrdering(iteration, infoGlobal.m_numIterations);
		}
	}

	{
		///solve all joint constraints
        leastSquaresResidual += resolveAllJointConstraints(iteration);

		if (iteration< infoGlobal.m_numIterations)
		{
            // this loop is only used for cone-twist constraints,
            // it would be nice to skip this loop if none of the constraints need it
            if ( m_useObsoleteJointConstraints )
            {
                for ( int j = 0; j<numConstraints; j++ )
                {
                    if ( constraints[ j ]->isEnabled() )
                    {
                        int bodyAid = getOrInitSolverBody( constraints[ j ]->getRigidBodyA(), infoGlobal.m_timeStep );
                        int bodyBid = getOrInitSolverBody( constraints[ j ]->getRigidBodyB(), infoGlobal.m_timeStep );
                        btSolverBody& bodyA = m_tmpSolverBodyPool[ bodyAid ];
                        btSolverBody& bodyB = m_tmpSolverBodyPool[ bodyBid ];
                        constraints[ j ]->solveConstraintObsolete( bodyA, bodyB, infoGlobal.m_timeStep );
                    }
                }
            }

			if (infoGlobal.m_solverMode & SOLVER_INTERLEAVE_CONTACT_AND_FRICTION_CONSTRAINTS)
			{
                // solve all contact, contact-friction, and rolling friction constraints interleaved
                leastSquaresResidual += resolveAllContactConstraintsInterleaved();
			}
			else//SOLVER_INTERLEAVE_CONTACT_AND_FRICTION_CONSTRAINTS
			{
                // don't interleave them
				// solve all contact constraints
                leastSquaresResidual += resolveAllContactConstraints();

				// solve all contact friction constraints
                leastSquaresResidual += resolveAllContactFrictionConstraints();

                // solve all rolling friction constraints
                leastSquaresResidual += resolveAllRollingFrictionConstraints();
			}
		}
	}
    return leastSquaresResidual;
}


btScalar btSequentialImpulseConstraintSolverMt::resolveMultipleJointConstraints( const btAlignedObjectArray<int>& consIndices, int batchBegin, int batchEnd, int iteration )
{
    btScalar leastSquaresResidual = 0.f;
    for ( int iiCons = batchBegin; iiCons < batchEnd; ++iiCons )
    {
        int iCons = consIndices[ iiCons ];
        const btSolverConstraint& constraint = m_tmpSolverNonContactConstraintPool[ iCons ];
        if ( iteration < constraint.m_overrideNumSolverIterations )
        {
            btSolverBody& bodyA = m_tmpSolverBodyPool[ constraint.m_solverBodyIdA ];
            btSolverBody& bodyB = m_tmpSolverBodyPool[ constraint.m_solverBodyIdB ];
            btScalar residual = resolveSingleConstraintRowGeneric( bodyA, bodyB, constraint );
            leastSquaresResidual += residual*residual;
        }
    }
    return leastSquaresResidual;
}


btScalar btSequentialImpulseConstraintSolverMt::resolveMultipleContactConstraints( const btAlignedObjectArray<int>& consIndices, int batchBegin, int batchEnd )
{
    btScalar leastSquaresResidual = 0.f;
    for ( int iiCons = batchBegin; iiCons < batchEnd; ++iiCons )
    {
        int iCons = consIndices[ iiCons ];
        const btSolverConstraint& solveManifold = m_tmpSolverContactConstraintPool[ iCons ];
        btSolverBody& bodyA = m_tmpSolverBodyPool[ solveManifold.m_solverBodyIdA ];
        btSolverBody& bodyB = m_tmpSolverBodyPool[ solveManifold.m_solverBodyIdB ];
        btScalar residual = resolveSingleConstraintRowLowerLimit( bodyA, bodyB, solveManifold );
        leastSquaresResidual += residual*residual;
    }
    return leastSquaresResidual;
}


btScalar btSequentialImpulseConstraintSolverMt::resolveMultipleContactFrictionConstraints( const btAlignedObjectArray<int>& consIndices, int batchBegin, int batchEnd )
{
    btScalar leastSquaresResidual = 0.f;
    for ( int iiCons = batchBegin; iiCons < batchEnd; ++iiCons )
    {
        int iContact = consIndices[ iiCons ];
        btScalar totalImpulse = m_tmpSolverContactConstraintPool[ iContact ].m_appliedImpulse;

        // apply sliding friction
        if ( totalImpulse > 0.0f )
        {
            int iBegin = iContact * m_numFrictionDirections;
            int iEnd = iBegin + m_numFrictionDirections;
            for ( int iFriction = iBegin; iFriction < iEnd; ++iFriction )
            {
                btSolverConstraint& solveManifold = m_tmpSolverContactFrictionConstraintPool[ iFriction++ ];
                btAssert( solveManifold.m_frictionIndex == iContact );

                solveManifold.m_lowerLimit = -( solveManifold.m_friction*totalImpulse );
                solveManifold.m_upperLimit = solveManifold.m_friction*totalImpulse;

                btSolverBody& bodyA = m_tmpSolverBodyPool[ solveManifold.m_solverBodyIdA ];
                btSolverBody& bodyB = m_tmpSolverBodyPool[ solveManifold.m_solverBodyIdB ];
                btScalar residual = resolveSingleConstraintRowGeneric( bodyA, bodyB, solveManifold );
                leastSquaresResidual += residual*residual;
            }
        }
    }
    return leastSquaresResidual;
}


btScalar btSequentialImpulseConstraintSolverMt::resolveMultipleContactRollingFrictionConstraints( const btAlignedObjectArray<int>& consIndices, int batchBegin, int batchEnd )
{
    btScalar leastSquaresResidual = 0.f;
    for ( int iiCons = batchBegin; iiCons < batchEnd; ++iiCons )
    {
        int iContact = consIndices[ iiCons ];
        int iFirstRollingFriction = m_rollingFrictionIndexTable[ iContact ];
        if ( iFirstRollingFriction >= 0 )
        {
            btScalar totalImpulse = m_tmpSolverContactConstraintPool[ iContact ].m_appliedImpulse;
            // apply rolling friction
            if ( totalImpulse > 0.0f )
            {
                int iBegin = iFirstRollingFriction;
                int iEnd = iBegin + 3;
                for ( int iRollingFric = iBegin; iRollingFric < iEnd; ++iRollingFric )
                {
                    btSolverConstraint& rollingFrictionConstraint = m_tmpSolverContactRollingFrictionConstraintPool[ iRollingFric ];
                    if ( rollingFrictionConstraint.m_frictionIndex != iContact )
                    {
                        break;
                    }
                    btScalar rollingFrictionMagnitude = rollingFrictionConstraint.m_friction*totalImpulse;
                    if ( rollingFrictionMagnitude > rollingFrictionConstraint.m_friction )
                    {
                        rollingFrictionMagnitude = rollingFrictionConstraint.m_friction;
                    }

                    rollingFrictionConstraint.m_lowerLimit = -rollingFrictionMagnitude;
                    rollingFrictionConstraint.m_upperLimit = rollingFrictionMagnitude;

                    btScalar residual = resolveSingleConstraintRowGeneric( m_tmpSolverBodyPool[ rollingFrictionConstraint.m_solverBodyIdA ], m_tmpSolverBodyPool[ rollingFrictionConstraint.m_solverBodyIdB ], rollingFrictionConstraint );
                    leastSquaresResidual += residual*residual;
                }
            }
        }
    }
    return leastSquaresResidual;
}


btScalar btSequentialImpulseConstraintSolverMt::resolveMultipleContactConstraintsInterleaved( const btAlignedObjectArray<int>& contactIndices,
                                                                                          int batchBegin,
                                                                                          int batchEnd
                                                                                          )
{
    btScalar leastSquaresResidual = 0.f;
    int numPoolConstraints = m_tmpSolverContactConstraintPool.size();

    for ( int iiCons = batchBegin; iiCons < batchEnd; iiCons++ )
    {
        btScalar totalImpulse = 0;
        int iContact = contactIndices[ iiCons ];
        // apply penetration constraint
        {
            const btSolverConstraint& solveManifold = m_tmpSolverContactConstraintPool[ iContact ];
            btScalar residual = resolveSingleConstraintRowLowerLimit( m_tmpSolverBodyPool[ solveManifold.m_solverBodyIdA ], m_tmpSolverBodyPool[ solveManifold.m_solverBodyIdB ], solveManifold );
            leastSquaresResidual += residual*residual;
            totalImpulse = solveManifold.m_appliedImpulse;
        }

        // apply sliding friction
        if ( totalImpulse > 0.0f )
        {
            int iBegin = iContact * m_numFrictionDirections;
            int iEnd = iBegin + m_numFrictionDirections;
            for ( int iFriction = iBegin; iFriction < iEnd; ++iFriction )
            {
                btSolverConstraint& solveManifold = m_tmpSolverContactFrictionConstraintPool[ iFriction ];
                btAssert( solveManifold.m_frictionIndex == iContact );

                solveManifold.m_lowerLimit = -( solveManifold.m_friction*totalImpulse );
                solveManifold.m_upperLimit = solveManifold.m_friction*totalImpulse;

                btSolverBody& bodyA = m_tmpSolverBodyPool[ solveManifold.m_solverBodyIdA ];
                btSolverBody& bodyB = m_tmpSolverBodyPool[ solveManifold.m_solverBodyIdB ];
                btScalar residual = resolveSingleConstraintRowGeneric( bodyA, bodyB, solveManifold );
                leastSquaresResidual += residual*residual;
            }
        }

        // apply rolling friction
        int iFirstRollingFriction = m_rollingFrictionIndexTable[ iContact ];
        if ( totalImpulse > 0.0f && iFirstRollingFriction >= 0)
        {
            int iBegin = iFirstRollingFriction;
            int iEnd = iBegin + 3;
            for ( int iRollingFric = iBegin; iRollingFric < iEnd; ++iRollingFric )
            {
                btSolverConstraint& rollingFrictionConstraint = m_tmpSolverContactRollingFrictionConstraintPool[ iRollingFric ];
                if ( rollingFrictionConstraint.m_frictionIndex != iContact )
                {
                    break;
                }
                btScalar rollingFrictionMagnitude = rollingFrictionConstraint.m_friction*totalImpulse;
                if ( rollingFrictionMagnitude > rollingFrictionConstraint.m_friction )
                {
                    rollingFrictionMagnitude = rollingFrictionConstraint.m_friction;
                }

                rollingFrictionConstraint.m_lowerLimit = -rollingFrictionMagnitude;
                rollingFrictionConstraint.m_upperLimit = rollingFrictionMagnitude;

                btScalar residual = resolveSingleConstraintRowGeneric( m_tmpSolverBodyPool[ rollingFrictionConstraint.m_solverBodyIdA ], m_tmpSolverBodyPool[ rollingFrictionConstraint.m_solverBodyIdB ], rollingFrictionConstraint );
                leastSquaresResidual += residual*residual;
            }
        }
    }
    return leastSquaresResidual;
}


void btSequentialImpulseConstraintSolverMt::randomizeBatchedConstraintOrdering( BatchedConstraints* batchedConstraints )
{
    BatchedConstraints& bc = *batchedConstraints;
    // randomize ordering of phases
    for ( int ii = 1; ii < bc.m_phaseOrder.size(); ++ii )
    {
        int iSwap = btRandInt2( ii + 1 );
        bc.m_phaseOrder.swap( ii, iSwap );
    }

    // for each batch,
    for ( int iBatch = 0; iBatch < bc.m_batches.size(); ++iBatch )
    {
        // randomize ordering of constraints within the batch
        const BatchedConstraints::Range& batch = bc.m_batches[ iBatch ];
        for ( int iiCons = batch.begin; iiCons < batch.end; ++iiCons )
        {
            int iSwap = batch.begin + btRandInt2( iiCons - batch.begin + 1 );
            btAssert(iSwap >= batch.begin && iSwap < batch.end);
            bc.m_constraintIndices.swap( iiCons, iSwap );
        }
    }
}


void btSequentialImpulseConstraintSolverMt::randomizeConstraintOrdering(int iteration, int numIterations)
{
    // randomize ordering of joint constraints
    randomizeBatchedConstraintOrdering( &m_batchedJointConstraints );

    //contact/friction constraints are not solved more than numIterations
    if ( iteration < numIterations )
    {
        randomizeBatchedConstraintOrdering( &m_batchedContactConstraints );
    }
}


struct JointSolverLoop : public btIParallelSumBody
{
    btSequentialImpulseConstraintSolverMt* m_solver;
    const BatchedConstraints* m_bc;
    int m_iteration;

    JointSolverLoop( btSequentialImpulseConstraintSolverMt* solver, const BatchedConstraints* bc, int iteration )
    {
        m_solver = solver;
        m_bc = bc;
        m_iteration = iteration;
    }
    btScalar sumLoop( int iBegin, int iEnd ) const BT_OVERRIDE
    {
        BT_PROFILE( "JointSolverLoop" );
        btScalar sum = 0;
        for ( int iBatch = iBegin; iBatch < iEnd; ++iBatch )
        {
            const BatchedConstraints::Range& batch = m_bc->m_batches[ iBatch ];
            sum += m_solver->resolveMultipleJointConstraints( m_bc->m_constraintIndices, batch.begin, batch.end, m_iteration );
        }
        return sum;
    }
};


btScalar btSequentialImpulseConstraintSolverMt::resolveAllJointConstraints(int iteration)
{
    BT_PROFILE( "resolveAllJointConstraints" );
    const BatchedConstraints& batchedCons = m_batchedJointConstraints;
    JointSolverLoop loop( this, &batchedCons, iteration );
    btScalar leastSquaresResidual = 0.f;
    for ( int iiPhase = 0; iiPhase < batchedCons.m_phases.size(); ++iiPhase )
    {
        int iPhase = batchedCons.m_phaseOrder[ iiPhase ];
        const BatchedConstraints::Range& phase = batchedCons.m_phases[ iPhase ];
        int grainSize = 1;
        leastSquaresResidual += btParallelSum( phase.begin, phase.end, grainSize, loop );
    }
    return leastSquaresResidual;
}


struct ContactSolverLoop : public btIParallelSumBody
{
    btSequentialImpulseConstraintSolverMt* m_solver;
    const BatchedConstraints* m_bc;

    ContactSolverLoop( btSequentialImpulseConstraintSolverMt* solver, const BatchedConstraints* bc )
    {
        m_solver = solver;
        m_bc = bc;
    }
    btScalar sumLoop( int iBegin, int iEnd ) const BT_OVERRIDE
    {
        BT_PROFILE( "ContactSolverLoop" );
        btScalar sum = 0;
        for ( int iBatch = iBegin; iBatch < iEnd; ++iBatch )
        {
            const BatchedConstraints::Range& batch = m_bc->m_batches[ iBatch ];
            sum += m_solver->resolveMultipleContactConstraints( m_bc->m_constraintIndices, batch.begin, batch.end );
        }
        return sum;
    }
};


btScalar btSequentialImpulseConstraintSolverMt::resolveAllContactConstraints()
{
    BT_PROFILE( "resolveAllContactConstraints" );
    const BatchedConstraints& batchedCons = m_batchedContactConstraints;
    ContactSolverLoop loop( this, &batchedCons );
    btScalar leastSquaresResidual = 0.f;
    for ( int iiPhase = 0; iiPhase < batchedCons.m_phases.size(); ++iiPhase )
    {
        int iPhase = batchedCons.m_phaseOrder[ iiPhase ];
        const BatchedConstraints::Range& phase = batchedCons.m_phases[ iPhase ];
        int grainSize = 2;
        leastSquaresResidual += btParallelSum( phase.begin, phase.end, grainSize, loop );
    }
    return leastSquaresResidual;
}


struct ContactFrictionSolverLoop : public btIParallelSumBody
{
    btSequentialImpulseConstraintSolverMt* m_solver;
    const BatchedConstraints* m_bc;

    ContactFrictionSolverLoop( btSequentialImpulseConstraintSolverMt* solver, const BatchedConstraints* bc )
    {
        m_solver = solver;
        m_bc = bc;
    }
    btScalar sumLoop( int iBegin, int iEnd ) const BT_OVERRIDE
    {
        BT_PROFILE( "ContactFrictionSolverLoop" );
        btScalar sum = 0;
        for ( int iBatch = iBegin; iBatch < iEnd; ++iBatch )
        {
            const BatchedConstraints::Range& batch = m_bc->m_batches[ iBatch ];
            sum += m_solver->resolveMultipleContactFrictionConstraints( m_bc->m_constraintIndices, batch.begin, batch.end );
        }
        return sum;
    }
};


btScalar btSequentialImpulseConstraintSolverMt::resolveAllContactFrictionConstraints()
{
    BT_PROFILE( "resolveAllContactFrictionConstraints" );
    const BatchedConstraints& batchedCons = m_batchedContactConstraints;
    ContactFrictionSolverLoop loop( this, &batchedCons );
    btScalar leastSquaresResidual = 0.f;
    for ( int iiPhase = 0; iiPhase < batchedCons.m_phases.size(); ++iiPhase )
    {
        int iPhase = batchedCons.m_phaseOrder[ iiPhase ];
        const BatchedConstraints::Range& phase = batchedCons.m_phases[ iPhase ];
        int grainSize = 2;
        leastSquaresResidual += btParallelSum( phase.begin, phase.end, grainSize, loop );
    }
    return leastSquaresResidual;
}


struct InterleavedContactSolverLoop : public btIParallelSumBody
{
    btSequentialImpulseConstraintSolverMt* m_solver;
    const BatchedConstraints* m_bc;

    InterleavedContactSolverLoop( btSequentialImpulseConstraintSolverMt* solver, const BatchedConstraints* bc )
    {
        m_solver = solver;
        m_bc = bc;
    }
    btScalar sumLoop( int iBegin, int iEnd ) const BT_OVERRIDE
    {
        BT_PROFILE( "InterleavedContactSolverLoop" );
        btScalar sum = 0;
        for ( int iBatch = iBegin; iBatch < iEnd; ++iBatch )
        {
            const BatchedConstraints::Range& batch = m_bc->m_batches[ iBatch ];
            sum += m_solver->resolveMultipleContactConstraintsInterleaved( m_bc->m_constraintIndices, batch.begin, batch.end );
        }
        return sum;
    }
};


btScalar btSequentialImpulseConstraintSolverMt::resolveAllContactConstraintsInterleaved()
{
    BT_PROFILE( "resolveAllContactConstraintsInterleaved" );
    const BatchedConstraints& batchedCons = m_batchedContactConstraints;
    InterleavedContactSolverLoop loop( this, &batchedCons );
    btScalar leastSquaresResidual = 0.f;
    for ( int iiPhase = 0; iiPhase < batchedCons.m_phases.size(); ++iiPhase )
    {
        int iPhase = batchedCons.m_phaseOrder[ iiPhase ];
        const BatchedConstraints::Range& phase = batchedCons.m_phases[ iPhase ];
        int grainSize = 1;
        leastSquaresResidual += btParallelSum( phase.begin, phase.end, grainSize, loop );
    }
    return leastSquaresResidual;
}


struct ContactRollingFrictionSolverLoop : public btIParallelSumBody
{
    btSequentialImpulseConstraintSolverMt* m_solver;
    const BatchedConstraints* m_bc;

    ContactRollingFrictionSolverLoop( btSequentialImpulseConstraintSolverMt* solver, const BatchedConstraints* bc )
    {
        m_solver = solver;
        m_bc = bc;
    }
    btScalar sumLoop( int iBegin, int iEnd ) const BT_OVERRIDE
    {
        BT_PROFILE( "ContactFrictionSolverLoop" );
        btScalar sum = 0;
        for ( int iBatch = iBegin; iBatch < iEnd; ++iBatch )
        {
            const BatchedConstraints::Range& batch = m_bc->m_batches[ iBatch ];
            sum += m_solver->resolveMultipleContactRollingFrictionConstraints( m_bc->m_constraintIndices, batch.begin, batch.end );
        }
        return sum;
    }
};


btScalar btSequentialImpulseConstraintSolverMt::resolveAllRollingFrictionConstraints()
{
    BT_PROFILE( "resolveAllRollingFrictionConstraints" );
    btScalar leastSquaresResidual = 0.f;
    //
    // We do not generate batches for rolling friction constraints. We assume that
    // one of two cases is true:
    //
    //  1. either most bodies in the simulation have rolling friction, in which case we can use the
    //     batches for contacts and use a lookup table to translate contact indices to rolling friction
    //     (ignoring any contact indices that don't map to a rolling friction constraint). As long as
    //     most contacts have a corresponding rolling friction constraint, this should parallelize well.
    //
    //  -OR-
    //
    //  2. few bodies in the simulation have rolling friction, so it is not worth trying to use the
    //     batches from contacts as most of the contacts won't have corresponding rolling friction
    //     constraints and most threads would end up doing very little work. Most of the time would
    //     go to threading overhead, so we don't bother with threading.
    //
    int numRollingFrictionPoolConstraints = m_tmpSolverContactRollingFrictionConstraintPool.size();
    if (numRollingFrictionPoolConstraints >= m_tmpSolverContactConstraintPool.size())
    {
        // use batching if there are many rolling friction constraints
        const BatchedConstraints& batchedCons = m_batchedContactConstraints;
        ContactRollingFrictionSolverLoop loop( this, &batchedCons );
        btScalar leastSquaresResidual = 0.f;
        for ( int iiPhase = 0; iiPhase < batchedCons.m_phases.size(); ++iiPhase )
        {
            int iPhase = batchedCons.m_phaseOrder[ iiPhase ];
            const BatchedConstraints::Range& phase = batchedCons.m_phases[ iPhase ];
            int grainSize = 1;
            leastSquaresResidual += btParallelSum( phase.begin, phase.end, grainSize, loop );
        }
    }
    else
    {
        // no batching, also ignores SOLVER_RANDMIZE_ORDER
        for ( int j = 0; j < numRollingFrictionPoolConstraints; j++ )
        {
            btSolverConstraint& rollingFrictionConstraint = m_tmpSolverContactRollingFrictionConstraintPool[ j ];
            if ( rollingFrictionConstraint.m_frictionIndex >= 0 )
            {
                btScalar totalImpulse = m_tmpSolverContactConstraintPool[ rollingFrictionConstraint.m_frictionIndex ].m_appliedImpulse;
                if ( totalImpulse > 0.0f )
                {
                    btScalar rollingFrictionMagnitude = rollingFrictionConstraint.m_friction*totalImpulse;
                    if ( rollingFrictionMagnitude > rollingFrictionConstraint.m_friction )
                        rollingFrictionMagnitude = rollingFrictionConstraint.m_friction;

                    rollingFrictionConstraint.m_lowerLimit = -rollingFrictionMagnitude;
                    rollingFrictionConstraint.m_upperLimit = rollingFrictionMagnitude;

                    btScalar residual = resolveSingleConstraintRowGeneric( m_tmpSolverBodyPool[ rollingFrictionConstraint.m_solverBodyIdA ], m_tmpSolverBodyPool[ rollingFrictionConstraint.m_solverBodyIdB ], rollingFrictionConstraint );
                    leastSquaresResidual += residual*residual;
                }
            }
        }
    }
    return leastSquaresResidual;
}


void btSequentialImpulseConstraintSolverMt::internalWarmstartingWriteContactPoints(int iBegin, int iEnd)
{
    BT_PROFILE("internalWarmstartingWriteContactPoints");
    for ( int iContact = iBegin; iContact < iEnd; ++iContact)
    {
        const btSolverConstraint& contactConstraint = m_tmpSolverContactConstraintPool[ iContact ];
        btManifoldPoint* pt = (btManifoldPoint*) contactConstraint.m_originalContactPoint;
        btAssert( pt );
        pt->m_appliedImpulse = contactConstraint.m_appliedImpulse;
        pt->m_appliedImpulseLateral1 = m_tmpSolverContactFrictionConstraintPool[ contactConstraint.m_frictionIndex ].m_appliedImpulse;
        if ( m_numFrictionDirections == 2 )
        {
            pt->m_appliedImpulseLateral2 = m_tmpSolverContactFrictionConstraintPool[ contactConstraint.m_frictionIndex + 1 ].m_appliedImpulse;
        }
    }
}


struct WarmstartingWriteContactPointsLoop : public btIParallelForBody
{
    btSequentialImpulseConstraintSolverMt* m_solver;

    WarmstartingWriteContactPointsLoop( btSequentialImpulseConstraintSolverMt* solver )
    {
        m_solver = solver;
    }
    void forLoop( int iBegin, int iEnd ) const BT_OVERRIDE
    {
        m_solver->internalWarmstartingWriteContactPoints( iBegin, iEnd );
    }
};


void btSequentialImpulseConstraintSolverMt::warmstartingWriteBackContacts(const btContactSolverInfo& infoGlobal)
{
	BT_PROFILE("warmstartingWriteBackContacts");
	int numPoolConstraints = m_tmpSolverContactConstraintPool.size();
    WarmstartingWriteContactPointsLoop loop( this );
    int grainSize = 500;
    btParallelFor( 0, numPoolConstraints, grainSize, loop );
}


