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


#include "btBatchedConstraints.h"

#include "LinearMath/btIDebugDraw.h"
#include "LinearMath/btMinMax.h"
#include "LinearMath/btStackAlloc.h"
#include "LinearMath/btQuickprof.h"

#include <string.h> //for memset

const int kUnassignedBatch = -1;
const int kNoMerge = -1;

bool btBatchedConstraints::s_debugDrawBatches = false;
btBatchedConstraints::BatchingMethod btBatchedConstraints::s_batchingMethod = btBatchedConstraints::BATCHING_METHOD_DIRECTIONAL;


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


bool btBatchedConstraints::validate(btConstraintArray* constraints, const btAlignedObjectArray<btSolverBody>& bodies) const
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


static void debugDrawSingleBatch( const btBatchedConstraints* bc,
    btConstraintArray* constraints,
    const btAlignedObjectArray<btSolverBody>& bodies,
    int iBatch,
    const btVector3& color,
    const btVector3& offset
    )
{
    if (bc && bc->m_debugDrawer && iBatch < bc->m_batches.size())
    {
        const btBatchedConstraints::Range& b = bc->m_batches[iBatch];
        for (int iiCon = b.begin; iiCon < b.end; ++iiCon)
        {
            int iCon = bc->m_constraintIndices[iiCon];
            const btSolverConstraint& con = constraints->at(iCon);
            int iBody0 = con.m_solverBodyIdA;
            int iBody1 = con.m_solverBodyIdB;
            btVector3 pos0 = bodies[iBody0].getWorldTransform().getOrigin() + offset;
            btVector3 pos1 = bodies[iBody1].getWorldTransform().getOrigin() + offset;
            bc->m_debugDrawer->drawLine(pos0, pos1, color);
        }
    }
}


static void debugDrawPhase( const btBatchedConstraints* bc,
    btConstraintArray* constraints,
    const btAlignedObjectArray<btSolverBody>& bodies,
    int iPhase,
    const btVector3& color0,
    const btVector3& color1,
    const btVector3& offset
    )
{
    BT_PROFILE( "debugDrawPhase" );
    if ( bc && bc->m_debugDrawer && iPhase < bc->m_phases.size() )
    {
        const btBatchedConstraints::Range& phase = bc->m_phases[iPhase];
        for (int iBatch = phase.begin; iBatch < phase.end; ++iBatch)
        {
            float tt = float(iBatch - phase.begin) / float(btMax(1, phase.end - phase.begin - 1));
            btVector3 col = lerp(color0, color1, tt);
            debugDrawSingleBatch(bc, constraints, bodies, iBatch, col, offset);
        }
    }
}


static void debugDrawAllBatches( const btBatchedConstraints* bc,
    btConstraintArray* constraints,
    const btAlignedObjectArray<btSolverBody>& bodies
    )
{
    BT_PROFILE( "debugDrawAllBatches" );
    if ( bc && bc->m_debugDrawer && bc->m_phases.size() > 0 )
    {
        btVector3 bboxMin(BT_LARGE_FLOAT, BT_LARGE_FLOAT, BT_LARGE_FLOAT);
        btVector3 bboxMax = -bboxMin;
        for (int iBody = 0; iBody < bodies.size(); ++iBody)
        {
            const btVector3& pos = bodies[iBody].getWorldTransform().getOrigin();
            bboxMin.setMin(pos);
            bboxMax.setMax(pos);
        }
        btVector3 bboxExtent = bboxMax - bboxMin;
        btVector3 offsetBase = btVector3( 0, bboxExtent.y()*1.1f, 0 );
        btVector3 offsetStep = btVector3( 0, 0, bboxExtent.z()*1.1f );
        btVector3 colorTable[] =
        {
            btVector3(1,0,0),  // R
            btVector3(0,1,0),  // G
            btVector3(0,0,1),  // B
            btVector3(0,1,1),  // C
            btVector3(1,0,1),  // M
            btVector3(1,1,0),  // Y
            btVector3(0,1,0.7),  // bluish green
            btVector3(1,0,0.7),  // pink
            btVector3(1,0.7,0),  // orange
            btVector3(0,0.7,1),  // greenish blue
            btVector3(0.7,0,1),  // purple
            btVector3(0.7,1,0),  // yellowish green
        };
        int numColors = sizeof(colorTable)/sizeof(colorTable[0]);
        int numPhases = bc->m_phases.size();
        for (int iPhase = 0; iPhase < numPhases; ++iPhase)
        {
            btVector3 color0 = colorTable[ iPhase % numColors ];
            btVector3 color1 = color0 * 0.5;
            btVector3 offset = offsetBase + offsetStep*(float(iPhase) - float(numPhases-1)*0.5);
            debugDrawPhase(bc, constraints, bodies, iPhase, color0, color1, offset);
        }
    }
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

    btAlignedObjectArray<bool> bodyConnectivityFlags;
    bodyConnectivityFlags.resize( numBodies, false );

    // for each unassigned constraint,
    for ( int iiCon = curConstraints.size() - 1; iiCon >= 0; --iiCon )
    {
        int iCon = curConstraints[ iiCon ];
        const btBatchedConstraintInfo& conInfo = conInfos[ iCon ];
        btAssert( constraintBatchIds[iCon] == kUnassignedBatch );
        int batchIds[ 2 ] = { kUnassignedBatch, kUnassignedBatch };
        bool isDynamic[ 2 ];
        bool isConnected[ 2 ];
        for ( int iiBody = 0; iiBody < 2; ++iiBody )
        {
            int iBody = conInfo.bodyIds[ iiBody ];
            isDynamic[ iiBody ] = bodyDynamicFlags[iBody];
            isConnected[ iiBody ] = bodyConnectivityFlags[ iBody ];
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
            if (!isConnected[ 1 ])
            {
                assignBatchId = batchIds[ 1 ];
            }
            if ( isDynamic[ 0 ] )
            {
                numBodiesToBeAdded = 1;
            }
        }
        else if ( batchIds[ 1 ] == kUnassignedBatch )
        {
            // one is unassigned, the other not
            if (!isConnected[ 0 ])
            {
                assignBatchId = batchIds[ 0 ];
            }
            if ( isDynamic[ 1 ] )
            {
                numBodiesToBeAdded = 1;
            }
        }
        else if ( batchIds[ 0 ] == batchIds[ 1 ] )
        {
            // both bodies in same batch already
            // if it won't make us exceed max batch size,
            //if ( batches[ batchIds[ 0 ] ].numConstraints < maxBatchSize )
            {
                assignBatchId = batchIds[ 0 ];
            }
        }
        else if ( allowMerging )
        {
            // bodies already assigned to different batches
            btAssert(isDynamic[0] && isDynamic[1]);
            // we could either merge batches, or postpone until next phase
            int batch0Size = batches[ batchIds[ 0 ] ].numConstraints;
            int batch1Size = batches[ batchIds[ 1 ] ].numConstraints;
            bool canMerge = true;
            if (isConnected[0] || isConnected[1])
            {
                canMerge = false;
            }

            btAssert( batches[ batchIds[ 0 ] ].phaseId == batches[ batchIds[ 1 ] ].phaseId );
            if ( canMerge && batch0Size < minBatchSize && batch1Size < minBatchSize && ( batch0Size + batch1Size ) < maxBatchSize )
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
                if ( isDynamic[ iiBody ] )
                {
                    int iBody = conInfo.bodyIds[ iiBody ];
                    if (batchIds[ iiBody ] == assignBatchId)
                    {
                        if (numBodiesToBeAdded > 0)
                        {
                            // prevent this body from causing any more bodies to be added this phase
                            bodyConnectivityFlags[iBody] = true;
                        }
                    }
                    else
                    {
                        bodyBatchIds[ iBody ] = assignBatchId;
                    }
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


inline bool BatchCompare(const btBatchedConstraints::Range& a, const btBatchedConstraints::Range& b)
{
    int lenA = a.end - a.begin;
    int lenB = b.end - b.begin;
    return lenA > lenB;
}


static void writeGrainSizes(btBatchedConstraints* bc)
{
    typedef btBatchedConstraints::Range Range;
    int numPhases = bc->m_phases.size();
    bc->m_phaseGrainSize.resizeNoInitialize(numPhases);
    int numThreads = btGetTaskScheduler()->getNumThreads();
    for (int iPhase = 0; iPhase < numPhases; ++iPhase)
    {
        const Range& phase = bc->m_phases[ iPhase ];
        int numBatches = phase.end - phase.begin;
        float grainSize = floor((0.25f*numBatches / float(numThreads)) + 0.0f);
        bc->m_phaseGrainSize[ iPhase ] = btMax(1, int(grainSize));
    }
}


static void writeOutBatchesSimple(btBatchedConstraints* bc,
    const int* constraintBatchIds,
    int numConstraints,
    const btBatchInfo* batches,
    int* batchWork,
    int numBatches
)
{
    BT_PROFILE("writeOutBatchesSimple");
    typedef btBatchedConstraints::Range Range;
    bc->m_constraintIndices.reserve( numConstraints );
    bc->m_batches.resizeNoInitialize( 0 );
    bc->m_phases.resizeNoInitialize( 0 );

    int maxNumBatches = numBatches;
    {
        int* constraintIdPerBatch = batchWork;  // for each batch, keep an index into the next available slot in the m_constraintIndices array
        int iConstraint = 0;
        int curPhaseBegin = 0;
        int curPhaseId = 0;
        for ( int i = 0; i < numBatches; ++i )
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
    writeGrainSizes(bc);
}


static void writeOutBatches(btBatchedConstraints* bc,
    const int* constraintBatchIds,
    int numConstraints,
    const btBatchInfo* batches,
    int* batchWork,
    int maxNumBatchesPerPhase,
    int numPhases
)
{
    BT_PROFILE("writeOutBatches");
    typedef btBatchedConstraints::Range Range;
    bc->m_constraintIndices.reserve( numConstraints );
    bc->m_batches.resizeNoInitialize( 0 );
    bc->m_phases.resizeNoInitialize( 0 );

    //int maxNumBatches = numPhases * maxNumBatchesPerPhase;
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
    writeGrainSizes(bc);
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
static void setupBatchesGreedyWithMerging(
    btBatchedConstraints* batchedConstraints,
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

    writeOutBatchesSimple(batchedConstraints, &constraintBatchIds[0], numConstraints, &batches[0], &batchWork[0], batches.size());
    btAssert(batchedConstraints->validate(constraints, bodies));
}


static void initTableRecursive(btAlignedObjectArray<int>& table, int numGroups, int x, int y, int step, int phase)
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
// setupBodyLookupMt
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
static void setupBodyLookupMt(
    btBatchedConstraints* batchedConstraints,
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

    writeOutBatchesSimple(batchedConstraints, &constraintBatchIds[0], numConstraints, &batches[0], &batchWork[0], batches.size());
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


static int makePhaseRemappingTable(char* phaseMappingTable, int* numConstraintsPerPhase, int numPhases, int iFlexPhase)
{
    BT_PROFILE("makePhaseRemappingTable");
    // join together small phases to try to create fewer more evenly-sized phases
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
            int iDest = 0;
            if ( iSrc == iFlexPhase )
            {
                int bestScore = maxConstraintsPerPhase;
                for ( int i = 0; i < numPhases; ++i )
                {
                    if (i != iSrc)
                    {
                        int n = numConstraintsPerPhase[ i ];
                        if ( n >= minAllowedConstraintsPerPhase && n < maxConstraintsPerPhase && n < bestScore )
                        {
                            bestScore = n;
                            iDest = i;
                        }
                    }
                }
            }
            else
            {
                float bestScore = BT_LARGE_FLOAT;
                for ( int i = 0; i < numPhases; ++i )
                {
                    int n = numConstraintsPerPhase[ i ];
                    float score = float( n );
                    if ( n >= minAllowedConstraintsPerPhase && n < maxConstraintsPerPhase && score < bestScore )
                    {
                        bestScore = n;
                        iDest = i;
                    }
                }
            }
            {
                // merge phases
                numConstraintsPerPhase[ iDest ] += numConstraintsPerPhase[ iSrc ];
                numConstraintsPerPhase[ iSrc ] = 0;
                phaseMappingTable[ iSrc ] = iDest;
            }
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


static void makePhaseOrderTable(char* phaseOrderTable, const int* numConstraintsPerPhase, int numPhases)
{
    BT_PROFILE("makePhaseOrderTable");
    // create a table that orders phases largest to smallest
    for (int iPhase = 0; iPhase < numPhases; ++iPhase)
    {
        phaseOrderTable[iPhase] = iPhase;
    }
    // insertion sort (numPhases is small so no need for more advanced sorting)
    for (int iPhase = 0; iPhase < numPhases; ++iPhase)
    {
        int maxN = numConstraintsPerPhase[phaseOrderTable[iPhase]];
        int iMaxElem = iPhase;
        for (int iSrc = iPhase + 1; iSrc < numPhases; ++iSrc)
        {
            int n = numConstraintsPerPhase[phaseOrderTable[iSrc]];
            if (n > maxN)
            {
                maxN = n;
                iMaxElem = iSrc;
            }
        }
        if (iMaxElem != iPhase)
        {
            btSwap(phaseOrderTable[iPhase], phaseOrderTable[iMaxElem]);
        }
    }
}


struct CreateBatchesParams
{
    int* constraintBatchIds;
    btBatchInfo* batches;
    const char* constraintPhaseIds;
    const char* phaseMappingTable;
    const btBatchedConstraintInfo* conInfos;
    const bool* bodyDynamicFlags;
    int numBodies;
    int numConstraints;
    int numGroups;
    int minNumBatchesPerPhase;
    int maxNumBatchesPerPhase;
    int minBatchSize;
    int maxBatchSize;
    btBatchedConstraints::CreateBatchesWork* workArray;

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
    btBatchedConstraints::CreateBatchesWork& work = params.workArray[iThread];

    btAlignedObjectArray<int>& curConstraints = work.m_curConstraints;
    curConstraints.reserve(numConstraints);

    btAlignedObjectArray<int>& curFlexConstraints = work.m_curFlexConstraints;  // flex constraints involve one non-dynamic body, so if the dynamic body is unassigned, it can be added to any batch
    curFlexConstraints.reserve(numConstraints);

    btAlignedObjectArray<int>& bodyBatchIds = work.m_bodyBatchIds;

    btUnionFind& unionFind = work.m_unionFind;

    unionFind.reset( numBodies );
    curConstraints.resize( 0 );
    curFlexConstraints.resize( 0 );
    btAssert(params.constraintPhaseIds);
    {
        const char* constraintPhaseIds = params.constraintPhaseIds;
        const char* phaseMappingTable = params.phaseMappingTable;
        for ( int iCon = 0; iCon < numConstraints; ++iCon )
        {
            char unmappedPhase = constraintPhaseIds[iCon];
            if ( unmappedPhase >= 0 )
            {
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
                    if (batches[iNewBatch].numConstraints > minBatchSize)
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
    const char* m_phaseOrderTable;

    CreateBatchesForPhaseLoop( const CreateBatchesParams& params, const char* phaseOrderTable )
    {
        m_params = &params;
        m_phaseOrderTable = phaseOrderTable;
    }
    void forLoop( int iBegin, int iEnd ) const BT_OVERRIDE
    {
        BT_PROFILE( "CreateBatchesForPhaseLoop" );
        for ( int iiPhase = iBegin; iiPhase < iEnd; ++iiPhase )
        {
            int iPhase = m_phaseOrderTable[ iiPhase ];
            createBatchesForPhase(iPhase, *m_params);
        }
    }
};


//
// setupBodyLookupGreedyHybridMt
//
// This method is a hybrid of the setupBodyLookupMt method and the greedyWithMerges method.
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
//   - usually generates 12 evenly sized phases each with a good number of evenly sized batches
//   - each of the ~12 phases is processed in parallel
//
// Cons:
//   - complicated
//   - tuned for a single use case (~1000 bodies in a stack, ~20K contact points)
//
static void setupBodyLookupGreedyHybridMt(
    btBatchedConstraints* batchedConstraints,
    btConstraintArray* constraints,
    const btAlignedObjectArray<btSolverBody>& bodies,
    int minBatchSize,
    int maxBatchSize,
    btBatchedConstraints::CreateBatchesWork* workArray
    )
{
    BT_PROFILE("setupBodyLookupGreedyHybridMt");
    const int numGroups = 16;
    const int groupMask = numGroups - 1;
    const int maxNumPhases = numGroups;
    int numPhases = maxNumPhases;
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

    btAlignedObjectArray<char> constraintPhaseIds;
    constraintPhaseIds.resizeNoInitialize(numConstraints);
    for ( int iCon = 0; iCon < numConstraints; ++iCon )
    {
        const btBatchedConstraintInfo& conInfo = conInfos[ iCon ];
        int body0 = conInfo.bodyIds[ 0 ];
        int body1 = conInfo.bodyIds[ 1 ];
        int groupA = body0 & groupMask;
        int groupB = body1 & groupMask;
        int iPhase = groupPhaseTable[ groupA + groupB*numGroups ];
        constraintPhaseIds[ iCon ] = iPhase;
    }

    char phaseMappingTable[maxNumPhases];
    for (int i = 0; i < maxNumPhases; ++i)
    {
        phaseMappingTable[ i ] = i;
    }

    CreateBatchesParams params;
    params.constraintBatchIds = &constraintBatchIds[0];
    params.batches = &batches[0];
    params.bodyDynamicFlags = &bodyDynamicFlags[0];
    params.conInfos = &conInfos[0];
    params.constraintPhaseIds = &constraintPhaseIds[0];
    params.phaseMappingTable = phaseMappingTable;
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
        char phaseOrderTable[maxNumPhases];
        for (int iPhase = 0; iPhase < numPhases; ++iPhase)
        {
            phaseOrderTable[iPhase] = iPhase;
        }
        // parallel batch creation (deterministic)
        CreateBatchesForPhaseLoop loop(params, phaseOrderTable);
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


static const FaceMap gCubeFaceMaps[ 6 ] =
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
    static const int NUM_AXES = 13;
    char m_data[6*N][N];
    btVector3 m_axes[NUM_AXES];

    CubeMap()
    {
        btVector3 axes[] =
        {
            // cube face axes
            btVector3( 1,0,0 ),
            btVector3( 0,1,0 ),
            btVector3( 0,0,1 ),
            // cube edge axes
            btVector3( 1,0,1 ),
            btVector3( -1,0,1 ),
            btVector3( 0,1,1 ),
            btVector3( 0,-1,1 ),
            btVector3( 1,1,0 ),
            btVector3( -1,1,0 ),
            // cube corner axes
            btVector3( 1,1,1 ),
            btVector3( -1,1,1 ),
            btVector3( 1,-1,1 ),
            btVector3( -1,-1,1 ),
        };
        for (int i = 0; i < NUM_AXES; ++i)
        {
            m_axes[i] = axes[i].normalized();
        }
        buildCubeMapLookupFromAxes(m_axes, NUM_AXES);
    }
    char lookupFromDirection(const btVector3& dir) const
    {
        btVector3 uv = normalizedDirectionToCubeUv( dir );
        int ix = int( uv.x() * N );
        int iy = int( uv.y() * N );
        btClamp(ix, 0, N-1);
        btClamp(iy, 0, 6*N-1);
        btAssert( ix >= 0 && ix < N );
        btAssert( iy >= 0 && iy < 6*N );
        char ret = m_data[ iy ][ ix ];
        return ret;
    }
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
                m_data[ iy ][ ix ] = axis;
            }
        }
    }
};


const int MAX_NUM_PHASES = CubeMap::NUM_AXES*2 + 1;

static CubeMap g_cubeMap;


struct AssignConstraintsToPhasesParams
{
    const CubeMap* cubeMap;
    btBatchedConstraintInfo* conInfos;
    const btSolverConstraint* constraints;
    const btVector3* bodyPositions;
    const bool* bodyDynamicFlags;
    int numAxes;
    char* constraintPhaseIds;
    int* constraintBatchIds;

    AssignConstraintsToPhasesParams(char* _constraintPhaseIds,
        int* _constraintBatchIds,
        btBatchedConstraintInfo* _conInfos,
        const btSolverConstraint* _constraints,
        const CubeMap* _cubeMap,
        const btVector3* _bodyPositions,
        const bool* _bodyDynamicFlags,
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
        bodyDynamicFlags = _bodyDynamicFlags;
    }
};


static void assignConstraintsToPhases(const AssignConstraintsToPhasesParams& params, int begin, int end)
{
    BT_PROFILE("assignConstraintsToPhases");
    for (int iCon = begin; iCon < end; ++iCon)
    {
        btBatchedConstraintInfo& conInfo = params.conInfos[iCon];
        const btSolverConstraint& con = params.constraints[iCon];
        int body0 = con.m_solverBodyIdA;
        int body1 = con.m_solverBodyIdB;
        conInfo.bodyIds[0] = body0;
        conInfo.bodyIds[1] = body1;
        int iPhase = MAX_NUM_PHASES - 1;
        if ( params.bodyDynamicFlags[ body0 ] && params.bodyDynamicFlags[ body1 ] )
        {
            btVector3 v = params.bodyPositions[ body1 ] - params.bodyPositions[ body0 ];
            btVector3 vDir = v.normalized();

            int axis = params.cubeMap->lookupFromDirection( vDir );
            btAssert( axis >= 0 && axis < params.numAxes );

            iPhase = axis * 2; // + ((body0)&1);
            btAssert( iPhase >= 0 && iPhase < MAX_NUM_PHASES );
        }
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


static int makePhaseRemappingTableDirectional(char* phaseMappingTable, int* numConstraintsPerPhase, int numPhases, int iFlexPhase)
{
    BT_PROFILE("makePhaseRemappingTableDirectional");
    // join together small phases to try to create fewer more evenly-sized phases
    btVector3 phaseAxis[ MAX_NUM_PHASES ];
    for (int i = 0; i < numPhases; ++i)
    {
        phaseMappingTable[i] = i;
        if (i == iFlexPhase)
        {
            // flex phase is for constraints with 1 dynamic body -- no axis
            phaseAxis[ i ] = btVector3(0,0,0);
        }
        else
        {
            int iAxis = i >> 1;
            btAssert( iAxis >= 0 && iAxis < CubeMap::NUM_AXES );
            phaseAxis[ i ] = g_cubeMap.m_axes[ iAxis ];
        }
    }
    int maxConstraintsPerPhase = 0;
    for (int i = 0; i < numPhases; ++i)
    {
        maxConstraintsPerPhase = btMax(maxConstraintsPerPhase, numConstraintsPerPhase[i]);
    }
    int minAllowedConstraintsPerPhase = maxConstraintsPerPhase/4;
    for (int iSrc = 0; iSrc < numPhases; ++iSrc)
    {
        if (numConstraintsPerPhase[iSrc] > 0 && numConstraintsPerPhase[iSrc] < minAllowedConstraintsPerPhase && iSrc != iFlexPhase)
        {
            btVector3 srcDir = phaseAxis[ iSrc ];
            // find dest phase to merge into
            // find smallest valid phase
            float bestDot = 0.0f;
            int iDest = 0;
            {
                float bestScore = BT_LARGE_FLOAT;
                for ( int i = 0; i < numPhases; ++i )
                {
                    int n = numConstraintsPerPhase[ i ];
                    btVector3 destDir = phaseAxis[ i ];
                    float dot = fabs( srcDir.dot( destDir ) );
                    float score = ( 1.0f - dot )*float( n );
                    if ( n >= minAllowedConstraintsPerPhase && n < maxConstraintsPerPhase && score < bestScore && i != iFlexPhase )
                    {
                        bestScore = n;
                        bestDot = dot;
                        iDest = i;
                    }
                }
            }
            const float dotThresh =
                0.8164f;
                //0.707f;
                //0.5773f;
            if (bestDot >= dotThresh)
            {
                // merge phases
                btVector3 srcDir = phaseAxis[ iSrc ];
                btVector3 destDir = phaseAxis[ iDest ];
                if (srcDir.dot(destDir) < 0.0f)
                {
                    srcDir = -srcDir;
                }
                btVector3 dir = srcDir * float( numConstraintsPerPhase[ iSrc ] ) + destDir *float( numConstraintsPerPhase[ iDest ] );
                dir.normalize();
                phaseAxis[ iDest ] = dir;
                numConstraintsPerPhase[ iDest ] += numConstraintsPerPhase[ iSrc ];
                numConstraintsPerPhase[ iSrc ] = 0;
                phaseMappingTable[ iSrc ] = iDest;
            }
        }
    }
    // now merge the flex phase into the smallest (non-zero) other phase
    if ( numConstraintsPerPhase[ iFlexPhase ] > 0 )
    {
        int bestScore = maxConstraintsPerPhase;
        int iDest = -1;
        bool foundDestPhase = false;
        for ( int i = 0; i < numPhases; ++i )
        {
            if ( i != iFlexPhase )
            {
                int n = numConstraintsPerPhase[ i ];
                if ( n > 0 && n < maxConstraintsPerPhase && n < bestScore )
                {
                    bestScore = n;
                    iDest = i;
                    foundDestPhase = true;
                }
            }
        }
        // if suitable dest is found,
        if (foundDestPhase )
        {
            int iSrc = iFlexPhase;
            numConstraintsPerPhase[ iDest ] += numConstraintsPerPhase[ iSrc ];
            numConstraintsPerPhase[ iSrc ] = 0;
            phaseMappingTable[ iSrc ] = iDest;
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
    btAssert(numActualPhases > 0);
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
// setupDirectionalBatchesMt -- we want a way to ensure that a single batch does not get too large
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
//
//    We pre-generate a small cubemap, which we can use to map from the direction
//    vector between bodies to initial phase index.
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
//    single batch. This issue has come up in testing.
//

/*

 New idea:

 So the problem is that by greedily combining constraints into batches, we are
 creating batches that are too large because they have too much connectivity
 (too much like a web) in places and end up touching most of the available bodies.
 Then there aren't enough constraints left over to form a useful number of evenly sized batches.

 When building batches, only allow each body to be connected to 2 other bodies.
 This is because we are looking to form batches that comprise a strand or filament-like
 structure, and in a strand, each body is only connected to, at most 2 others (with the
 bodies on the ends only connected to one other).
 The hope is that by limiting the body connectivity, the batches will naturally
 form into stands.
 Recall that with contact constraints it is not unusual to have 3 or 4 of them connecting
 the same pair of bodies, and we want to group those all in the same batch to avoid needing
 more phases than necessary.
 So the rule is, if a constraint comes along that would add too much connectivity, then we
 transfer that constraint to another phase.

*/
static void setupDirectionalBatchesMt(
    btBatchedConstraints* batchedConstraints,
    btAlignedObjectArray<char>* scratchMemory,
    btConstraintArray* constraints,
    const btAlignedObjectArray<btSolverBody>& bodies,
    int minBatchSize,
    int maxBatchSize,
    btBatchedConstraints::CreateBatchesWork* workArray
    )
{
    BT_PROFILE("setupDirectionalBatchesMt");
    const int numAxes = CubeMap::NUM_AXES;
    const CubeMap& cubeMap = g_cubeMap;

    const int numPhases = MAX_NUM_PHASES;
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

    bool assignInParallel = true;
    {
        AssignConstraintsToPhasesParams params( &constraintPhaseIds[ 0 ],
            constraintBatchIds,
            &conInfos[ 0 ],
            &constraints->at(0),
            &cubeMap,
            &bodyPositions[ 0 ],
            &bodyDynamicFlags[ 0 ],
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

    int numConstraintsPerPhase[numPhases];
    for (int i = 0; i < numPhases; ++i)
    {
        numConstraintsPerPhase[i] = 0;
    }

    for (int iCon = 0; iCon < numConstraints; ++iCon)
    {
        int iPhase = constraintPhaseIds[iCon];
        numConstraintsPerPhase[iPhase]++;
    }
    char phaseMappingTable[numPhases];
    int numActualPhases = makePhaseRemappingTableDirectional(phaseMappingTable, numConstraintsPerPhase, numPhases, MAX_NUM_PHASES-1);

    CreateBatchesParams params;
    params.constraintBatchIds = &constraintBatchIds[0];
    params.batches = &batches[0];
    params.bodyDynamicFlags = &bodyDynamicFlags[0];
    params.conInfos = &conInfos[0];
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
        // dispatch phases in sorted order (largest to smallest) because it improves thread utilization for some task schedulers
        char phaseOrderTable[numPhases];
        makePhaseOrderTable(phaseOrderTable, numConstraintsPerPhase, numPhases);
        // parallel batch creation (deterministic)
        CreateBatchesForPhaseLoop loop(params, phaseOrderTable);
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
    btBatchedConstraints* bc,
    int numConstraints
)
{
    BT_PROFILE("setupSingleBatch");
    typedef btBatchedConstraints::Range Range;

    bc->m_constraintIndices.resize( numConstraints );
    for ( int i = 0; i < numConstraints; ++i )
    {
        bc->m_constraintIndices[ i ] = i;
    }

    bc->m_batches.resizeNoInitialize( 0 );
    bc->m_phases.resizeNoInitialize( 0 );
    bc->m_phaseOrder.resizeNoInitialize( 0 );
    bc->m_phaseGrainSize.resizeNoInitialize( 0 );

    if (numConstraints > 0)
    {
        bc->m_batches.push_back( Range( 0, numConstraints ) );
        bc->m_phases.push_back( Range( 0, 1 ) );
        bc->m_phaseOrder.push_back(0);
        bc->m_phaseGrainSize.push_back(1);
    }
}


void btBatchedConstraints::setup(
    btConstraintArray* constraints,
    const btAlignedObjectArray<btSolverBody>& bodies,
    int minBatchSize,
    int maxBatchSize,
    btAlignedObjectArray<char>* scratchMemory,
    CreateBatchesWork* workArray
    )
{
    if (constraints->size() >= minBatchSize*4)
    {
        switch (s_batchingMethod)
        {
        case BATCHING_METHOD_GREEDY:
            setupBatchesGreedyWithMerging( this, constraints, bodies, minBatchSize, maxBatchSize );
            break;

        case BATCHING_METHOD_BODY_LOOKUP:
            setupBodyLookupMt( this, constraints, bodies, minBatchSize, maxBatchSize );
            break;

        case BATCHING_METHOD_BODY_LOOKUP_GREEDY_HYBRID:
            setupBodyLookupGreedyHybridMt( this, constraints, bodies, minBatchSize, maxBatchSize, workArray );
            break;

        case BATCHING_METHOD_DIRECTIONAL:
            setupDirectionalBatchesMt( this, scratchMemory, constraints, bodies, minBatchSize, maxBatchSize, workArray );
            break;
        }
        if (s_debugDrawBatches)
        {
            debugDrawAllBatches( this, constraints, bodies );
        }
    }
    else
    {
        setupSingleBatch( this, constraints->size() );
    }
}


