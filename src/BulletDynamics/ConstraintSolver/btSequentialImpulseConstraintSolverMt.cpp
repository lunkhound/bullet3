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
#include "BulletCollision/NarrowPhaseCollision/btPersistentManifold.h"

#include "LinearMath/btIDebugDraw.h"
//#include "LinearMath/btCpuFeatureUtility.h"

//#include "btJacobianEntry.h"
#include "LinearMath/btMinMax.h"
#include "BulletDynamics/ConstraintSolver/btTypedConstraint.h"
#include <new>
#include "LinearMath/btStackAlloc.h"
#include "LinearMath/btQuickprof.h"
//#include "btSolverBody.h"
//#include "btSolverConstraint.h"
#include "LinearMath/btAlignedObjectArray.h"
#include <string.h> //for memset
#include "BulletDynamics/Dynamics/btRigidBody.h"


bool btSequentialImpulseConstraintSolverMt::sAllowNestedParallelForLoops = false;  // some task schedulers don't like nested loops
int btSequentialImpulseConstraintSolverMt::sMinimumContactManifoldsForBatching = 1500;
int btSequentialImpulseConstraintSolverMt::sMinBatchSize = 20;
int btSequentialImpulseConstraintSolverMt::sMaxBatchSize = 60;


struct btBatchedConstraintInfo
{
    int constraintId;
    int batchId;
    int phaseId;
    int bodyIds[2];

    class SortPredicate
    {
    public:
        bool operator() ( const btBatchedConstraintInfo& lhs, const btBatchedConstraintInfo& rhs ) const
        {
            return lhs.batchId < rhs.batchId;
        }
    };
};


struct btBatchInfo
{
    int numConstraints;
    int mergeIndex;
    int phaseId;

    btBatchInfo() : numConstraints(0), mergeIndex(-1), phaseId(-1) {}
};

void btSequentialImpulseConstraintSolverMt::BatchedConstraints::setup(
    btConstraintArray* constraints,
    const btAlignedObjectArray<btSolverBody>& bodies,
    int minBatchSize,
    int maxBatchSize
    )
{
    int numConstraints = constraints->size();

    const int kNonDynamicBodyId = -1;
    const int kUnassignedBatch = -1;
    const int kUnassignedPhase = -1;
    const int kNoMerge = -1;

    btAlignedObjectArray<btBatchedConstraintInfo> conInfos;
    conInfos.resize(numConstraints);
    for (int i = 0; i < numConstraints; ++i)
    {
        btBatchedConstraintInfo& conInfo = conInfos[i];
        conInfo.constraintId = i;
        conInfo.batchId = kUnassignedBatch;
        conInfo.phaseId = kUnassignedPhase;
        const btSolverConstraint& con = constraints->at( i );
        conInfo.bodyIds[0] = con.m_solverBodyIdA;
        conInfo.bodyIds[1] = con.m_solverBodyIdB;
        for (int iiBody = 0; iiBody < 2; iiBody++)
        {
            int iBody = conInfo.bodyIds[iiBody];
            const btSolverBody& body = bodies[ iBody ];
            // if not dynamic body
            if ( body.internalGetInvMass().x() == btScalar( 0 ) )
            {
                // static/kinematic bodies should not be batched -- because they are not mutated, they
                // can appear in more than one batch
                conInfo.bodyIds[iiBody] = kNonDynamicBodyId;
            }
        }
    }

    btAlignedObjectArray<int> bodyBatchId;
    btAlignedObjectArray<btBatchInfo> curBatches;
    curBatches.reserve( numConstraints/4 ); // can never have more batches in a phase than num bodies
    btAlignedObjectArray<int> curConstraints;
    curConstraints.resize(numConstraints);
    for (int i = 0; i < numConstraints; ++i)
    {
        curConstraints[i] = i;
    }

    int curBatchId = 0;
    int curPhaseId = 0;
    while (curConstraints.size() > 0)
    {
        // begin new phase
        // mark all bodies unassigned
        bodyBatchId.resize(0);
        bodyBatchId.resize( bodies.size(), kUnassignedBatch );
        int curPhaseBegin = curBatches.size();
        // for each unassigned constraint,
        for (int iiCon = curConstraints.size() - 1; iiCon >= 0; --iiCon)
        {
            int iCon = curConstraints[iiCon];
            btBatchedConstraintInfo& conInfo = conInfos[ iCon ];
            btAssert( conInfo.batchId == kUnassignedBatch );
            btAssert( conInfo.phaseId == kUnassignedPhase );
            int batchIds[2] = {kUnassignedBatch, kUnassignedBatch};
            for ( int iiBody = 0; iiBody < 2; ++iiBody )
            {
                int iBody = conInfo.bodyIds[iiBody];
                if (iBody != kNonDynamicBodyId)
                {
                    int iBatch = bodyBatchId[ iBody ];
                    // check for merge remapping
                    while (iBatch != kUnassignedBatch && curBatches[iBatch].mergeIndex != kNoMerge)
                    {
                        iBatch = curBatches[iBatch].mergeIndex;
                    }
                    batchIds[ iiBody ] = iBatch;
                }
            }
            int assignBatchId = kUnassignedBatch;
            int numBodiesToBeAdded = 0;
            // if both unassigned,
            if ( batchIds[0] == kUnassignedBatch && batchIds[1] == kUnassignedBatch )
            {
                // create a new batch
                assignBatchId = curBatchId++;
                btBatchInfo batch;
                batch.phaseId = curPhaseId;
                curBatches.push_back(batch);
                btAssert( assignBatchId < curBatches.size() );
            }
            else if ( batchIds[ 0 ] == kUnassignedBatch )
            {
                // one is unassigned, the other not
                assignBatchId = batchIds[ 1 ];
                numBodiesToBeAdded = 1;
            }
            else if ( batchIds[ 1 ] == kUnassignedBatch )
            {
                // one is unassigned, the other not
                assignBatchId = batchIds[ 0 ];
                numBodiesToBeAdded = 1;
            }
            else if (batchIds[ 0 ] == batchIds[ 1 ])
            {
                // both bodies in same batch already
                // if it won't make us exceed max batch size,
                if (curBatches[batchIds[ 0 ]].numConstraints < maxBatchSize)
                {
                    assignBatchId = batchIds[ 0 ];
                }
            }
            else
            {
                // bodies already assigned to different batches
                // we could either merge batches, or postpone until next phase
                int batch0Size = curBatches[ batchIds[ 0 ] ].numConstraints;
                int batch1Size = curBatches[ batchIds[ 1 ] ].numConstraints;
                btAssert( curBatches[ batchIds[ 0 ] ].phaseId == curBatches[ batchIds[ 1 ] ].phaseId );
                if (batch0Size < minBatchSize && batch1Size < minBatchSize && (batch0Size+batch1Size) < maxBatchSize)
                {
                    // merge higher index batch into lower index batch
                    assignBatchId = btMin( batchIds[0], batchIds[1] );
                    int otherBatch = btMax( batchIds[0], batchIds[1] );
                    curBatches[assignBatchId].numConstraints += curBatches[otherBatch].numConstraints;
                    btAssert(curBatches[assignBatchId].mergeIndex == kNoMerge);
                    curBatches[otherBatch].numConstraints = 0;
                    curBatches[otherBatch].mergeIndex = assignBatchId;
                }
            }
            // if batch has reached minBatchSize and this constraint would add a body to the batch,
            if (assignBatchId != kUnassignedBatch && numBodiesToBeAdded > 0 && curBatches[assignBatchId].numConstraints >= minBatchSize)
            {
                // don't assign
                assignBatchId = kUnassignedBatch;
            }
            if ( assignBatchId != kUnassignedBatch )
            {
                conInfo.batchId = assignBatchId;
                conInfo.phaseId = curPhaseId;
                curBatches[ assignBatchId ].numConstraints++;
                // swap and pop
                curConstraints.swap(iiCon, curConstraints.size()-1);
                curConstraints.pop_back();
                // make sure bodies are assigned
                for ( int iiBody = 0; iiBody < 2; ++iiBody )
                {
                    if ( batchIds[ iiBody ] != assignBatchId )
                    {
                        int iBody = conInfo.bodyIds[ iiBody ];
                        if ( iBody != kNonDynamicBodyId )
                        {
                            bodyBatchId[ iBody ] = assignBatchId;
                        }
                    }
                }
            }
        }
        // merge small batches into others
        for (int iBatch = curBatches.size()-1; iBatch >= curPhaseBegin; --iBatch)
        {
            btBatchInfo& batch = curBatches[ iBatch ];
            if (batch.mergeIndex == kNoMerge && batch.numConstraints < minBatchSize)
            {
                for ( int iDestBatch = iBatch - 1; iDestBatch >= curPhaseBegin; --iDestBatch )
                {
                    btBatchInfo& destBatch = curBatches[ iBatch ];
                    if (destBatch.mergeIndex == kNoMerge && (destBatch.numConstraints + batch.numConstraints) < maxBatchSize)
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
        for (int iBatch = curPhaseBegin; iBatch < curBatches.size(); ++iBatch)
        {
            btBatchInfo& batch = curBatches[ iBatch ];
            if (batch.mergeIndex != kNoMerge)
            {
                int iMergeDest = curBatches[ batch.mergeIndex ].mergeIndex;
                // follow trail of merges to the end
                while ( iMergeDest != kNoMerge )
                {
                    int iNext = curBatches[ iMergeDest ].mergeIndex;
                    if (iNext == kNoMerge)
                    {
                        batch.mergeIndex = iMergeDest;
                        break;
                    }
                    iMergeDest = iNext;
                }
            }
        }
        curPhaseId++;
    }
    // all constraints have been assigned a batchId and phaseId
    // update batchIds to account for merges
    for (int i = 0; i < numConstraints; ++i)
    {
        btBatchedConstraintInfo& conInfo = conInfos[i];
        int iBatch = conInfo.batchId;
        // if this constraint references a batch that was merged into another batch
        if (curBatches[iBatch].mergeIndex != kNoMerge)
        {
            // update batchId
            conInfo.batchId = curBatches[iBatch].mergeIndex;
        }
    }
    // sort them by phase and batch
    conInfos.quickSort(btBatchedConstraintInfo::SortPredicate());
    mConstraintIndices.reserve(numConstraints);
    mConstraintIndices.resize(0);
    mBatches.resize(0);
    mPhases.resize(0);

    int curBatchBegin = 0;
    int curBatchEnd = 0;
    int curPhaseBegin = 0;
    curBatchId = 0;
    curPhaseId = 0;
    for (int iiCon = 0; iiCon < conInfos.size(); ++iiCon)
    {
        const btBatchedConstraintInfo& conInfo = conInfos[ iiCon ];
        if (conInfo.batchId != curBatchId)
        {
            // output batch
            mBatches.push_back( Range( curBatchBegin, curBatchEnd ) );
            curBatchId = conInfo.batchId;
            curBatchBegin = curBatchEnd = mConstraintIndices.size();
            if (conInfo.phaseId != curPhaseId)
            {
                mPhases.push_back( Range( curPhaseBegin, mBatches.size() ) );
                curPhaseBegin = mBatches.size();
                curPhaseId++;
            }
            btAssert( conInfo.batchId == curBatchId );
            btAssert( conInfo.phaseId == curPhaseId );
        }
        curBatchEnd++;
        mConstraintIndices.push_back(conInfo.constraintId);
    }
    {
        // output last batch
        mBatches.push_back( Range( curBatchBegin, curBatchEnd ) );
        mPhases.push_back( Range( curPhaseBegin, mBatches.size() ) );
    }
    mPhaseOrder.resize(mPhases.size());
    for (int i = 0; i < mPhases.size(); ++i)
    {
        mPhaseOrder[i] = i;
    }
#ifdef _DEBUG
    // verify coloring of bodies, that no body is touched by more than one batch in any given phase
    for (int iPhase = 0; iPhase < mPhases.size(); ++iPhase)
    {
        bodyBatchId.resize(0);
        bodyBatchId.resize( bodies.size(), kUnassignedBatch );
        const Range& phase = mPhases[iPhase];
        for (int iBatch = phase.begin; iBatch < phase.end; ++iBatch)
        {
            const Range& batch = mBatches[iBatch];
            for (int iiCons = batch.begin; iiCons < batch.end; ++iiCons)
            {
                int iCons = mConstraintIndices[iiCons];
                const btSolverConstraint& cons = constraints->at(iCons);
                const btSolverBody& bodyA = bodies[cons.m_solverBodyIdA];
                const btSolverBody& bodyB = bodies[cons.m_solverBodyIdB];
                if (! bodyA.internalGetInvMass().isZero())
                {
                    if (bodyBatchId[cons.m_solverBodyIdA] == kUnassignedBatch)
                    {
                        bodyBatchId[cons.m_solverBodyIdA] = iBatch;
                    }
                    btAssert(bodyBatchId[cons.m_solverBodyIdA] == iBatch);
                }
                if (! bodyB.internalGetInvMass().isZero())
                {
                    if (bodyBatchId[cons.m_solverBodyIdB] == kUnassignedBatch)
                    {
                        bodyBatchId[cons.m_solverBodyIdB] = iBatch;
                    }
                    btAssert(bodyBatchId[cons.m_solverBodyIdB] == iBatch);
                }
            }
        }
    }
#endif
}


void btSequentialImpulseConstraintSolverMt::BatchedConstraints::setup2(
    btConstraintArray* constraints,
    const btAlignedObjectArray<btSolverBody>& bodies,
    int minBatchSize,
    int maxBatchSize
    )
{
    int numConstraints = constraints->size();

    btAlignedObjectArray<int> bodyBatchId;
    mConstraintIndices.reserve(numConstraints);
    mConstraintIndices.resize(0);
    mBatches.resize(0);
    mPhases.resize(0);
    btAlignedObjectArray<int> constraintIndices1;
    btAlignedObjectArray<int> constraintIndices2;
    btAlignedObjectArray<int>* curConstraints = &constraintIndices1;
    btAlignedObjectArray<int>* pendingConstraints = &constraintIndices2;
    curConstraints->resize(numConstraints);
    pendingConstraints->reserve(numConstraints);
    for (int i = 0; i < numConstraints; ++i)
    {
        curConstraints->at(i) = i;
    }

    int curBatchId = 0;
    int curBatchBegin = 0;
    int curBatchEnd = 0;
    int phaseBegin = 0;
    const int kUnassignedBatch = -1;

    while (curConstraints->size() > 0)
    {
        // mark all bodies unassigned
        bodyBatchId.resize(0);
        bodyBatchId.resize( bodies.size(), kUnassignedBatch );
        // for each constraint in curConstraints:
        for ( int iiCon = 0; iiCon < curConstraints->size(); ++iiCon )
        {
            bool constraintReady = true;
            int iCon = curConstraints->at(iiCon);
            const btSolverConstraint& con = constraints->at(iCon);
            int conBodyIndices[2] = {con.m_solverBodyIdA, con.m_solverBodyIdB};
            int bodiesNotInCurrentBatch = 0;
            for (int iiBody = 0; iiBody < 2; ++iiBody)
            {
                int iBody = conBodyIndices[iiBody];
                const btSolverBody& body = bodies[iBody];
                // if dynamic body
                if (body.internalGetInvMass().x() != btScalar(0))
                {
                    if (bodyBatchId[iBody] == kUnassignedBatch)
                    {
                        bodiesNotInCurrentBatch++;
                    }
                    else if (bodyBatchId[iBody] != curBatchId)
                    {
                        // constraint is not ready because this body has been assigned to a previous batch
                        constraintReady = false;
                        break;
                    }
                }
            }
            if ( constraintReady )
            {
                for ( int iiBody = 0; iiBody < 2; ++iiBody )
                {
                    int iBody = conBodyIndices[ iiBody ];
                    const btSolverBody& body = bodies[ iBody ];
                    // if dynamic body
                    if ( bodyBatchId[ iBody ] == kUnassignedBatch && body.internalGetInvMass().x() != btScalar( 0 ) )
                    {
                        bodyBatchId[ iBody ] = curBatchId;
                    }
                }
                mConstraintIndices.push_back( iCon );
                curBatchEnd++;
                int curBatchSize = curBatchEnd - curBatchBegin;
                if (curBatchSize >= maxBatchSize)
                {
                    // output batch
                    mBatches.push_back(Range(curBatchBegin, curBatchEnd));
                    curBatchId++;
                    curBatchBegin = curBatchEnd = mConstraintIndices.size();
                }
            }
            else
            {
                pendingConstraints->push_back( iCon );
            }
        }
        if (curBatchEnd != curBatchBegin)
        {
            // output batch
            mBatches.push_back( Range( curBatchBegin, curBatchEnd ) );
            curBatchId++;
            curBatchBegin = curBatchEnd = mConstraintIndices.size();
        }
        // output phase
        mPhases.push_back( Range( phaseBegin, mBatches.size() ) );
        phaseBegin = mBatches.size();
        btSwap(pendingConstraints, curConstraints);
        pendingConstraints->resize(0);
    }
    mPhaseOrder.resize(mPhases.size());
    for (int i = 0; i < mPhases.size(); ++i)
    {
        mPhaseOrder[i] = i;
    }
}


btSequentialImpulseConstraintSolverMt::btSequentialImpulseConstraintSolverMt()
{
    m_numFrictionDirections = 1;
    m_useBatching = false;
}

btSequentialImpulseConstraintSolverMt::~btSequentialImpulseConstraintSolverMt()
{
}

void btSequentialImpulseConstraintSolverMt::setupBatchedConstraints()
{
    BT_PROFILE("setupBatchedConstraints");
    int minBatchSize = sMinBatchSize;
    int maxBatchSize = sMaxBatchSize;
    m_batchedNonContactConstraints.setup( &m_tmpSolverNonContactConstraintPool, m_tmpSolverBodyPool, minBatchSize, maxBatchSize );
    m_batchedContactConstraints.setup( &m_tmpSolverContactConstraintPool, m_tmpSolverBodyPool, minBatchSize, maxBatchSize );
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
    btSequentialImpulseConstraintSolver::solveGroupCacheFriendlySetup( bodies,
                                                                       numBodies,
                                                                       manifoldPtr,
                                                                       numManifolds,
                                                                       constraints,
                                                                       numConstraints,
                                                                       infoGlobal,
                                                                       debugDrawer
                                                                       );
    m_useBatching = false;
    if ( numManifolds >= sMinimumContactManifoldsForBatching &&
        (sAllowNestedParallelForLoops || !btThreadsAreRunning())
        )
    {
        m_useBatching = true;
    }
    if ( m_useBatching )
    {
        setupBatchedConstraints();
    }
    return 0.0f;
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
        leastSquaresResidual += resolveAllNonContactConstraints(iteration);

		if (iteration< infoGlobal.m_numIterations)
		{
            // this loop is only used for cone-twist constraints,
            // it would be nice to skip this loop if none of the constraints need it
            if ( true )
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
                // solve all contact and contact-friction constraints interleaved
                leastSquaresResidual += resolveAllContactConstraintsInterleaved();
			}
			else//SOLVER_INTERLEAVE_CONTACT_AND_FRICTION_CONSTRAINTS
			{
                // don't interleave them
				// solve all contact constraints
                leastSquaresResidual += resolveAllContactConstraints();

				// solve all contact friction constraints
                leastSquaresResidual += resolveAllContactFrictionConstraints();
			}

            // rolling friction
            leastSquaresResidual += resolveAllRollingFrictionConstraints();
		}
	}
    return leastSquaresResidual;
}

btScalar btSequentialImpulseConstraintSolverMt::resolveMultipleNonContactConstraints( const btAlignedObjectArray<int>& consIndices, int batchBegin, int batchEnd, int iteration )
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
            btScalar residual = resolveSingleConstraintRowGenericSIMD( bodyA, bodyB, constraint );
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
        btScalar residual = resolveSingleConstraintRowLowerLimitSIMD( bodyA, bodyB, solveManifold );
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

        if ( totalImpulse > btScalar( 0 ) )
        {
            int iFriction = iContact * m_numFrictionDirections;
            for ( int ii = m_numFrictionDirections; ii > 0; --ii )
            {
                btSolverConstraint& solveManifold = m_tmpSolverContactFrictionConstraintPool[ iFriction++ ];
                btAssert( solveManifold.m_frictionIndex == iContact );

                solveManifold.m_lowerLimit = -( solveManifold.m_friction*totalImpulse );
                solveManifold.m_upperLimit = solveManifold.m_friction*totalImpulse;

                btSolverBody& bodyA = m_tmpSolverBodyPool[ solveManifold.m_solverBodyIdA ];
                btSolverBody& bodyB = m_tmpSolverBodyPool[ solveManifold.m_solverBodyIdB ];
                btScalar residual = resolveSingleConstraintRowGenericSIMD( bodyA, bodyB, solveManifold );
                leastSquaresResidual += residual*residual;
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
        {
            const btSolverConstraint& solveManifold = m_tmpSolverContactConstraintPool[ iContact ];
            btScalar residual = resolveSingleConstraintRowLowerLimitSIMD( m_tmpSolverBodyPool[ solveManifold.m_solverBodyIdA ], m_tmpSolverBodyPool[ solveManifold.m_solverBodyIdB ], solveManifold );
            leastSquaresResidual += residual*residual;
            totalImpulse = solveManifold.m_appliedImpulse;
        }
        if ( totalImpulse > btScalar( 0 ) )
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
                btScalar residual = resolveSingleConstraintRowGenericSIMD( bodyA, bodyB, solveManifold );
                leastSquaresResidual += residual*residual;
            }
        }
    }
    return leastSquaresResidual;
}


void btSequentialImpulseConstraintSolverMt::randomizeConstraintOrdering(int iteration, int numIterations)
{
    // randomize ordering of joint constraints
    int numNonContactPool = m_tmpSolverNonContactConstraintPool.size();
    for ( int j = 0; j < numNonContactPool; ++j )
    {
        int iSwap = btRandInt2( j + 1 );
        m_orderNonContactConstraintPool.swap( j, iSwap );
    }

    //contact/friction constraints are not solved more than numIterations
    if ( iteration < numIterations )
    {
        BatchedConstraints& batchedCons = m_batchedContactConstraints;
        // randomize ordering of phases
        for ( int ii = 0; ii < batchedCons.mPhaseOrder.size(); ++ii )
        {
            int iSwap = btRandInt2( ii + 1 );
            batchedCons.mPhaseOrder.swap( ii, iSwap );
        }

        // for each batch,
        for ( int iBatch = 0; iBatch < batchedCons.mBatches.size(); ++iBatch )
        {
            // randomize ordering of constraints within the batch
            const Range& batch = batchedCons.mBatches[ iBatch ];
            int batchSize = batch.end - batch.begin;
            for ( int iiCons = batch.begin; iiCons < batch.end; ++iiCons )
            {
                int iSwap = batch.begin + btRandInt2( iiCons - batch.begin + 1 );
                batchedCons.mConstraintIndices.swap( iiCons, iSwap );
            }
        }
    }
}


struct NonContactSolverLoop : public btIParallelSumBody
{
    btSequentialImpulseConstraintSolverMt* m_solver;
    const btSequentialImpulseConstraintSolverMt::BatchedConstraints* m_bc;
    int m_iteration;

    NonContactSolverLoop( btSequentialImpulseConstraintSolverMt* solver, const btSequentialImpulseConstraintSolverMt::BatchedConstraints* bc, int iteration )
    {
        m_solver = solver;
        m_bc = bc;
        m_iteration = iteration;
    }
    btScalar sumLoop( int iBegin, int iEnd ) const BT_OVERRIDE
    {
        BT_PROFILE( "NonContactSolverLoop" );
        btScalar sum = 0;
        for ( int iBatch = iBegin; iBatch < iEnd; ++iBatch )
        {
            const btSequentialImpulseConstraintSolverMt::Range& batch = m_bc->mBatches[ iBatch ];
            sum += m_solver->resolveMultipleNonContactConstraints( m_bc->mConstraintIndices, batch.begin, batch.end, m_iteration );
        }
        return sum;
    }
};

btScalar btSequentialImpulseConstraintSolverMt::resolveAllNonContactConstraints(int iteration)
{
    BT_PROFILE( "resolveAllNonContactConstraints" );
    const BatchedConstraints& batchedCons = m_batchedNonContactConstraints;
    NonContactSolverLoop loop( this, &batchedCons, iteration );
    btScalar leastSquaresResidual = 0.f;
    for ( int iiPhase = 0; iiPhase < batchedCons.mPhases.size(); ++iiPhase )
    {
        int iPhase = batchedCons.mPhaseOrder[ iiPhase ];
        const Range& phase = batchedCons.mPhases[ iPhase ];
        int grainSize = 1;
        leastSquaresResidual += btParallelSum( phase.begin, phase.end, grainSize, loop );
    }
    return leastSquaresResidual;
}


struct ContactSolverLoop : public btIParallelSumBody
{
    btSequentialImpulseConstraintSolverMt* m_solver;
    const btSequentialImpulseConstraintSolverMt::BatchedConstraints* m_bc;

    ContactSolverLoop( btSequentialImpulseConstraintSolverMt* solver, const btSequentialImpulseConstraintSolverMt::BatchedConstraints* bc )
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
            const btSequentialImpulseConstraintSolverMt::Range& batch = m_bc->mBatches[ iBatch ];
            sum += m_solver->resolveMultipleContactConstraints( m_bc->mConstraintIndices, batch.begin, batch.end );
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
    for ( int iiPhase = 0; iiPhase < batchedCons.mPhases.size(); ++iiPhase )
    {
        int iPhase = batchedCons.mPhaseOrder[ iiPhase ];
        const Range& phase = batchedCons.mPhases[ iPhase ];
        int grainSize = 1;
        leastSquaresResidual += btParallelSum( phase.begin, phase.end, grainSize, loop );
    }
    return leastSquaresResidual;
}


struct ContactFrictionSolverLoop : public btIParallelSumBody
{
    btSequentialImpulseConstraintSolverMt* m_solver;
    const btSequentialImpulseConstraintSolverMt::BatchedConstraints* m_bc;

    ContactFrictionSolverLoop( btSequentialImpulseConstraintSolverMt* solver, const btSequentialImpulseConstraintSolverMt::BatchedConstraints* bc )
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
            const btSequentialImpulseConstraintSolverMt::Range& batch = m_bc->mBatches[ iBatch ];
            sum += m_solver->resolveMultipleContactFrictionConstraints( m_bc->mConstraintIndices, batch.begin, batch.end );
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
    for ( int iiPhase = 0; iiPhase < batchedCons.mPhases.size(); ++iiPhase )
    {
        int iPhase = batchedCons.mPhaseOrder[ iiPhase ];
        const Range& phase = batchedCons.mPhases[ iPhase ];
        int grainSize = 1;
        leastSquaresResidual += btParallelSum( phase.begin, phase.end, grainSize, loop );
    }
    return leastSquaresResidual;
}


struct InterleavedContactSolverLoop : public btIParallelSumBody
{
    btSequentialImpulseConstraintSolverMt* m_solver;
    const btSequentialImpulseConstraintSolverMt::BatchedConstraints* m_bc;

    InterleavedContactSolverLoop( btSequentialImpulseConstraintSolverMt* solver, const btSequentialImpulseConstraintSolverMt::BatchedConstraints* bc )
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
            const btSequentialImpulseConstraintSolverMt::Range& batch = m_bc->mBatches[ iBatch ];
            sum += m_solver->resolveMultipleContactConstraintsInterleaved( m_bc->mConstraintIndices, batch.begin, batch.end );
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
    for ( int iiPhase = 0; iiPhase < batchedCons.mPhases.size(); ++iiPhase )
    {
        int iPhase = batchedCons.mPhaseOrder[ iiPhase ];
        const Range& phase = batchedCons.mPhases[ iPhase ];
        int grainSize = 1;
        leastSquaresResidual += btParallelSum( phase.begin, phase.end, grainSize, loop );
    }
    return leastSquaresResidual;
}


btScalar btSequentialImpulseConstraintSolverMt::resolveAllRollingFrictionConstraints()
{
    BT_PROFILE( "resolveAllRollingFrictionConstraints" );
    // TODO: batchify if needed (need a test case)
    int numRollingFrictionPoolConstraints = m_tmpSolverContactRollingFrictionConstraintPool.size();
    btScalar leastSquaresResidual = 0.f;
    for ( int j = 0; j<numRollingFrictionPoolConstraints; j++ )
    {
        btSolverConstraint& rollingFrictionConstraint = m_tmpSolverContactRollingFrictionConstraintPool[ j ];
        btScalar totalImpulse = m_tmpSolverContactConstraintPool[ rollingFrictionConstraint.m_frictionIndex ].m_appliedImpulse;
        if ( totalImpulse>btScalar( 0 ) )
        {
            btScalar rollingFrictionMagnitude = rollingFrictionConstraint.m_friction*totalImpulse;
            if ( rollingFrictionMagnitude > rollingFrictionConstraint.m_friction )
                rollingFrictionMagnitude = rollingFrictionConstraint.m_friction;

            rollingFrictionConstraint.m_lowerLimit = -rollingFrictionMagnitude;
            rollingFrictionConstraint.m_upperLimit = rollingFrictionMagnitude;

            btScalar residual = resolveSingleConstraintRowGenericSIMD( m_tmpSolverBodyPool[ rollingFrictionConstraint.m_solverBodyIdA ], m_tmpSolverBodyPool[ rollingFrictionConstraint.m_solverBodyIdB ], rollingFrictionConstraint );
            leastSquaresResidual += residual*residual;
        }
    }
    return leastSquaresResidual;
}



