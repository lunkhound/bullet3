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

#ifndef BT_SEQUENTIAL_IMPULSE_CONSTRAINT_SOLVER_MT_H
#define BT_SEQUENTIAL_IMPULSE_CONSTRAINT_SOLVER_MT_H

#include "btSequentialImpulseConstraintSolver.h"
#include "LinearMath/btThreads.h"


ATTRIBUTE_ALIGNED16(class) btSequentialImpulseConstraintSolverMt : public btSequentialImpulseConstraintSolver
{
public:
	//virtual void convertContacts(btPersistentManifold** manifoldPtr, int numManifolds, const btContactSolverInfo& infoGlobal);
	//virtual void solveGroupCacheFriendlySplitImpulseIterations(btCollisionObject** bodies,int numBodies,btPersistentManifold** manifoldPtr, int numManifolds,btTypedConstraint** constraints,int numConstraints,const btContactSolverInfo& infoGlobal,btIDebugDraw* debugDrawer);
	//virtual btScalar solveGroupCacheFriendlyFinish(btCollisionObject** bodies,int numBodies,const btContactSolverInfo& infoGlobal);
	virtual btScalar solveSingleIteration(int iteration, btCollisionObject** bodies ,int numBodies,btPersistentManifold** manifoldPtr, int numManifolds,btTypedConstraint** constraints,int numConstraints,const btContactSolverInfo& infoGlobal,btIDebugDraw* debugDrawer) BT_OVERRIDE;
	virtual btScalar solveGroupCacheFriendlySetup(btCollisionObject** bodies,int numBodies,btPersistentManifold** manifoldPtr, int numManifolds,btTypedConstraint** constraints,int numConstraints,const btContactSolverInfo& infoGlobal,btIDebugDraw* debugDrawer) BT_OVERRIDE;
	//virtual btScalar solveGroupCacheFriendlyIterations(btCollisionObject** bodies,int numBodies,btPersistentManifold** manifoldPtr, int numManifolds,btTypedConstraint** constraints,int numConstraints,const btContactSolverInfo& infoGlobal,btIDebugDraw* debugDrawer);

    struct Range
    {
        int begin;
        int end;

        Range() : begin(0), end(0) {}
        Range(int _beg, int _end) : begin(_beg), end(_end) {}
    };
    struct BatchedConstraints
    {
        btAlignedObjectArray<int> mConstraintIndices;
        btAlignedObjectArray<Range> mBatches;  // each batch is a range of indices in the mConstraintIndices array
        btAlignedObjectArray<Range> mPhases;  // each phase is range of indices in the mBatches array
        btAlignedObjectArray<int> mPhaseOrder;

        BatchedConstraints() {}
        void setup(btConstraintArray* constraints, const btAlignedObjectArray<btSolverBody>& bodies, int minBatchSize, int maxBatchSize);
        void setup2(btConstraintArray* constraints, const btAlignedObjectArray<btSolverBody>& bodies, int minBatchSize, int maxBatchSize);
    };
    // parameters to control batching
    static bool sAllowNestedParallelForLoops;
    static int sMinimumContactManifoldsForBatching;
    static int sMinBatchSize;
    static int sMaxBatchSize;

protected:
    BatchedConstraints m_batchedNonContactConstraints;
    BatchedConstraints m_batchedContactConstraints;
    int m_numFrictionDirections;
    bool m_useBatching;

    virtual void randomizeConstraintOrdering( int iteration, int numIterations );
    virtual btScalar resolveAllNonContactConstraints( int iteration );
    virtual btScalar resolveAllContactConstraints();
    virtual btScalar resolveAllContactFrictionConstraints();
    virtual btScalar resolveAllContactConstraintsInterleaved();
    virtual btScalar resolveAllRollingFrictionConstraints();

    virtual void setupBatchedConstraints();

public:

	BT_DECLARE_ALIGNED_ALLOCATOR();
	
	btSequentialImpulseConstraintSolverMt();
	virtual ~btSequentialImpulseConstraintSolverMt();

    btScalar resolveMultipleNonContactConstraints( const btAlignedObjectArray<int>& consIndices, int batchBegin, int batchEnd, int iteration );
    btScalar resolveMultipleContactConstraints( const btAlignedObjectArray<int>& consIndices, int batchBegin, int batchEnd );
    btScalar resolveMultipleContactFrictionConstraints( const btAlignedObjectArray<int>& consIndices, int batchBegin, int batchEnd );
    btScalar resolveMultipleContactConstraintsInterleaved( const btAlignedObjectArray<int>& consIndices, int batchBegin, int batchEnd );
	//virtual btScalar solveGroup(btCollisionObject** bodies,int numBodies,btPersistentManifold** manifold,int numManifolds,btTypedConstraint** constraints,int numConstraints,const btContactSolverInfo& info, btIDebugDraw* debugDrawer,btDispatcher* dispatcher);
};




#endif //BT_SEQUENTIAL_IMPULSE_CONSTRAINT_SOLVER_MT_H

