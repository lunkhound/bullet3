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

#include "LinearMath/btQuickprof.h"

#include "BulletCollision/NarrowPhaseCollision/btPersistentManifold.h"

#include "BulletDynamics/ConstraintSolver/btTypedConstraint.h"
#include "BulletDynamics/Dynamics/btRigidBody.h"



bool btSequentialImpulseConstraintSolverMt::s_allowNestedParallelForLoops = false;  // some task schedulers don't like nested loops
int btSequentialImpulseConstraintSolverMt::s_minimumContactManifoldsForBatching = 1500;
int btSequentialImpulseConstraintSolverMt::s_minBatchSize = 80;
int btSequentialImpulseConstraintSolverMt::s_maxBatchSize = 100;


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
    m_batchedContactConstraints.setup( &m_tmpSolverContactConstraintPool, m_tmpSolverBodyPool, s_minBatchSize, s_maxBatchSize, &m_scratchMemory, m_createBatchesWorkArray );
}


void btSequentialImpulseConstraintSolverMt::setupBatchedJointConstraints()
{
    BT_PROFILE("setupBatchedJointConstraints");
    m_batchedJointConstraints.setup( &m_tmpSolverNonContactConstraintPool, m_tmpSolverBodyPool, s_minBatchSize, s_maxBatchSize, &m_scratchMemory, m_createBatchesWorkArray );
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
    const btBatchedConstraints* m_bc;
    const btContactSolverInfo* m_infoGlobal;

    SetupContactConstraintsLoop( btSequentialImpulseConstraintSolverMt* solver, const btBatchedConstraints* bc, const btContactSolverInfo& infoGlobal )
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
            const btBatchedConstraints::Range& batch = m_bc->m_batches[ iBatch ];
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
        const btBatchedConstraints& batchedCons = m_batchedContactConstraints;
        SetupContactConstraintsLoop loop( this, &batchedCons, infoGlobal );
        for ( int iiPhase = 0; iiPhase < batchedCons.m_phases.size(); ++iiPhase )
        {
            int iPhase = batchedCons.m_phaseOrder[ iiPhase ];
            const btBatchedConstraints::Range& phase = batchedCons.m_phases[ iPhase ];
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
    if (false)
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
    if ( numManifolds >= s_minimumContactManifoldsForBatching &&
        (s_allowNestedParallelForLoops || !btThreadsAreRunning())
        )
    {
        m_useBatching = true;
        m_batchedContactConstraints.m_debugDrawer = debugDrawer;
        m_batchedJointConstraints.m_debugDrawer = debugDrawer;
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
    const btBatchedConstraints* m_bc;

    ContactSplitPenetrationImpulseSolverLoop( btSequentialImpulseConstraintSolverMt* solver, const btBatchedConstraints* bc )
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
            const btBatchedConstraints::Range& batch = m_bc->m_batches[ iBatch ];
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
                const btBatchedConstraints& batchedCons = m_batchedContactConstraints;
                ContactSplitPenetrationImpulseSolverLoop loop( this, &batchedCons );
                btScalar leastSquaresResidual = 0.f;
                for ( int iiPhase = 0; iiPhase < batchedCons.m_phases.size(); ++iiPhase )
                {
                    int iPhase = batchedCons.m_phaseOrder[ iiPhase ];
                    const btBatchedConstraints::Range& phase = batchedCons.m_phases[ iPhase ];
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


void btSequentialImpulseConstraintSolverMt::randomizeBatchedConstraintOrdering( btBatchedConstraints* batchedConstraints )
{
    btBatchedConstraints& bc = *batchedConstraints;
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
        const btBatchedConstraints::Range& batch = bc.m_batches[ iBatch ];
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
    const btBatchedConstraints* m_bc;
    int m_iteration;

    JointSolverLoop( btSequentialImpulseConstraintSolverMt* solver, const btBatchedConstraints* bc, int iteration )
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
            const btBatchedConstraints::Range& batch = m_bc->m_batches[ iBatch ];
            sum += m_solver->resolveMultipleJointConstraints( m_bc->m_constraintIndices, batch.begin, batch.end, m_iteration );
        }
        return sum;
    }
};


btScalar btSequentialImpulseConstraintSolverMt::resolveAllJointConstraints(int iteration)
{
    BT_PROFILE( "resolveAllJointConstraints" );
    const btBatchedConstraints& batchedCons = m_batchedJointConstraints;
    JointSolverLoop loop( this, &batchedCons, iteration );
    btScalar leastSquaresResidual = 0.f;
    for ( int iiPhase = 0; iiPhase < batchedCons.m_phases.size(); ++iiPhase )
    {
        int iPhase = batchedCons.m_phaseOrder[ iiPhase ];
        const btBatchedConstraints::Range& phase = batchedCons.m_phases[ iPhase ];
        int grainSize = 1;
        leastSquaresResidual += btParallelSum( phase.begin, phase.end, grainSize, loop );
    }
    return leastSquaresResidual;
}


struct ContactSolverLoop : public btIParallelSumBody
{
    btSequentialImpulseConstraintSolverMt* m_solver;
    const btBatchedConstraints* m_bc;

    ContactSolverLoop( btSequentialImpulseConstraintSolverMt* solver, const btBatchedConstraints* bc )
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
            const btBatchedConstraints::Range& batch = m_bc->m_batches[ iBatch ];
            sum += m_solver->resolveMultipleContactConstraints( m_bc->m_constraintIndices, batch.begin, batch.end );
        }
        return sum;
    }
};


btScalar btSequentialImpulseConstraintSolverMt::resolveAllContactConstraints()
{
    BT_PROFILE( "resolveAllContactConstraints" );
    const btBatchedConstraints& batchedCons = m_batchedContactConstraints;
    ContactSolverLoop loop( this, &batchedCons );
    btScalar leastSquaresResidual = 0.f;
    for ( int iiPhase = 0; iiPhase < batchedCons.m_phases.size(); ++iiPhase )
    {
        int iPhase = batchedCons.m_phaseOrder[ iiPhase ];
        const btBatchedConstraints::Range& phase = batchedCons.m_phases[ iPhase ];
        int grainSize = 1;
        leastSquaresResidual += btParallelSum( phase.begin, phase.end, grainSize, loop );
    }
    return leastSquaresResidual;
}


struct ContactFrictionSolverLoop : public btIParallelSumBody
{
    btSequentialImpulseConstraintSolverMt* m_solver;
    const btBatchedConstraints* m_bc;

    ContactFrictionSolverLoop( btSequentialImpulseConstraintSolverMt* solver, const btBatchedConstraints* bc )
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
            const btBatchedConstraints::Range& batch = m_bc->m_batches[ iBatch ];
            sum += m_solver->resolveMultipleContactFrictionConstraints( m_bc->m_constraintIndices, batch.begin, batch.end );
        }
        return sum;
    }
};


btScalar btSequentialImpulseConstraintSolverMt::resolveAllContactFrictionConstraints()
{
    BT_PROFILE( "resolveAllContactFrictionConstraints" );
    const btBatchedConstraints& batchedCons = m_batchedContactConstraints;
    ContactFrictionSolverLoop loop( this, &batchedCons );
    btScalar leastSquaresResidual = 0.f;
    for ( int iiPhase = 0; iiPhase < batchedCons.m_phases.size(); ++iiPhase )
    {
        int iPhase = batchedCons.m_phaseOrder[ iiPhase ];
        const btBatchedConstraints::Range& phase = batchedCons.m_phases[ iPhase ];
        int grainSize = 1;
        leastSquaresResidual += btParallelSum( phase.begin, phase.end, grainSize, loop );
    }
    return leastSquaresResidual;
}


struct InterleavedContactSolverLoop : public btIParallelSumBody
{
    btSequentialImpulseConstraintSolverMt* m_solver;
    const btBatchedConstraints* m_bc;

    InterleavedContactSolverLoop( btSequentialImpulseConstraintSolverMt* solver, const btBatchedConstraints* bc )
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
            const btBatchedConstraints::Range& batch = m_bc->m_batches[ iBatch ];
            sum += m_solver->resolveMultipleContactConstraintsInterleaved( m_bc->m_constraintIndices, batch.begin, batch.end );
        }
        return sum;
    }
};


btScalar btSequentialImpulseConstraintSolverMt::resolveAllContactConstraintsInterleaved()
{
    BT_PROFILE( "resolveAllContactConstraintsInterleaved" );
    const btBatchedConstraints& batchedCons = m_batchedContactConstraints;
    InterleavedContactSolverLoop loop( this, &batchedCons );
    btScalar leastSquaresResidual = 0.f;
    for ( int iiPhase = 0; iiPhase < batchedCons.m_phases.size(); ++iiPhase )
    {
        int iPhase = batchedCons.m_phaseOrder[ iiPhase ];
        const btBatchedConstraints::Range& phase = batchedCons.m_phases[ iPhase ];
        int grainSize = 1;
        leastSquaresResidual += btParallelSum( phase.begin, phase.end, grainSize, loop );
    }
    return leastSquaresResidual;
}


struct ContactRollingFrictionSolverLoop : public btIParallelSumBody
{
    btSequentialImpulseConstraintSolverMt* m_solver;
    const btBatchedConstraints* m_bc;

    ContactRollingFrictionSolverLoop( btSequentialImpulseConstraintSolverMt* solver, const btBatchedConstraints* bc )
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
            const btBatchedConstraints::Range& batch = m_bc->m_batches[ iBatch ];
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
        const btBatchedConstraints& batchedCons = m_batchedContactConstraints;
        ContactRollingFrictionSolverLoop loop( this, &batchedCons );
        btScalar leastSquaresResidual = 0.f;
        for ( int iiPhase = 0; iiPhase < batchedCons.m_phases.size(); ++iiPhase )
        {
            int iPhase = batchedCons.m_phaseOrder[ iiPhase ];
            const btBatchedConstraints::Range& phase = batchedCons.m_phases[ iPhase ];
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


