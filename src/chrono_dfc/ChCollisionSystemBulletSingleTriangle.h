// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Authors: Mariusz Warzecha
// =============================================================================
//
// Custom collision system based on Bullet library.
// Reports contact only for one triangle in sphere - mesh contact scenario.
// =============================================================================

// Use the main namespace of Chrono, and other chrono namespaces
using namespace chrono;

class ChCollisionSystemBulletSingleTriangle : public ChCollisionSystemBullet {
public:
    ChCollisionSystemBulletSingleTriangle() {}
    ~ChCollisionSystemBulletSingleTriangle() {}
    void PostProcess() override {
        int num_manifolds = bt_collision_world->getDispatcher()->getNumManifolds();
        int num_contacts = 0;
        int num_contacts_compare = 0;
        int triangle_index = 0;
        double max_penetration = 0;
        for (int i = 0; i < num_manifolds; ++i) {
            cbtPersistentManifold* contact_manifold = bt_collision_world->getDispatcher()->getManifoldByIndexInternal(i);
            num_contacts = contact_manifold->getNumContacts();
            if (num_contacts != 0) {
                const cbtCollisionObject* obA = contact_manifold->getBody0();
                const cbtCollisionObject* obB = contact_manifold->getBody1();
                cbtManifoldPoint& contact_point = contact_manifold->getContactPoint(0);
                double max_penetration = std::abs(contact_point.getDistance());
                for (int j = i + 1; j < num_manifolds; ++j) {
                    cbtPersistentManifold* contact_manifold_compare = bt_collision_world->getDispatcher()->getManifoldByIndexInternal(j);
                    num_contacts_compare = contact_manifold_compare->getNumContacts();
                    if (num_contacts_compare != 0) {
                        const cbtCollisionObject* obA_compare = contact_manifold_compare->getBody0();
                        const cbtCollisionObject* obB_compare = contact_manifold_compare->getBody1();
                        if (obA == obA_compare && obB == obB_compare) {
                            //assume only one contact in manifold (expected for a triangle)
                            contact_point = contact_manifold_compare->getContactPoint(0);
                            double penetration = std::abs(contact_point.getDistance());
                            if (penetration > max_penetration) {
                                bt_collision_world->getDispatcher()->clearManifold(bt_collision_world->getDispatcher()->getManifoldByIndexInternal(triangle_index));
                                triangle_index = j;
                                max_penetration = penetration;
                            } else {
                                bt_collision_world->getDispatcher()->clearManifold(contact_manifold);
                            }
                        }
                    }
                }
            }
        }
    }    
};
