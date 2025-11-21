// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2014 projectchrono.org
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Authors: Mariusz Warzecha
// =============================================================================
//
// Simple contact model between a sphere and a mesh surface.
// It is intended for tests of contact detection between Chrono native contact
// geometries and triangles building mesh.
// For now it uses Hook contact force model for simplicity.
// =============================================================================

#include "chrono/collision/ChCollisionSystem.h"
#include "chrono/collision/bullet/ChCollisionSystemBullet.h"
#include "chrono_irrlicht/ChVisualSystemIrrlicht.h"
#include "chrono/physics/ChSystemSMC.h"
#include "chrono/physics/ChBodyEasy.h"
#include "chrono/utils/ChUtilsGeometry.h"
#include "chrono/utils/ChUtilsGenerators.h"
#include "chrono/utils/ChUtilsInputOutput.h"
#include "chrono/assets/ChVisualShapeTriangleMesh.h"
#include "chrono/assets/ChTexture.h"
#include "chrono/geometry/ChTriangleMeshConnected.h"
#include "chrono/physics/ChContactContainerSMC.h"
#include "chrono_dfc/ChCollisionSystemBulletSingleTriangle.h"


// Use the main namespace of Chrono, and other chrono namespaces
using namespace chrono;
using namespace chrono::irrlicht;

class ContactContainerPrint : public ChContactContainerSMC {
public:
    ContactContainerPrint() {}
    ~ContactContainerPrint() {}
    void PrintContacts() {
        auto iter = contactlist_6_6.begin();
        while (iter != contactlist_6_6.end()) {
            ChVector3d p1 = (*iter)->GetContactP1();
            ChVector3d p2 = (*iter)->GetContactP2();
            ChVector3d force = (*iter)->GetContactForce();
            double contact_penetration = (*iter)->GetContactPenetration();
            ChVector3d contact_normal = (*iter)->GetContactNormal();
            std::cout << "Coordinates of contact point P1, in absolute coordinates: " << p1 << "\n";
            std::cout << "Coordinates of contact point P2, in absolute coordinates: " << p2 << "\n";
            std::cout << "Contact force, in contact coordinate system: " << force << "\n";
            std::cout << "Contact penetration, positive if there is overlap: " << contact_penetration << "\n";
            std::cout << "Contact normal, in absolute coordinates: " << contact_normal << "\n\n";
            ++iter;
        }
    }
};

int main(int argc, char* argv[]) {
    std::cout << "Copyright (c) 2017 projectchrono.org\nChrono version: " << CHRONO_VERSION << "\n\n";
	SetChronoDataPath(CHRONO_DATA_DIR);
    std::string mesh_dir = std::string(CHRONO_DATA_DIR) + "models"; 

// Create a Chrono system
    ChSystemSMC sys;
	sys.SetGravitationalAcceleration(ChVector3d(0, 0, 0));
	sys.SetNumThreads(1,1,1);
    auto single_triangle_detection = chrono_types::make_shared<ChCollisionSystemBulletSingleTriangle>();
    sys.SetCollisionSystem(single_triangle_detection);  //  inherited from BULLET,

// Create mesh and sphere
    auto floor_mat = chrono_types::make_shared<ChContactMaterialSMC>();
	floor_mat->SetYoungModulus(1e2);
    floor_mat->SetFriction(0.5f);
    floor_mat->SetRestitution(0.0f);
    floor_mat->SetAdhesion(0);
    auto meshed_surface = chrono_types::make_shared<ChBodyEasyMesh>(
        (mesh_dir+"/plane_xy_30mm_grid1mm.obj").c_str(),  // file name for OBJ Wavefront mesh
        10,                                         // density of the body
        false,                                         // automatic evaluation of mass, COG position, inertia tensor
        true,                                         // attach visualization asset
        true,                                         // enable the collision detection
        floor_mat,                                     // surface contact material
        0  // radius of 'inflating' of mesh (for more robust collision detection)
        );	
	meshed_surface->SetPos(ChVector3d(0,0,0));
    meshed_surface->SetFixed(true);
    sys.AddBody(meshed_surface);
    double radius = 5;
    double mass = 0.001*4/3*pow(radius, 3)*CH_PI;
    auto sphere = chrono_types::make_shared<ChBody>();
    sphere->SetInertiaXX((2.0 / 5.0) * mass * pow(radius, 2) * ChVector3d(1, 1, 1));
    sphere->SetMass(mass);
    sphere->SetPos(ChVector3d(0, 0, 4.9));
    sphere->SetPosDt(ChVector3d(0, 0, 0));
	auto sphere_shape = chrono_types::make_shared<ChCollisionShapeSphere>(floor_mat, radius);
    sphere->AddCollisionShape(sphere_shape);    
    sphere->EnableCollision(true);
    sys.AddBody(sphere);
    sys.SetContactForceModel(ChSystemSMC::ContactForceModel::Hooke);
    sys.SetAdhesionForceModel(ChSystemSMC::AdhesionForceModel::Constant);
    auto container = chrono_types::make_shared<ContactContainerPrint>();
    sys.SetContactContainer(container);	
    while (sys.GetChTime() < 1.0E-2) {		
        sys.DoStepDynamics(2.0E-2);
        ChVector3d sphere_pos = sphere->GetPos();
        ChVector3d sphere_lin_vel = sphere->GetPosDt();
        ChVector3d sphere_ang_vel = sphere->GetAngVelLocal();
        ChVector3d sphere_resultant_force = sphere->GetAppliedForce();
        std::cout << "time: " << sys.GetChTime() << "\n";
        std::cout << "  Sphere position:         " << sphere_pos << "\n";
        std::cout << "  Sphere linear velocity:  " << sphere_lin_vel << "\n";
        std::cout << "  Sphere angular velocity: " << sphere_ang_vel << "\n";
        std::cout << "  Resultant force acting on sphere: " << sphere_resultant_force << "\n\n";
        container->PrintContacts();
    }
    return 0;
}
