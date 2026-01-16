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
// Reference simulation for GPU model. It simulates layers of ordered DFC particles
// in a container. Number of layers, particles and its initial velocity may vary.
// The intention is to have a reference solution to exactly the same model run on
// GPU.
// =============================================================================

#include <string>
#include <cstring>
#include <filesystem>
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
#include "chrono_dfc/FreshConcreteContact.cpp"


// Use the main namespace of Chrono, and other chrono namespaces
using namespace chrono;
using namespace chrono::irrlicht;
using namespace std::filesystem;

// Argument parser (run with or without visualization
struct AppOptions {
    bool visualization = false;     // default: headless
    int  fps_limit = 60;            // opcjonalnie
};

AppOptions parse_args(int argc, char* argv[]) {
    AppOptions opt;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--vis" || a == "--visualization") {
            opt.visualization = true;
        } else if (a == "--no-vis") {
            opt.visualization = false;
        } else if (a.rfind("--fps=", 0) == 0) {
            opt.fps_limit = std::max(1, std::stoi(a.substr(6)));
        }
    }
    return opt;
}

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
            ChContactable* objA = (*iter)->GetObjA();
            ChContactable* objB = (*iter)->GetObjB();
            std::cout << "Coordinates of contact point P1, in absolute coordinates: " << p1 << "\n";
            std::cout << "Coordinates of contact point P2, in absolute coordinates: " << p2 << "\n";
            std::cout << "Contact force, in contact coordinate system: " << force << "\n";
            std::cout << "Contact penetration, positive if there is overlap: " << contact_penetration << "\n";
            std::cout << "Contact normal, in absolute coordinates: " << contact_normal << "\n";
            std::cout << "Contactable objects A and B: " << objA << "  " << objB << "\n\n";
            ++iter;
        }
    }
};

std::vector<std::shared_ptr<ChBody>> create_container(ChSystemSMC& sys,
                                                      const std::shared_ptr<ChMaterialFCM>& DFC_material,
                                                      int size_in_particles,
                                                      float sphere_radius) {
    // container size given in mm, wall thickness is fixed in function code
    // box is positioned is such way, that first particle in first layer has coordinates (0, 0, 0)
    // that ensures the same position of particles as in GPU simulation
    
    // dummy material, only to create the bodies, it is not used in contact force calculation
    float density = 1.0;  // not relevant as container is fixed
    float wall_thickness = 10;  // assumed constant, does not matter if size of container will change
    std::shared_ptr<ChContactMaterialSMC> container_wall_material = chrono_types::make_shared<ChContactMaterialSMC>();
    container_wall_material->SetYoungModulus(1e2);
    container_wall_material->SetFriction(0.5);
    container_wall_material->SetRestitution(0.2);
    float h = DFC_material->Get_mortar_h();
    float size = (size_in_particles - 1) * (2 * sphere_radius - h) + 2*sphere_radius;

    std::shared_ptr<ChBodyEasyBox> bottom_wall = chrono_types::make_shared<ChBodyEasyBox>(size, size, wall_thickness, // x, y, z size
                                                                                          density, // density
                                                                                          true, // visualise
                                                                                          true, // collide
                                                                                          container_wall_material); // material
    bottom_wall->SetPos(ChVector3d(size/2 - sphere_radius, size/2 - sphere_radius, -wall_thickness/2));
    bottom_wall->SetFixed(true);
    bottom_wall->GetVisualShape(0)->SetOpacity(0.25);
    bottom_wall->EnableCollision(true);
    sys.Add(bottom_wall);

    std::shared_ptr<ChBodyEasyBox> side_wall_1 = chrono_types::make_shared<ChBodyEasyBox>(size, wall_thickness, size, // x, y, z size
                                                                                          density, // density
                                                                                          true, // visualise
                                                                                          true, // collide
                                                                                          container_wall_material); // material
    side_wall_1->SetPos(ChVector3d(size/2 - sphere_radius, size - sphere_radius + wall_thickness/2, size/2));
    side_wall_1->SetFixed(true);
    side_wall_1->GetVisualShape(0)->SetOpacity(0.25);
    side_wall_1->EnableCollision(true);
    sys.Add(side_wall_1);

    std::shared_ptr<ChBodyEasyBox> side_wall_2 = chrono_types::make_shared<ChBodyEasyBox>(size, wall_thickness, size, // x, y, z size
                                                                                          density, // density
                                                                                          true, // visualise
                                                                                          true, // collide
                                                                                          container_wall_material); // material
    side_wall_2->SetPos(ChVector3d(size/2 - sphere_radius, -(sphere_radius + wall_thickness/2), size/2));
    side_wall_2->SetFixed(true);
    side_wall_2->GetVisualShape(0)->SetOpacity(0.25);
    side_wall_2->EnableCollision(true);
    sys.Add(side_wall_2);

    std::shared_ptr<ChBodyEasyBox> side_wall_3 = chrono_types::make_shared<ChBodyEasyBox>(wall_thickness, size, size, // x, y, z size
                                                                                          density, // density
                                                                                          true, // visualise
                                                                                          true, // collide
                                                                                          container_wall_material); // material
    side_wall_3->SetPos(ChVector3d(-(sphere_radius + wall_thickness/2), size/2 - sphere_radius, size/2));
    side_wall_3->SetFixed(true);
    side_wall_3->GetVisualShape(0)->SetOpacity(0.25);
    side_wall_3->EnableCollision(true);
    sys.Add(side_wall_3);

    std::shared_ptr<ChBodyEasyBox> side_wall_4 = chrono_types::make_shared<ChBodyEasyBox>(wall_thickness, size, size, // x, y, z size
                                                                                          density, // density
                                                                                          true, // visualise
                                                                                          true, // collide
                                                                                          container_wall_material); // material
    side_wall_4->SetPos(ChVector3d((size - sphere_radius + wall_thickness/2), size/2 - sphere_radius, size/2));
    side_wall_4->SetFixed(true);
    side_wall_4->GetVisualShape(0)->SetOpacity(0.25);
    side_wall_4->EnableCollision(true);
    sys.Add(side_wall_4);

    std::vector<std::shared_ptr<ChBody>> walls;
    walls.push_back(bottom_wall);
    walls.push_back(side_wall_1);
    walls.push_back(side_wall_2);
    walls.push_back(side_wall_3);
    walls.push_back(side_wall_4);
    return walls;
}

std::shared_ptr<ChMaterialFCM> create_concrete_material() {
    const float E_Nm = 2.0;
    const float E_Na = 100.0;
    const float h = 3.0;
    const float sigma_tau0 = 1E-7;
    const float alpha_a = 0.25;
    const float beta = 0.5;
    const float n = 1.0;
    const float sigma_tau = 2.0E-4;
    const float kappa_0 = 100.0;
    const float eta_inf = 2.0E-5;
    std::shared_ptr<ChMaterialFCM> DFC_concrete = chrono_types::make_shared<ChMaterialFCM>(E_Nm, E_Na, h, alpha_a, beta, n, sigma_tau, sigma_tau0, kappa_0, eta_inf);
    DFC_concrete->Set_flocbeta(0.01);
    DFC_concrete->Set_flocm(1.0);
    DFC_concrete->Set_flocTcr(100.);
    DFC_concrete->Set_lambda_init(2.0);
    DFC_concrete->Set_ThixOnFlag(false);
    return DFC_concrete;
}

std::shared_ptr<ChBody> create_sphere(ChSystemSMC& sys, const std::shared_ptr<ChMaterialFCM>& material, const ChVector3d& pos, const ChVector3d vel, float radius, float density){
    // dummy material, only to create the bodies, it is not used in contact force calculation
    // not 100 % sure
    std::shared_ptr<ChContactMaterialSMC> dummy_sphere_material = chrono_types::make_shared<ChContactMaterialSMC>();
    dummy_sphere_material->SetYoungModulus(1e2);
    dummy_sphere_material->SetFriction(0.5);
    dummy_sphere_material->SetRestitution(0.2);
    dummy_sphere_material->SetAdhesion(0.0);
    float mass = density*4/3*pow(radius, 3)*CH_PI;
    std::shared_ptr<ChBody> sphere = chrono_types::make_shared<ChBody>();
    sphere->SetInertiaXX((2.0 / 5.0) * mass * pow(radius, 2) * ChVector3d(1, 1, 1));
    sphere->SetMass(mass);
    sphere->SetPos(pos);
    sphere->SetPosDt(vel);
    std::shared_ptr<ChCollisionShapeSphere> sphere_shape = chrono_types::make_shared<ChCollisionShapeSphere>(dummy_sphere_material, radius);
    sphere->AddCollisionShape(sphere_shape);    
    sphere->EnableCollision(true);
    std::shared_ptr<ChVisualShapeSphere> sphere_mor = chrono_types::make_shared<ChVisualShapeSphere>(radius);
	sphere_mor->SetColor(ChColor(128.f/255, 128.f/255, 128.f/255));
	sphere_mor->SetOpacity(0.25f);
	sphere->AddVisualShape(sphere_mor);
    float h = material->Get_mortar_h();
    std::shared_ptr<ChVisualShapeSphere> sphere_agg = chrono_types::make_shared<ChVisualShapeSphere>(radius - h);
	sphere_agg->SetColor(ChColor(5.f/255, 48.f/255, 173.f/255));	
	sphere->AddVisualShape(sphere_agg);
    sys.AddBody(sphere);
    return sphere;
}

std::vector<std::shared_ptr<ChBody>> create_sphere_layers(ChSystemSMC& sys,
                                                          const std::shared_ptr<ChMaterialFCM>& material,
                                                          const ChVector3d vel,
                                                          float radius,
                                                          float density,
                                                          int layer_size_in_spheres){
    std::vector<std::shared_ptr<ChBody>> created_spheres;
    float h = material->Get_mortar_h();
    for (int i = 0; i < layer_size_in_spheres; ++i) {
        float layer_shift = (i % 2 == 0) ? 0 : radius / 2;
        for (int j = 0; j < layer_size_in_spheres; ++j) {
            for (int k = 0; k < layer_size_in_spheres; ++k) {
                float x = k * (2 * radius - h) + layer_shift;
                float y = j * (2 * radius - h) + layer_shift;
                float z = i * (2 * radius - h) + 2*radius;
                std::shared_ptr<ChBody> temp_sphere = create_sphere(sys, material, ChVector3d(x, y, z), vel, radius, density);
                created_spheres.push_back(temp_sphere);
            }
        }
    }
    return created_spheres;
}

void write_sphere_file(const std::string& file_name, std::vector<std::shared_ptr<ChBody>> bodies_to_write){
    std::ofstream csv_file(file_name);
    csv_file << "X,Y,Z,r,absv\n";
    for (std::shared_ptr<ChBody> body : bodies_to_write) {
        ChVector3f pos = body->GetPos();
        ChVector3f vel = body->GetPosDt();
        float abs_vel = vel.Length();
        std::shared_ptr<ChCollisionShape> contact_shape = body->GetCollisionModel()->GetShapeInstance(0).first;
        std::shared_ptr<ChCollisionShapeSphere> contact_sphere = std::static_pointer_cast<ChCollisionShapeSphere>(contact_shape);
        csv_file << pos.x() << "," << pos.y() << "," << pos.z() << "," << contact_sphere-> GetRadius() << "," << abs_vel << "\n";
    }
}

int main(int argc, char* argv[]) {
    std::cout << "Application linked against Chrono version: " << CHRONO_VERSION << "\n\n";
	SetChronoDataPath(CHRONO_DATA_DIR);
    // parse arguments
    const auto opt = parse_args(argc, argv);
    
    float simulation_time = 1.0;
    float sphere_radius = 5;
    int number_of_paricles_in_layer = 5;
    float initial_velocity = -50.0;
    path out_dir = current_path();
    out_dir += "/experiment_test_CPU";
    out_dir += "/simulation_time_" + std::to_string(int(simulation_time)) + "s";
    out_dir += "/box_size_in_spheres_" + std::to_string(number_of_paricles_in_layer);
    out_dir += "/initial_velocity_" + std::to_string(int(initial_velocity));
    std::cout << out_dir;
    create_directories(out_dir);

// Create a Chrono system
    ChSystemSMC sys;
	sys.SetGravitationalAcceleration(ChVector3d(0, 0, 0));
	sys.SetNumThreads(1,1,1);
    sys.SetCollisionSystemType(ChCollisionSystem::Type::BULLET);
    auto Contact_force_algorithm = chrono_types::make_unique<FCContactForce>();
	auto materialFCM = create_concrete_material();    
    Contact_force_algorithm->Set_Material(materialFCM);
    sys.SetContactForceTorqueAlgorithm(std::move(Contact_force_algorithm));
    //
    auto container = chrono_types::make_shared<MyContactContainer>();
    sys.SetContactContainer(container);
    std::vector<std::shared_ptr<ChBody>> container_walls = create_container(sys, materialFCM, number_of_paricles_in_layer, sphere_radius);
    std::shared_ptr<ChVisualSystemIrrlicht>  vis = nullptr;
    if (opt.visualization) {
        vis = chrono_types::make_shared<ChVisualSystemIrrlicht>();
        vis->SetWindowSize(800, 600);
        vis->SetWindowTitle("Simulation of 3D particle chains for comparison with GPU");
        vis->Initialize();
        vis->AddLogo();
        vis->AddSkyBox();
        vis->AddCamera(ChVector3d(250, 250, 500));
        vis->AddTypicalLights();
    }
    std::vector<std::shared_ptr<ChBody>> my_spheres = create_sphere_layers(sys, materialFCM, ChVector3d(0, 0, initial_velocity), sphere_radius, 1.4858e-9, number_of_paricles_in_layer);
    bool closed_window = false;
    if (opt.visualization) {
        vis->AttachSystem(&sys);
    }
    float time_step = 2.0E-5;
    int number_of_files_to_save = 100;
    int save_dt_in_steps = static_cast<int>(simulation_time / time_step / number_of_files_to_save);
    int simulation_steps = 0;
    int total_save_counter = 0;
    
    while (sys.GetChTime() < simulation_time && !closed_window) {		

        if (simulation_steps % save_dt_in_steps == 0) {
            char file_name[200];
            ++total_save_counter;
            sprintf(file_name, "%s/Case_%04d.csv", out_dir.c_str(), total_save_counter);
            write_sphere_file(std::string(file_name), my_spheres);
            std::cout << "Current simulation time: " << sys.GetChTime() << " " << total_save_counter-1  << " % of simulation has been completed.\n";
        }
        
        sys.DoStepDynamics(time_step);
        if (opt.visualization) {
            vis->BeginScene();
            vis->Render();
            vis->EndScene();
            if (!vis->Run()) {
                closed_window = true;
            }
        }
        simulation_steps += 1;
    }
    return 0;
}

