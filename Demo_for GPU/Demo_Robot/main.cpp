// Robot G-code execution

#include "chrono/physics/ChSystemNSC.h"
#include "chrono/physics/ChContactMaterialNSC.h"
#include "chrono/physics/ChBodyEasy.h"
#include "chrono/physics/ChLinkMate.h"
#include "chrono/physics/ChLinkMotorRotationAngle.h"
#include "chrono/core/ChRealtimeStep.h"
#include "chrono_irrlicht/ChVisualSystemIrrlicht.h"
#include "chrono/serialization/ChArchiveJSON.h"
#include "chrono/physics/ChLinkMotionImposed.h"
#include "chrono/functions/ChFunctionPositionLine.h"
#include "chrono/functions/ChFunctionRotationBSpline.h"

#include "chrono/physics/ChLinkMotorRotationSpeed.h"

#include "chrono/physics/ChParticleCloud.h"
#include "chrono/particlefactory/ChParticleEmitter.h"
#include "chrono/particlefactory/ChParticleRemover.h"
#include "chrono/geometry/ChBox.h"
#include "chrono/collision/ChCollisionModel.h"


#include "chrono/timestepper/ChAssemblyAnalysis.h"

#include "ChParserGCode.h"
#include "chrono_thirdparty/rapidjson/prettywriter.h"
#include "chrono_thirdparty/rapidjson/stringbuffer.h"
#include "chrono_thirdparty/filesystem/path.h"

#include <filesystem>

using namespace chrono;
using namespace chrono::irrlicht;
using namespace chrono::particlefactory;



std::shared_ptr<ChLinkMotorRotationAngle> CreateMotor(std::shared_ptr<ChBody> body1, std::shared_ptr<ChBody> body2, ChFramed frame) {
	auto motor = chrono_types::make_shared<ChLinkMotorRotationAngle>();
	motor->Initialize(body1, body2, frame);

	auto motfun = chrono_types::make_shared<ChFunctionSetpoint>();
	motor->SetMotorFunction(motfun);

	return motor;
}

std::shared_ptr<ChLinkMateRevolute> CreateRevolute(std::shared_ptr<ChBody> body1, std::shared_ptr<ChBody> body2, ChFramed frame) {
	auto revolute = chrono_types::make_shared<ChLinkMateRevolute>();
	revolute->Initialize(body1, body2, frame);

	return revolute;
}



int main(int argc, char* argv[]) {
	SetChronoDataPath(CHRONO_DATA_DIR);


	//////old
	filesystem::path cwd1 = filesystem::path::getcwd();
    std::cout << "Current working directory: " << cwd1 << std::endl;
	
	/////new
	/// Get the current diectory
    ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////////
	/*
    std::string current_dir(argv[0]);
    int pos = current_dir.find_last_of("/\\");
    current_dir=current_dir.substr(0, pos-5); 
	*/

	//std::string gcode_file = "bin/data/gcode.txt";
	
	std::filesystem::path cwd = std::filesystem::current_path();
	std::filesystem::path gcode_file = cwd / "data/gcode_Buildability.txt";
	std::string gcode_file_str = gcode_file.string(); // Convert to std::string if needed

	//std::string gcode_file = cwd + "/data/gcode.txt";
	if (argc > 1) {
		gcode_file = argv[1];
	}

	// System
	ChSystemNSC sys;
	//sys.SetCollisionSystemType(ChCollisionSystem::Type::BULLET);

	sys.SetCollisionSystemType(chrono::ChCollisionSystem::Type::BULLET);
	ChCollisionModel::SetDefaultSuggestedEnvelope(0.001);
	sys.SetGravitationalAcceleration(ChVector3d(0, -9.81, 0));


	// Import CAD
	try {
	    //std::ifstream filei("bin/data/robot_export.json");
		std::ifstream filei( cwd / "data/Robot_export_j2.json");
	    ChArchiveInJSON archivein(filei);
	    archivein >> CHNVP(sys);
        filesystem::path shapes_path = "Robot_export_j2_shapes";


        bool filepathExists = shapes_path.exists() && shapes_path.is_directory();
        if (!filepathExists) {
            std::cerr << "CAD shapes cannot be found. Must be placed in: " << cwd << "\\" << shapes_path << std::endl;
			//std::cerr << "CAD shapes cannot be found. Must be placed in: " << pos << "\\" << shapes_path << std::endl;
            return 1;
        }

	    sys.RemoveAllLinks();
	}
	catch (const std::exception& exc) {
	    std::cerr << "Unable to load CAD." << std::endl;
	    return 1;
	}
	sys.ShowHierarchy(std::cout);

	// Robot
	//auto base = std::dynamic_pointer_cast<ChBodyAuxRef>(sys.SearchBody("base-1"));
	//auto shoulder = std::dynamic_pointer_cast<ChBodyAuxRef>(sys.SearchBody("shoulder-1"));
	//auto biceps = std::dynamic_pointer_cast<ChBodyAuxRef>(sys.SearchBody("biceps-1"));
	//auto elbow = std::dynamic_pointer_cast<ChBodyAuxRef>(sys.SearchBody("elbow-1"));
	//auto forearm = std::dynamic_pointer_cast<ChBodyAuxRef>(sys.SearchBody("forearm-1"));
	//auto wrist = std::dynamic_pointer_cast<ChBodyAuxRef>(sys.SearchBody("wrist-1"));
	//auto end_effector = std::dynamic_pointer_cast<ChBodyAuxRef>(sys.SearchBody("end_effector-1"));
	//std::vector<std::shared_ptr<ChBodyAuxRef>> robot_bodylist = { base, shoulder, biceps, elbow, forearm, wrist, end_effector };
	//IRB6700-MH3_150-320_IRC5_rev04_LINK05_CAD.step-1
	//beam_011025-1

	auto base = std::dynamic_pointer_cast<ChBodyAuxRef>(sys.SearchBody("IRB6700-MH3_150-320_IRC5_rev01_BASE_CAD-1"));
	auto shoulder = std::dynamic_pointer_cast<ChBodyAuxRef>(sys.SearchBody("IRB6700_155-285_IRC5_rev03_LINK01_CAD-1"));
	auto biceps = std::dynamic_pointer_cast<ChBodyAuxRef>(sys.SearchBody("IRB6700_155-285_IRC5_rev01_LINK02_CAD-1"));
	auto elbow = std::dynamic_pointer_cast<ChBodyAuxRef>(sys.SearchBody("IRB6700_155-285_IRC5_rev02_LINK03_CAD-1"));
	auto forearm = std::dynamic_pointer_cast<ChBodyAuxRef>(sys.SearchBody("IRB6700_155-285_IRC5_rev01_LINK04_CAD-1"));
	auto wrist = std::dynamic_pointer_cast<ChBodyAuxRef>(sys.SearchBody("IRB6700-MH3_150-320_IRC5_rev04_LINK05_CAD.step-1"));
	auto end_effector = std::dynamic_pointer_cast<ChBodyAuxRef>(sys.SearchBody("IRB6700_155-285_IRC5_rev01_LINK06_CAD-1"));
	auto beam = std::dynamic_pointer_cast<ChBodyAuxRef>(sys.SearchBody("beam_011025-1"));
	auto augerout = std::dynamic_pointer_cast<ChBodyAuxRef>(sys.SearchBody("Auger_002_RB_out-1"));
	auto augerin = std::dynamic_pointer_cast<ChBodyAuxRef>(sys.SearchBody("Auger_In-1"));
	std::vector<std::shared_ptr<ChBodyAuxRef>> robot_bodylist = { base, shoulder, biceps, elbow, forearm, wrist, end_effector, beam, augerout, augerin };

	//auto marker1 = sys.SearchMarker("MARKER_1");
	//auto marker2 = sys.SearchMarker("MARKER_2");
	//auto marker3 = sys.SearchMarker("MARKER_3");
	//auto marker4 = sys.SearchMarker("MARKER_4");
	//auto marker5 = sys.SearchMarker("MARKER_5");
	//auto marker6 = sys.SearchMarker("MARKER_6");
	//auto markerTCP = sys.SearchMarker("MARKER_TCP");
	//std::vector<std::shared_ptr<ChMarker>> robot_markerlist = { marker1, marker2, marker3, marker4, marker5, marker6, markerTCP };

	auto marker1 = sys.SearchMarker("Link01");
	auto marker2 = sys.SearchMarker("Link02");
	auto marker3 = sys.SearchMarker("Link03");
	auto marker4 = sys.SearchMarker("Link04");
	auto marker5 = sys.SearchMarker("Link05");
	auto marker6 = sys.SearchMarker("Link06");
	auto marker7 = sys.SearchMarker("Link07");
	auto marker8 = sys.SearchMarker("Link08");
	auto marker9 = sys.SearchMarker("Link09");
	auto markerTCP = sys.SearchMarker("Link10");
	std::vector<std::shared_ptr<ChMarker>> robot_markerlist = { marker1, marker2, marker3, marker4, marker5, marker6, marker7, marker8, marker9 };

	std::vector<std::shared_ptr<ChLinkMateRevolute>> robot_linklist = {};
	for (auto i = 0; i < robot_bodylist.size() - 1; ++i) {
		robot_linklist.push_back(CreateRevolute(robot_bodylist[i + 1], robot_bodylist[i], robot_markerlist[i]->GetAbsFrame()));
	    sys.Add(robot_linklist[i]);
	}
	end_effector->GetVisualShape(0)->SetOpacity(0.5f);

	// Ground
	auto ground = chrono_types::make_shared<ChBody>();
	ground->SetFixed(true);
	sys.Add(ground);

	//
	auto augerout_ref = augerout->GetFrameRefToAbs();
	auto augerin_ref  = augerin->GetFrameRefToAbs();
	
	auto p_out = augerout_ref.GetPos();
	auto q_out = augerout_ref.GetRot();
	
	auto p_in  = augerin_ref.GetPos();
	auto q_in  = augerin_ref.GetRot();


	// Load and scale mesh
	auto trimesh = chrono_types::make_shared<ChTriangleMeshConnected>();
	double scale_ratio = 1.0;
	trimesh->LoadWavefrontMesh("body_10_1.obj", false, true);
	//trimesh->Transform(ChVector3d(0, 0, 0), ChMatrix33<>(scale_ratio));
	trimesh->Transform(p_out, ChMatrix33<>(scale_ratio));
	//trimesh->RepairDuplicateVertexes(1e-9);  // optional but recommended

	std::cout << "Mesh triangles: " << trimesh->GetNumTriangles() << "\n";
	std::cout << "Mesh vertices : " << trimesh->GetNumVertices()  << "\n";

	// Create contact material
	auto mat = chrono_types::make_shared<ChContactMaterialNSC>();
	mat->SetFriction(0.6f);

	// Enable collision and assign mesh
	auto coll_shape = chrono_types::make_shared<ChCollisionShapeTriangleMesh>(mat, trimesh, false, false, 0.001);
	augerout->AddCollisionShape(coll_shape, ChFrame<>()); 
	augerout->EnableCollision(true);

	std::cout << "augerout collision enabled: " << augerout->IsCollisionEnabled() << "\n";
	std::cout << "augerout num collision shapes: " << augerout->GetCollisionModel()->GetNumShapes() << "\n";

	// //test to place a box 
	// auto box = chrono_types::make_shared<ChBodyEasyBox>(
	// 	0.1, 0.1, 0.1,    // box size (X, Y, Z)
	// 	1000,             // density
	// 	true,             // visualization (true to create a visual shape automatically)
	// 	true              // collision enabled
	// );

	// box->SetPos(augerout->GetPos()+ ChVector3d(0, 0.3, -0.1));
	// box->SetRot(augerout->GetRot());
	// box->SetName("CollisionBox");
	// box->SetFixed(true); // If you want it fixed in space, set to true
	// box->EnableCollision(true);

	// // Add to system
	// sys.Add(box);

	// Create a basic particle emitter
	ChParticleEmitter emitter;

	emitter.ParticlesPerSecond() = 50;  // 10 particles/sec
	emitter.SetUseParticleReservoir(true);
	emitter.ParticleReservoirAmount() = 500;

	// Set outlet near the beam
	auto emitter_positions = chrono_types::make_shared<ChRandomParticlePositionRectangleOutlet>();
	emitter_positions->Outlet() = ChCoordsys<>(augerout->GetPos() + ChVector3d(0, 0.3, -0.1), augerout->GetRot());
	emitter_positions->OutletWidth() = 0.1;
	emitter_positions->OutletHeight() = 0.1;
	emitter.SetParticlePositioner(emitter_positions);
	std::cout << "augerout pos: " << augerout->GetPos() << "\n";
	std::cout << "emitter outlet pos: " << emitter_positions->Outlet().pos << "\n";

	// Use constant downward velocity
	auto mvelo = chrono_types::make_shared<ChRandomParticleVelocityConstantDirection>();
	mvelo->SetDirection(-VECT_Y);  // Downward
	mvelo->SetModulusDistribution(0.5); 
	emitter.SetParticleVelocity(mvelo);


	// // Create box-shaped particles
	// auto creator = chrono_types::make_shared<ChRandomShapeCreatorBoxes>();
	// creator->SetXsizeDistribution(chrono_types::make_shared<ChConstantDistribution>(0.2));
	// creator->SetSizeRatioZDistribution(chrono_types::make_shared<ChConstantDistribution>(1.0));
	// creator->SetSizeRatioYZDistribution(chrono_types::make_shared<ChConstantDistribution>(1.0));
	// creator->SetDensityDistribution(chrono_types::make_shared<ChConstantDistribution>(1000));
	// emitter.SetParticleCreator(creator);

	// Create the contact material
	auto particle_mat = chrono_types::make_shared<ChContactMaterialNSC>();
	particle_mat->SetFriction(0.5f);
	particle_mat->SetRestitution(0.1f);

	// Create the sphere creator
	auto creator = chrono_types::make_shared<ChRandomShapeCreatorSpheres>();
	creator->SetDiameterDistribution(chrono_types::make_shared<ChConstantDistribution>(0.01));  // Sphere diameter
	creator->SetDensityDistribution(chrono_types::make_shared<ChConstantDistribution>(1000));

	// Create a box to show as the emitter 
	auto outlet_marker = chrono_types::make_shared<chrono::ChBodyEasyBox>(0.1, 0.05, 0.1, 1000, true, true);
	outlet_marker->SetPos(augerout->GetPos() + ChVector3d(0, 0.4, -0.1));
	outlet_marker->SetRot(augerout->GetRot());
	outlet_marker->SetFixed(true);
	outlet_marker->EnableCollision(false);

	sys.Add(outlet_marker);

	// G-code parser
	ChParserGCode gparser(gcode_file.string());
	auto gparser_path = gparser.GetPath();
	std::cout << "Continuity error: " << gparser_path->GetContinuityMaxError() << std::endl;

	ground->AddVisualModel(gparser.GetVisualModel(markerTCP->GetAbsFrame()));

	// Imposed trajectory
	auto funpos = chrono_types::make_shared<ChFunctionPositionLine>();
	funpos->SetSpaceFunction(chrono_types::make_shared<ChFunctionConst>(0)); // Hold at start
	funpos->SetLine(gparser_path);
	//auto motfun = chrono_types::make_shared<ChFunctionRamp>(0, 3);
	//funpos->SetSpaceFunction(gparser.GetSpaceFunction());

	
	std::vector<ChQuaterniond> rots = { QUNIT, QUNIT };
	auto rotfun = chrono_types::make_shared<ChFunctionRotationBSpline>(1, rots);

	auto motionimposed = chrono_types::make_shared<ChLinkMotionImposed>();
	motionimposed->Initialize(beam, base, markerTCP->GetAbsFrame());
	motionimposed->SetPositionFunction(funpos);
	motionimposed->SetRotationFunction(rotfun);
	sys.Add(motionimposed);

	// Fix the beam to the end_effector
	auto fix_beam_to_effector = chrono_types::make_shared<ChLinkMateGeneric>();
	fix_beam_to_effector->SetConstrainedCoords(true, true, true, true, true, true);  // Fix all 6 DOF
	fix_beam_to_effector->Initialize(beam, end_effector, beam->GetFrameRefToAbs());
	sys.Add(fix_beam_to_effector);

	// Fix the augerout to the beam
	auto fix_augerout_to_beam = chrono_types::make_shared<ChLinkMateGeneric>();
	fix_augerout_to_beam->SetConstrainedCoords(true, true, true, true, true, true);  // Fix all 6 DOF
	fix_augerout_to_beam->Initialize(augerout, beam, augerout->GetFrameRefToAbs());
	sys.Add(fix_augerout_to_beam);

	////////REFEENCE POINT///////////////////////////////////////
    auto rodRefBody = chrono_types::make_shared<ChBody>();
    rodRefBody->SetPos(augerout->GetPos());
    //rodRefBody->SetFixed(true);
    sys.AddBody(rodRefBody);

	auto motor = chrono_types::make_shared<ChLinkMotorRotationSpeed>();
    motor->SetName("Auger");
    motor->Initialize(augerin, rodRefBody, ChFrame<>(augerin->GetPos(),chrono::QuatFromAngleAxis(CH_PI/2.0, VECT_X)));
	auto msp = chrono_types::make_shared<ChFunctionConst>(0);
	motor->SetSpeedFunction(msp); 
    //motor->SetSpeedFunction(chrono_types::make_shared<ChFunctionConst>(100));
    sys.Add(motor);

	//add a floor to check 
	auto floor = chrono_types::make_shared<ChBodyEasyBox>(20.0, 0.1, 20.0, 1000, true, true, particle_mat);
	floor->SetPos(ChVector3d(0, -0.2, 0));
	floor->SetFixed(true);
	auto floor_vmat = chrono_types::make_shared<ChVisualMaterial>();
	floor_vmat->SetDiffuseColor(ChColor(0.16f, 0.16f, 0.16f));  // concrete-like dark grey floor
	if (auto floor_vshape = floor->GetVisualShape(0)) {
		floor_vshape->SetMaterial(0, floor_vmat);
	}
	sys.Add(floor);
	std::cout << "floor pos: " << floor->GetPos() << " size=(20,0.1,20)\n";


	// Irrlicht
	ChVisualSystemIrrlicht vis;
	vis.AttachSystem(&sys);
	vis.SetWindowSize(800, 600);
	vis.SetWindowTitle("Chrono::Engine");	
	vis.Initialize();
	vis.AddLogo();
	vis.AddSkyBox();
	vis.AddTypicalLights();
	vis.AddCamera(ChVector3d(2, 2, 0), ChVector3d(0, 0, 0));

	vis.EnableCollisionShapeDrawing(true);
	vis.EnableBodyFrameDrawing(true);
	
	vis.BindItem(floor);

	// Visual callback
	class MyEmitterCallback : public ChRandomShapeCreator::AddBodyCallback {
	public:
		MyEmitterCallback(ChVisualSystem* vis, ChCollisionSystem* collsys) : vis(vis), collsys(collsys) {}

		void OnAddBody(std::shared_ptr<ChBody> body,
					ChCoordsys<>,
					ChRandomShapeCreator&) override {
			vis->BindItem(body);  // Bind visuals

			auto vmat = chrono_types::make_shared<ChVisualMaterial>();
			vmat->SetDiffuseColor(ChColor(0.3f, 0.6f, 1.0f));

			// Chrono 9: no GetNumVisualShapes(); just check pointer
			if (auto vshape = body->GetVisualShape(0)) {
				vshape->SetMaterial(0, vmat);
			}

			if (auto cm = body->GetCollisionModel()) {
				body->EnableCollision(true);
				collsys->Add(cm);
			}
		}

		ChVisualSystem* vis;
		ChCollisionSystem* collsys;
	};

	auto callback = chrono_types::make_shared<MyEmitterCallback>(&vis, sys.GetCollisionSystem().get());
	creator->RegisterAddBodyCallback(callback);



	// 	ChVisualSystem* vis;
	// 	ChCollisionSystem* coll;
	// 	std::shared_ptr<ChContactMaterial> mat;
	// };

	// // Create and configure callback
	// auto mcreation_callback = chrono_types::make_shared<MyCreatorForAll>();
	// mcreation_callback->vis = &vis;
	// mcreation_callback->coll = sys.GetCollisionSystem().get();
	// mcreation_callback->mat = particle_mat;
	// emitter.RegisterAddBodyCallback(mcreation_callback);

	// Connect creator to the emitter
	emitter.SetParticleCreator(creator);

	// Simulation loop
	ChRealtimeStepTimer realtime_timer;
	double timestep = 0.001;

	sys.SetSolverType(ChSolver::Type::MINRES);
	sys.GetSolver()->AsIterative()->SetTolerance(1e-5);
	sys.GetSolver()->AsIterative()->SetMaxIterations(200);

	//sys.SetSolverType(ChSolver::Type::BARZILAIBORWEIN);
    //sys.GetSolver()->AsIterative()->SetMaxIterations(400);

	//sys.DoAssembly(ChAssemblyAnalysis::AssemblyLevel::POSITION);
	//sys.DoAssembly(chrono::ChAssemblyAnalysis::POSITION);  // If POSITION is a direct argument
	//sys.DoAssembly(AssemblyLevel::FULL);

	double simulation_time = 0.0;
	double pause_duration = 10.0;  // pause for 10 seconds
	int step_count = 0;
	//double timestep = 0.001;       // your timestep value

	while (vis.Run()) {
		vis.BeginScene();
		vis.Render();
		tools::drawGrid(&vis, 0.5, 0.5, 10, 10, ChCoordsysd(VNULL, Q_ROTATE_Z_TO_Y), ChColor(0.7f, 0.7f, 0.7f), true);
		vis.EndScene();

		if (simulation_time >= pause_duration) {
		// After 10 seconds: Enable motion
		funpos->SetSpaceFunction(gparser.GetSpaceFunction());
		auto msp = chrono_types::make_shared<ChFunctionConst>(10);
	} else {
		// Freeze motion before 10s
		funpos->SetSpaceFunction(chrono_types::make_shared<ChFunctionConst>(0));
		auto msp = chrono_types::make_shared<ChFunctionConst>(0);
		emitter.EmitParticles(sys, timestep);
		std::cout << "num bodies in system: " << sys.GetNumBodies() << "\n";

	}
	sys.DoStepDynamics(timestep);
		if ((step_count % 100) == 0) {
			std::cout << "[step " << step_count << "] bodies=" << sys.GetNumBodies()
			          << " contacts=" << sys.GetContactContainer()->GetNumContacts()
			          << "\n";
		}
	realtime_timer.Spin(timestep);
	simulation_time += timestep;
	++step_count;

	}

	return 0;



	// 	while (vis.Run()) {
	// 		vis.BeginScene();
	// 		vis.Render();
	// 		tools::drawGrid(&vis, 0.5, 0.5, 10, 10, ChCoordsysd(VNULL, Q_ROTATE_Z_TO_Y), ChColor(0.7f, 0.7f, 0.7f), true);
	// 		vis.EndScene();

	// 		if (!vis.GetUtilityFlag()) {
	// 			emitter.EmitParticles(sys, timestep);
	// 			// Use the processor to count particle flow in the rectangle section:
	// 			//processor_flowcount.ProcessParticles(sys);
	// 			//std::cout << "Particles being flown across rectangle:" << counter->counter << "\n";
				
	// 			sys.DoStepDynamics(timestep);
	// 			realtime_timer.Spin(timestep);
	// 		}
	// 	}

}
