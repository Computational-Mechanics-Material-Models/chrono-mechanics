
// ================================================================================
// CHRONO WORKBENCH - github.com/Concrete-Chrono-Development/chrono-preprocessor
//
// Copyright (c) 2023 
// All rights reserved. 
//
// Use of the code that generated this file is governed by a BSD-style license that
// can be found in the LICENSE file at the top level of the distribution and at
// github.com/Concrete-Chrono-Development/chrono-preprocessor/blob/main/LICENSE
//
// ================================================================================
// Chrono Input File
// ================================================================================
//
//
// ================================================================================
#include <math.h>
#include <filesystem>


#include "chrono/collision/ChCollisionModel.h"

#include "chrono/core/ChQuaternion.h"
#include "chrono/core/ChCoordsys.h"
#include "chrono/core/ChMatrix.h"
#include "chrono/core/ChFrame.h"

#include "chrono_thirdparty/rapidjson/prettywriter.h"
#include "chrono_thirdparty/rapidjson/stringbuffer.h"
//#include "chrono_thirdparty/filesystem/path.h"

// #include "chrono_irrlicht/ChVisualSystemIrrlicht.h"

#include <chrono/physics/ChBody.h>
#include <chrono/physics/ChSystemSMC.h>
#include <chrono/physics/ChLinkMate.h>
#include <chrono/physics/ChBodyEasy.h>
// #include <chrono/solver/ChIterativeSolverLS.h>
// #include <chrono/solver/ChDirectSolverLS.h>
#include <chrono/timestepper/ChTimestepper.h>
#include <chrono/timestepper/ChTimestepperHHT.h>
#include "chrono_pardisomkl/ChSolverPardisoMKL.h"
//#include "chrono/solver/ChSolverADMM.h"
#include "chrono/physics/ChLinkMotionImposed.h"



#include "chrono/utils/ChUtilsSamplers.h"
#include "chrono/utils/ChUtilsCreators.h"
#include "chrono/utils/ChUtilsGenerators.h"
#include "chrono/utils/ChUtilsGeometry.h"

#include "chrono/fea/ChElementSpring.h"
#include "chrono/fea/ChLinkNodeFrame.h"
#include "chrono/fea/ChLinkNodeSlopeFrame.h"
#include "chrono/fea/ChMesh.h"

#include "chrono/geometry/ChCylinder.h"
#include "chrono/geometry/ChBox.h"

#include "chrono/timestepper/ChTimestepper.h"

#include "chrono/solver/ChIterativeSolverLS.h"

#include "chrono/assets/ChVisualShapePointPoint.h"
#include "chrono/assets/ChTexture.h"

#include "chrono/physics/ChLinkMotorLinearPosition.h"
#include "chrono/physics/ChLinkMotorRotationSpeed.h"



#include <chrono/functions/ChFunctionPositionXYZFunctions.h>
#include "chrono/physics/ChLinkMotionImposed.h"
      

#include <chrono/physics/ChBody.h>
#include <chrono/physics/ChSystemSMC.h>
#include <chrono/physics/ChLinkMate.h>
//#include <chrono/physics/ChBodyEasy.h>
#include <chrono/solver/ChIterativeSolverLS.h>
#include <chrono/solver/ChDirectSolverLS.h>
#include <chrono/timestepper/ChTimestepper.h>
#include <chrono/timestepper/ChTimestepperHHT.h>
#include "chrono_pardisomkl/ChSolverPardisoMKL.h"
#include "chrono/solver/ChSolverADMM.h"



#include "chrono_ldpm/ChElementLDPM.h"
#include "chrono_ldpm/ChBuilderLDPM.h"
#include "chrono_ldpm/ChMaterialVECT.h"
#include "chrono_ldpm/ChSectionLDPM.h"

//#include "ChElementLDPM.h"
//#include "ChBuilderLDPM.h"
//#include "ChMaterialVECT.h"
//#include "ChSectionLDPM.h"

#include "chrono/fea/ChElementTetraCorot_4.h"
#include "chrono/fea/ChMeshFileLoader.h"

#include "chrono/fea/ChNodeFEAxyzrot.h"
#include "chrono/fea/ChElementBeamEuler.h"
#include "chrono/fea/ChBuilderBeam.h"
#include "chrono/fea/ChMesh.h"
#include "chrono/fea/ChLinkNodeFrame.h"
#include "chrono/fea/ChLinkNodeSlopeFrame.h"
#include "chrono/assets/ChVisualShapeFEA.h"
//#include <chrono_irrlicht/ChVisualSystemIrrlicht.h>



using namespace chrono;
//using namespace chrono::geometry;
using namespace chrono::fea;
//using namespace chrono::irrlicht;
using namespace chrono::ldpm;
//using namespace irr;

/*
void WriteMesh(std::shared_ptr<ChMesh> mesh, const std::string& mesh_filename) {
    std::ofstream out_stream;
    out_stream.open(mesh_filename, std::ios::out);
    out_stream.precision(7);
    out_stream << std::scientific;

    std::vector<std::vector<int>> CableElemNodes;
    std::vector<std::vector<int>> ShellElemNodes;
    std::vector<std::vector<int>> BrickElemNodes;
    std::vector<std::vector<int>> BeamElemNodes;
    std::vector<std::vector<int>> TetElemNodes;

    std::vector<std::shared_ptr<ChNodeFEAbase>> myvector;
    myvector.resize(mesh->GetNumNodes());

    for (unsigned int i = 0; i < mesh->GetNumNodes(); i++) {
        myvector[i] = std::dynamic_pointer_cast<ChNodeFEAbase>(mesh->GetNode(i));
    }

    int numCables = 0;
    int numShells = 0;
    int numBricks = 0;
    int numBeams = 0;
    int numTets = 0;

    for (unsigned int iele = 0; iele < mesh->GetNumElements(); iele++) {
        if (std::dynamic_pointer_cast<ChElementCableANCF>(mesh->GetElement(iele)))
            numCables++;
        if (std::dynamic_pointer_cast<ChElementShellANCF_3423>(mesh->GetElement(iele)))
            numShells++;
        if (std::dynamic_pointer_cast<ChElementHexaANCF_3813>(mesh->GetElement(iele)))
            numBricks++;
        if (std::dynamic_pointer_cast<ChElementBeamEuler>(mesh->GetElement(iele)))
            numBeams++;
        //if (std::dynamic_pointer_cast<ChElementCSL>(mesh->GetElement(iele)))
            //    numBeams++;
        if (std::dynamic_pointer_cast<ChElementLDPM>(mesh->GetElement(iele)))
            numTets++;
        if (std::dynamic_pointer_cast<ChElementTetraCorot_4>(mesh->GetElement(iele)))
            numTets++;
    }
    out_stream << "\nCELLS " << mesh->GetNumElements() << " "
        << (unsigned int)(numCables * 3 + numShells * 5 + numBricks * 9 + numBeams * 3 + numTets * 5) << "\n";

    for (unsigned int iele = 0; iele < mesh->GetNumElements(); iele++) {
        std::vector<int> mynodes;

        if (auto elementC = std::dynamic_pointer_cast<ChElementCableANCF>(mesh->GetElement(iele))) {
            mynodes.resize(2);
            out_stream << "2 ";
            int nodeOrder[] = { 0, 1 };
            mynodes[0] = elementC->GetNode(nodeOrder[0])->GetIndex();
            mynodes[1] = elementC->GetNode(nodeOrder[1])->GetIndex();
            CableElemNodes.push_back(mynodes);
            for (int myNodeN = 0; myNodeN < mynodes.size(); myNodeN++) {
                auto nodeA = (elementC->GetNode(nodeOrder[myNodeN]));
                std::vector<std::shared_ptr<ChNodeFEAbase>>::iterator it;
                it = find(myvector.begin(), myvector.end(), nodeA);
                if (it == myvector.end()) {
                    // name not in vector
                }
                else {
                    auto index = std::distance(myvector.begin(), it);
                    out_stream << (unsigned int)index << " ";
                }
            }
            out_stream << "\n";
        }
        else if (auto elementBm = std::dynamic_pointer_cast<ChElementBeamEuler>(mesh->GetElement(iele))) {
            mynodes.resize(2);
            out_stream << "2 ";
            int nodeOrder[] = { 0, 1 };
            mynodes[0] = elementBm->GetNode(nodeOrder[0])->GetIndex();
            mynodes[1] = elementBm->GetNode(nodeOrder[1])->GetIndex();
            BeamElemNodes.push_back(mynodes);
            for (int myNodeN = 0; myNodeN < mynodes.size(); myNodeN++) {
                auto nodeA = (elementBm->GetNode(nodeOrder[myNodeN]));
                std::vector<std::shared_ptr<ChNodeFEAbase>>::iterator it;
                it = find(myvector.begin(), myvector.end(), nodeA);
                if (it == myvector.end()) {
                    // name not in vector
                }
                else {
                    auto index = std::distance(myvector.begin(), it);
                    out_stream << (unsigned int)index << " ";
                }
            }
            out_stream << "\n";
        }
        else if (auto elementS = std::dynamic_pointer_cast<ChElementShellANCF_3423>(mesh->GetElement(iele))) {
            mynodes.resize(4);
            out_stream << "4 ";
            int nodeOrder[] = { 0, 1, 2, 3 };
            mynodes[0] = elementS->GetNode(nodeOrder[0])->GetIndex();
            mynodes[1] = elementS->GetNode(nodeOrder[1])->GetIndex();
            mynodes[2] = elementS->GetNode(nodeOrder[2])->GetIndex();
            mynodes[3] = elementS->GetNode(nodeOrder[3])->GetIndex();
            ShellElemNodes.push_back(mynodes);
            for (int myNodeN = 0; myNodeN < mynodes.size(); myNodeN++) {
                auto nodeA = (elementS->GetNode(nodeOrder[myNodeN]));
                std::vector<std::shared_ptr<ChNodeFEAbase>>::iterator it;
                it = find(myvector.begin(), myvector.end(), nodeA);
                if (it == myvector.end()) {
                    // name not in vector
                }
                else {
                    auto index = std::distance(myvector.begin(), it);
                    out_stream << (unsigned int)index << " ";
                }
            }
            out_stream << "\n";
        }
        else if (auto elementB = std::dynamic_pointer_cast<ChElementHexaANCF_3813>(mesh->GetElement(iele))) {
            mynodes.resize(8);
            out_stream << "8 ";
            int nodeOrder[] = { 0, 1, 2, 3, 4, 5, 6, 7 };
            mynodes[0] = elementB->GetNode(nodeOrder[0])->GetIndex();
            mynodes[1] = elementB->GetNode(nodeOrder[1])->GetIndex();
            mynodes[2] = elementB->GetNode(nodeOrder[2])->GetIndex();
            mynodes[3] = elementB->GetNode(nodeOrder[3])->GetIndex();
            mynodes[4] = elementB->GetNode(nodeOrder[4])->GetIndex();
            mynodes[5] = elementB->GetNode(nodeOrder[5])->GetIndex();
            mynodes[6] = elementB->GetNode(nodeOrder[6])->GetIndex();
            mynodes[7] = elementB->GetNode(nodeOrder[7])->GetIndex();
            BrickElemNodes.push_back(mynodes);
            for (int myNodeN = 0; myNodeN < mynodes.size(); myNodeN++) {
                auto nodeA = (elementB->GetNode(nodeOrder[myNodeN]));
                std::vector<std::shared_ptr<ChNodeFEAbase>>::iterator it;
                it = find(myvector.begin(), myvector.end(), nodeA);
                if (it == myvector.end()) {
                    // name not in vector
                }
                else {
                    auto index = std::distance(myvector.begin(), it);
                    out_stream << (unsigned int)index << " ";
                }
            }
            out_stream << "\n";
        }
        else if (auto elementB = std::dynamic_pointer_cast<ChElementLDPM>(mesh->GetElement(iele))) {
            mynodes.resize(4);
            out_stream << "4 ";
            int nodeOrder[] = { 0, 1, 2, 3 };
            mynodes[0] = elementB->GetNode(nodeOrder[0])->GetIndex();
            mynodes[1] = elementB->GetNode(nodeOrder[1])->GetIndex();
            mynodes[2] = elementB->GetNode(nodeOrder[2])->GetIndex();
            mynodes[3] = elementB->GetNode(nodeOrder[3])->GetIndex();
            TetElemNodes.push_back(mynodes);
            for (int myNodeN = 0; myNodeN < mynodes.size(); myNodeN++) {
                auto nodeA = (elementB->GetNode(nodeOrder[myNodeN]));
                std::vector<std::shared_ptr<ChNodeFEAbase>>::iterator it;
                it = find(myvector.begin(), myvector.end(), nodeA);
                if (it == myvector.end()) {
                    // name not in vector
                }
                else {
                    auto index = std::distance(myvector.begin(), it);
                    out_stream << (unsigned int)index << " ";
                }
            }
            out_stream << "\n";
        }
        else if (auto elementB = std::dynamic_pointer_cast<ChElementTetraCorot_4>(mesh->GetElement(iele))) {
            mynodes.resize(4);
            out_stream << "4 ";
            int nodeOrder[] = { 0, 1, 2, 3 };
            mynodes[0] = elementB->GetNode(nodeOrder[0])->GetIndex();
            mynodes[1] = elementB->GetNode(nodeOrder[1])->GetIndex();
            mynodes[2] = elementB->GetNode(nodeOrder[2])->GetIndex();
            mynodes[3] = elementB->GetNode(nodeOrder[3])->GetIndex();
            TetElemNodes.push_back(mynodes);
            for (int myNodeN = 0; myNodeN < mynodes.size(); myNodeN++) {
                auto nodeA = (elementB->GetNode(nodeOrder[myNodeN]));
                std::vector<std::shared_ptr<ChNodeFEAbase>>::iterator it;
                it = find(myvector.begin(), myvector.end(), nodeA);
                if (it == myvector.end()) {
                    // name not in vector
                }
                else {
                    auto index = std::distance(myvector.begin(), it);
                    out_stream << (unsigned int)index << " ";
                }
            }
            out_stream << "\n";
        }
    }

    out_stream << "\nCELL_TYPES " << mesh->GetNumElements() << "\n";

    for (unsigned int iele = 0; iele < mesh->GetNumElements(); iele++) {
        if (std::dynamic_pointer_cast<ChElementCableANCF>(mesh->GetElement(iele)))
            out_stream << "3\n";
        else if (std::dynamic_pointer_cast<ChElementBeamEuler>(mesh->GetElement(iele)))
            out_stream << "3\n";
        //else if (std::dynamic_pointer_cast<ChElementCSL>(mesh->GetElement(iele)))
             //   out_stream << "3\n";
        else if (std::dynamic_pointer_cast<ChElementShellANCF_3423>(mesh->GetElement(iele)))
            out_stream << "9\n";
        else if (std::dynamic_pointer_cast<ChElementHexaANCF_3813>(mesh->GetElement(iele)))
            out_stream << "12\n";
        else if (std::dynamic_pointer_cast<ChElementLDPM>(mesh->GetElement(iele)))
            out_stream << "10\n";
        else if (std::dynamic_pointer_cast<ChElementTetraCorot_4>(mesh->GetElement(iele)))
            out_stream << "10\n";
    }

    out_stream.close();
}

void WriteFrame(std::shared_ptr<ChMesh> mesh,
    const std::string& mesh_filename,
    const std::string& vtk_filename) {
    std::ofstream out_stream;
    out_stream.open(vtk_filename, std::ios::trunc);

    out_stream << "# vtk DataFile Version 2.0" << std::endl;
    out_stream << "Unstructured Grid Example" << std::endl;
    out_stream << "ASCII" << std::endl;
    out_stream << "DATASET UNSTRUCTURED_GRID" << std::endl;


    out_stream << "POINTS " << mesh->GetNumNodes() << " double\n";

    for (unsigned int i = 0; i < mesh->GetNumNodes(); i++) {
        if (auto node = std::dynamic_pointer_cast<ChNodeFEAxyzrot>(mesh->GetNode(i))) {
            out_stream << node->GetPos().x() << " " << node->GetPos().y() << " " << node->GetPos().z() << "\n";
        }

        if (auto node = std::dynamic_pointer_cast<ChNodeFEAxyz>(mesh->GetNode(i))) {
            out_stream << node->GetPos().x() << " " << node->GetPos().y() << " " << node->GetPos().z() << "\n";
        }
    }


    std::ifstream in_stream(mesh_filename);
    out_stream << in_stream.rdbuf();

    int numCell = 0;
    for (unsigned int iele = 0; iele < mesh->GetNumElements(); iele++) {
        if (std::dynamic_pointer_cast<ChElementCableANCF>(mesh->GetElement(iele)))
            numCell++;
        else if (std::dynamic_pointer_cast<ChElementBeamEuler>(mesh->GetElement(iele)))
            numCell++;
        //else if (std::dynamic_pointer_cast<ChElementCSL>(mesh->GetElement(iele)))
            //    numCell++;
        else if (std::dynamic_pointer_cast<ChElementShellANCF_3423>(mesh->GetElement(iele)))
            numCell++;
        else if (std::dynamic_pointer_cast<ChElementLDPM>(mesh->GetElement(iele)))
            numCell++;
        else if (std::dynamic_pointer_cast<ChElementTetraCorot_4>(mesh->GetElement(iele)))
            numCell++;
    }

    out_stream << "\nCELL_DATA " << numCell << "\n";
    out_stream << "SCALARS Force double\n";
    out_stream << "LOOKUP_TABLE default\n";

    double scalar = 0;
    ChVector3d Fforce;
    ChVector3d Mtorque;
    for (unsigned int iele = 0; iele < mesh->GetNumElements(); iele++) {
        if (auto elementC = std::dynamic_pointer_cast<ChElementCableANCF>(mesh->GetElement(iele)))
            scalar = elementC->GetCurrLength() - elementC->GetRestLength();
        //  else if (auto elementC = std::dynamic_pointer_cast<ChElementBeamEuler>(mesh->GetElement(iele)))
        //       elementC->EvaluateSectionForceTorque(0.0, Fforce, Mtorque);          
        else if (auto elementS = std::dynamic_pointer_cast<ChElementShellANCF_3423>(mesh->GetElement(iele)))
            elementS->EvaluateDeflection(scalar);

        out_stream << Fforce[0] + 1e-40 << "\n";
    }

    //out_stream << "\nCELL_DATA " << numCell << "\n";
    out_stream << "VECTORS Strain double\n";
    //out_stream << "LOOKUP_TABLE default\n";
    ChVector3d StrainV;
    for (unsigned int iele = 0; iele < mesh->GetNumElements(); iele++) {
        if (auto elementC = std::dynamic_pointer_cast<ChElementCableANCF>(mesh->GetElement(iele)))
            elementC->EvaluateSectionStrain(0.0, StrainV);
        else if (auto elementC = std::dynamic_pointer_cast<ChElementBeamEuler>(mesh->GetElement(iele)))
            elementC->EvaluateSectionStrain(0.0, StrainV);
        //	else if (auto elementC = std::dynamic_pointer_cast<ChElementCSL>(mesh->GetElement(iele)))
                //elementC->EvaluateSectionStrain(0.0, StrainV);            
        else if (auto elementS = std::dynamic_pointer_cast<ChElementShellANCF_3423>(mesh->GetElement(iele))) {
            const ChStrainStress3D strainStressOut =
                elementS->EvaluateSectionStrainStress(ChVector3d(0, 0, 0), 0);
            StrainV.Set(strainStressOut.strain[0], strainStressOut.strain[1], strainStressOut.strain[3]);
        }
        StrainV += ChVector3d(1e-40);
        out_stream << StrainV.x() << " " << StrainV.y() << " " << StrainV.z() << "\n";
    }


    out_stream << "\nPOINT_DATA " << mesh->GetNumNodes() << "\n";
    out_stream << "VECTORS Displacement double\n";
    //out_stream << "LOOKUP_TABLE default\n";
    for (unsigned int i = 0; i < mesh->GetNumNodes(); i++) {
        if (auto node = std::dynamic_pointer_cast<ChNodeFEAxyz>(mesh->GetNode(i))) {
            ChVector3d disp = std::dynamic_pointer_cast<ChNodeFEAxyz>(mesh->GetNode(i))->GetPos();
            ChVector3d disp0 = std::dynamic_pointer_cast<ChNodeFEAxyz>(mesh->GetNode(i))->GetX0();
            disp -= disp0; //ChVector3d(1e-40);
            out_stream << (double)disp.x() << " " << (double)disp.y() << " " << (double)disp.z() << "\n";
        }
        else if (auto node = std::dynamic_pointer_cast<ChNodeFEAxyzrot>(mesh->GetNode(i))) {
            ChVector3d disp = std::dynamic_pointer_cast<ChNodeFEAxyzrot>(mesh->GetNode(i))->Frame().GetPos();
            ChVector3d disp0 = std::dynamic_pointer_cast<ChNodeFEAxyzrot>(mesh->GetNode(i))->GetX0().GetPos();
            disp -= disp0; //ChVector3d(1e-40);
            out_stream << (double)disp.x() << " " << (double)disp.y() << " " << (double)disp.z() << "\n";
        }

    }

    out_stream << "\nPOINT_DATA " << mesh->GetNumNodes() << "\n";
    out_stream << "VECTORS Velocity double\n";
    //out_stream << "LOOKUP_TABLE default\n";
    for (unsigned int i = 0; i < mesh->GetNumNodes(); i++) {
        if (auto node = std::dynamic_pointer_cast<ChNodeFEAxyz>(mesh->GetNode(i))) {
            ChVector3d vel = std::dynamic_pointer_cast<ChNodeFEAxyz>(mesh->GetNode(i))->GetPosDt();
            vel += ChVector3d(1e-40);
            out_stream << (double)vel.x() << " " << (double)vel.y() << " " << (double)vel.z() << "\n";
        }
        else if (auto node = std::dynamic_pointer_cast<ChNodeFEAxyzrot>(mesh->GetNode(i))) {
            ChVector3d vel = std::dynamic_pointer_cast<ChNodeFEAxyzrot>(mesh->GetNode(i))->Frame().GetPosDt();
            vel += ChVector3d(1e-40);
            out_stream << (double)vel.x() << " " << (double)vel.y() << " " << (double)vel.z() << "\n";
        }

    }

    out_stream << "\nPOINT_DATA " << mesh->GetNumNodes() << "\n";
    out_stream << "VECTORS Acceleration double\n";
    //out_stream << "LOOKUP_TABLE default\n";

    for (unsigned int i = 0; i < mesh->GetNumNodes(); i++) {
        if (auto node = std::dynamic_pointer_cast<ChNodeFEAxyz>(mesh->GetNode(i))) {
            ChVector3d acc = std::dynamic_pointer_cast<ChNodeFEAxyz>(mesh->GetNode(i))->GetPosDt2();
            acc += ChVector3d(1e-40);
            out_stream << (double)acc.x() << " " << (double)acc.y() << " " << (double)acc.z() << "\n";
        }
        else if (auto node = std::dynamic_pointer_cast<ChNodeFEAxyzrot>(mesh->GetNode(i))) {
            ChVector3d acc = std::dynamic_pointer_cast<ChNodeFEAxyzrot>(mesh->GetNode(i))->Frame().GetPosDt2();
            acc += ChVector3d(1e-40);
            out_stream << (double)acc.x() << " " << (double)acc.y() << " " << (double)acc.z() << "\n";
        }

    }

    out_stream.close();
}

void WriteMesh1(std::shared_ptr<ChMesh> mesh, const std::string& mesh_filename) {
    std::ofstream out_stream;
    out_stream.open(mesh_filename, std::ios::out);
    out_stream.precision(7);
    out_stream << std::scientific;

    std::vector<std::vector<int>> TetElemNodes;

    std::vector<std::shared_ptr<ChNodeFEAbase>> myvector;
    myvector.resize(mesh->GetNumNodes());

    for (unsigned int i = 0; i < mesh->GetNumNodes(); i++) {
        myvector[i] = std::dynamic_pointer_cast<ChNodeFEAbase>(mesh->GetNode(i));
    }


    int numTets = 0;

    for (unsigned int iele = 0; iele < mesh->GetNumElements(); iele++) {

        //if (std::dynamic_pointer_cast<ChElementCSL>(mesh->GetElement(iele)))
            //    numBeams++;
        if (std::dynamic_pointer_cast<ChElementLDPM>(mesh->GetElement(iele)))
            numTets++;
    }
    out_stream << "\nCELLS " << (unsigned int)(numTets * 12) << " "
        << (unsigned int)(numTets * 48) << "\n";

    for (unsigned int iele = 0; iele < numTets * 12; iele++) {
        std::vector<int> mynodes;

        mynodes.resize(3);
        out_stream << "3 ";
        mynodes[0] = (iele + 1) * 3 - 3;
        mynodes[1] = (iele + 1) * 3 - 2;
        mynodes[2] = (iele + 1) * 3 - 1;

        TetElemNodes.push_back(mynodes);

        out_stream << mynodes[0] << " " << mynodes[1] << " " << mynodes[2] << " ";
        out_stream << "\n";
    }

    out_stream << "\nCELL_TYPES " << (unsigned int)(numTets * 12) << "\n";

    for (unsigned int iele = 0; iele < numTets * 12; iele++) {
        out_stream << "5\n";

    }

    out_stream.close();
}

void WriteFrame1(std::shared_ptr<ChMesh> mesh,
    const std::string& mesh_filename,
    const std::string& vtk_filename) {
    std::ofstream out_stream;
    out_stream.open(vtk_filename, std::ios::trunc);

    out_stream << "# vtk DataFile Version 2.0" << std::endl;
    out_stream << "Unstructured Grid Example" << std::endl;
    out_stream << "ASCII" << std::endl;
    out_stream << "DATASET UNSTRUCTURED_GRID" << std::endl;

    int nelements = mesh->GetNumElements();
    out_stream << "POINTS " << nelements * 36 << " double\n";

    for (int i = 0; i < nelements; ++i) {
        auto elem = std::dynamic_pointer_cast<ChElementLDPM>(mesh->GetElement(i));

        for (int j = 0; j < 12; j++) {
            auto vertices = elem->GetVertNodeVec(j);
            auto pC = vertices[0]->GetX0().GetPos();
            auto pA = vertices[1]->GetX0().GetPos();
            auto pB = vertices[2]->GetX0().GetPos();
            out_stream << pC << "\n";
            out_stream << pA << "\n";
            out_stream << pB << "\n";
        }

    }



    std::ifstream in_stream(mesh_filename);
    out_stream << in_stream.rdbuf();

    int numCell = 0;
    for (unsigned int iele = 0; iele < mesh->GetNumElements(); iele++) {

        if (std::dynamic_pointer_cast<ChElementLDPM>(mesh->GetElement(iele)))
            numCell++;

    }




    out_stream << "\nCELL_DATA " << numCell * 12 << "\n";
    out_stream << "SCALARS Crack float\n";
    out_stream << "LOOKUP_TABLE default\n";
    //out_stream << "LOOKUP_TABLE default\n";

    for (unsigned int iele = 0; iele < mesh->GetNumElements(); iele++) {

        auto elem = std::dynamic_pointer_cast<ChElementLDPM>(mesh->GetElement(iele));

        for (auto facet : elem->GetSection()) {
            auto statev = facet->Get_StateVar();

            out_stream << statev(11) << "\n";
        }

    }


    out_stream.close();
}
*/



/*
class MyCallback : public ChStaticNonLinearIncremental::LoadIncrementCallback {
          public:
            // Perform updates on the model. This is called before each load scaling.
            // Here we will update all "external" relevan loads.
            virtual void OnLoadScaling(const double load_scaling,              // ranging from 0 to 1
                                       const int iteration_n,                  // actual number of load step
                                       ChStaticNonLinearIncremental* analysis  // back-pointer to this analysis
            ) {
                // Scale the external loads. In our example, just two forces. 
                // Note: if gravity is used, consider scaling also gravity effect, e.g: 
                //    sys.Set_G_acc(load_scaling * ChVector<>(0,-9.8,0))
                cb_loaded_node_1->SetForce(load_scaling * cb_F_node_1);
                //cb_loaded_node_2->SetForce(load_scaling * cb_F_node_2);
            }
            // helper data for the callback
            ChVector<> cb_F_node_1;            
            std::shared_ptr<ChNodeFEAxyzrot> cb_loaded_node_1;
        };
*/
// Instead of using sys.DoStaticNonLinear(), which is quite basic, we will use ChStaticNonLinearIncremental.
// This requires a custom callback for incrementing the external loads:



/*ChVector<double> calculate_Force(std::shared_ptr<ChMesh> mesh){
    unsigned int icons=0;
    ChVector<double> Force(0,0,0);
    for (unsigned int i = 0; i < mesh->GetNnodes(); i++) {
        auto node = std::dynamic_pointer_cast<ChNodeFEAxyzrot>(mesh->GetNode(i)); 
        auto cy=node->GetPos().y();
        if (cy==49.99){         	
        	auto fn=node->GetForce(); 
        	//std::cout<<"fn "<<fn<<"\n";   
        	Force+=fn;
        }       
    }
    
    return Force;

}*/

/*
ChVector<double> calculate_Force(std::vector< std::shared_ptr<ChLinkMateGeneric> > const_list ){
    unsigned int icons=0;   
    ChVector<double> Force(0.,0.,0.); 
    for (auto constraint:const_list) { 
        	auto fn=constraint->Get_react_force(); 
        	//std::cout<<"fn "<<fn<<"\n";   
        	Force+=fn;
    }
    
    return Force;

}
*/



//collision::ChCollisionSystemType collision_type = collision::ChCollisionSystemType::BULLET;

int main(int argc, char** argv) {
    std::cout<<"Starting..." << std::endl;
    SetChronoDataPath(CHRONO_DATA_DIR);
    //
    // Create a Chrono::Engine physical system
    ChSystemSMC sys;
    //sys.SetNumThreads(1,0,0);
    //
    //    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    ///
    /// Give the path and file name of LDPM data
    ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    //std::string current_dir(argv[0]);
    // int pos = current_dir.find_last_of("/\\");
    // current_dir=current_dir.substr(0, pos-5);  
    // std::cout<<"current_dir "<<current_dir;
    std::string current_dir="";  
    std::cout<<"current_dir "<<current_dir<<"\n";
    //  
    std::string LDPM_data_path=current_dir+"LDPMgeo000Box000/";
    std::string LDPM_GeoName="LDPM_debugRegTet";
    //
    std::string out_dir=current_dir+"out/";
    if (!std::filesystem::create_directory(std::filesystem::path(out_dir))) {
        std::cerr << "Error creating directory " << out_dir << std::endl;
        return 1;
    }
    	
    std::string history_filename="hist.dat"; 
    std::ofstream histfile;
    histfile.open(out_dir+history_filename, std::ios::out);
    //
    //	
    // Create ground:   
    // 
    auto mtruss = chrono_types::make_shared<ChBody>();
    mtruss->SetFixed(true);
    sys.Add(mtruss);
    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    ///
    /// Read LDPM Freecad outputs and insert into mesh object
    ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    //
    // Create a mesh object
    //
    auto my_mesh = chrono_types::make_shared<ChMesh>();
    //
    // Create  vectorial material for LDPM
    //
    auto my_mat = chrono_types::make_shared<ChMaterialVECT>();
    my_mat->Set_density(2.4E-9);
    my_mat->Set_E0(30000);
    my_mat->Set_alpha(0.2);
    my_mat->Set_sigmat(3);
    my_mat->Set_sigmas(7.5);
    my_mat->Set_nt(0.2);
    my_mat->Set_lt(250);
    my_mat->Set_Ed(30000);
    my_mat->Set_sigmac0(120);
    my_mat->Set_Hc0(9900);
    my_mat->Set_Hc1(3000);
    my_mat->Set_beta(0);
    my_mat->Set_kc0(3);
    my_mat->Set_kc1(0.5);
    my_mat->Set_kc2(5);
    my_mat->Set_kc3(0.1);
    my_mat->Set_mu0(0.4);
    my_mat->Set_muinf(0);
    my_mat->Set_sigmaN0(600);
    my_mat->Set_kt(0);
    //
    //read_LDPM_info(my_mesh, nodesFilename, elemFilename, facetFilename, tetsFilename, verticesFilename);
    ChBuilderLDPM builder;
    builder.read_LDPM_info(my_mesh, my_mat, LDPM_data_path, LDPM_GeoName);
    //	
    sys.Add(my_mesh);
	
    //unsigned int elnum=my_mesh->GetNumElements();
    //std::cout << "elnum " << elnum<<std::endl;
    
    
	
    // Apply forces on a group of nodes
	//ChVector <> F_node_1(0, 1000, 0);
	std::vector< std::shared_ptr<ChNodeFEAxyzrot> > loaded_nodes;	
	auto loaded_node0 = std::dynamic_pointer_cast<ChNodeFEAxyzrot>(my_mesh->GetNode(0));
    auto loaded_node1 = std::dynamic_pointer_cast<ChNodeFEAxyzrot>(my_mesh->GetNode(1));
    auto loaded_node2 = std::dynamic_pointer_cast<ChNodeFEAxyzrot>(my_mesh->GetNode(2));
    auto loaded_node3 = std::dynamic_pointer_cast<ChNodeFEAxyzrot>(my_mesh->GetNode(3));
    loaded_nodes.push_back(loaded_node0);
    loaded_nodes.push_back(loaded_node1);
    loaded_nodes.push_back(loaded_node2);
    loaded_nodes.push_back(loaded_node3);

	/*auto constr_loaded_node=chrono_types::make_shared<ChLinkMateGeneric>(false, true, true, true, true, true);
	constr_loaded_node->Initialize(loaded_node, mtruss, false, loaded_node->Frame(), loaded_node->Frame());    	     	
	sys.Add(constr_loaded_node);*/
   	
	
    // Apply constraint on a group of nodes
    std::vector< std::shared_ptr<ChLinkMateGeneric> > const_list;	   
    unsigned int icons=0;

  
    //std::vector<int> fixed_node = {0, 1, 2 };
    //for (unsigned int i = 0; i < 3; i++) {
        //auto node = std::dynamic_pointer_cast<ChNodeFEAxyzrot>(my_mesh->GetNode(fixed_node[i]));
        //auto node = std::dynamic_pointer_cast<ChNodeFEAxyzrot>(my_mesh->GetNode(i)); 
        /*auto constr_tie=chrono_types::make_shared<ChLinkMateGeneric>(true, true, true, true, true, true);
        constr_tie->Initialize(node, mtruss, false, node->Frame(), node->Frame());
        sys.Add(constr_tie);
        const_list.push_back(constr_tie);*/
       // node->SetFixed(true);
    //}

 
    
	
	
	
    // We do not want gravity effect on FEA elements in this demo
    my_mesh->SetAutomaticGravity(false);

    // Remember to add the mesh to the system!
    //sys.Add(my_mesh);
    //sys.Set_G_acc(ChVector<>(0, 0, 0));
	
	/*
    auto mvisualizebeamA = chrono_types::make_shared<ChVisualShapeFEA>(my_mesh);
    //mvisualizebeamA->SetFEMdataType(ChVisualShapeFEA::DataType::ELEM_BEAM_MZ);   
    mvisualizebeamA->SetFEMdataType(ChVisualShapeFEA::DataType::NODE_DISP_Y); 
    mvisualizebeamA->SetColorscaleMinMax(-5., 0.);
    mvisualizebeamA->SetSmoothFaces(true);
    mvisualizebeamA->SetWireframe(false);
    my_mesh->AddVisualShapeFEA(mvisualizebeamA);
	

    auto mvisualizemeshcoll = chrono_types::make_shared<ChVisualShapeFEA>(my_mesh);
    mvisualizemeshcoll->SetFEMdataType(ChVisualShapeFEA::DataType::CONTACTSURFACES);
    mvisualizemeshcoll->SetWireframe(true);
    mvisualizemeshcoll->SetDefaultMeshColor(ChColor(1, 0.5, 0));
    my_mesh->AddVisualShapeFEA(mvisualizemeshcoll);
	
	
	
    auto mvisualizebeamC = chrono_types::make_shared<ChVisualShapeFEA>(my_mesh);
    mvisualizebeamC->SetFEMglyphType(ChVisualShapeFEA::GlyphType::NODE_CSYS);
    mvisualizebeamC->SetFEMdataType(ChVisualShapeFEA::DataType::NONE);
    mvisualizebeamC->SetSymbolsThickness(0.006);
    mvisualizebeamC->SetSymbolsScale(0.01);
    mvisualizebeamC->SetZbufferHide(false);
    my_mesh->AddVisualShapeFEA(mvisualizebeamC);

    // Create the Irrlicht visualization system
	
    auto vis = chrono_types::make_shared<ChVisualSystemIrrlicht>();
    vis->SetWindowSize(800, 600);
    vis->SetWindowTitle("Frame analysis");
    vis->Initialize();
    vis->AddLogo();
    vis->AddSkyBox();
    vis->AddTypicalLights();
    vis->AddCamera(ChVector<>(75.0, 75.0, 75.));
    vis->AttachSystem(&sys);
	
    */
    
    // THE SIMULATION LOOP

    // Create a Chrono solver and set solver settings    
	
	/*
    //auto solver = chrono_types::make_shared<ChSolverSparseQR>();   
    auto solver = chrono_types::make_shared<ChSolverSparseLU>();
    sys.SetSolver(solver);
    solver->UseSparsityPatternLearner(true);
    solver->LockSparsityPattern(true);
    //solver->SetVerbose(true);
    
    */
    
    

    auto solver = chrono_types::make_shared<ChSolverPardisoMKL>();    
    sys.SetSolver(solver); 
    solver->UseSparsityPatternLearner(true);
    solver->LockSparsityPattern(true);   
    solver->SetVerbose(false);
    
    
	
    
    /*
    auto solver = chrono_types::make_shared<ChSolverMINRES>();
    sys.SetSolver(solver);
    solver->SetMaxIterations(200);
    solver->SetTolerance(1e-10);
    solver->EnableDiagonalPreconditioner(true);
    solver->EnableWarmStart(true);  // Enable for better convergence when using Euler implicit linearized

    sys.SetSolverForceTolerance(1e-14);
    */
      
     
     
    // Dont forget to comment this part later
    //sys.SetNumThreads(1,0,0);
    //std::cout<<"sys.GetNumThreadsChrono()\t"<<sys.GetNumThreadsChrono()<<std::endl;
    
    
    //sys.SetTimestepperType(ChTimestepper::Type::EULER_EXPLICIT);
    //sys.SetTimestepperType(ChTimestepper::Type::EULER_IMPLICIT);
    //sys.SetTimestepperType(ChTimestepper::Type::RUNGEKUTTA45);
    

    auto mystepper = chrono_types::make_shared<ChTimestepperHHT>(&sys);
    //if (mystepper==ChTimestepper::Type::HHT){
    mystepper->SetAlpha(-0.05); // alpha=-0.2 default value
    mystepper->SetMaxIters(50);
    mystepper->SetAbsTolerances(1e-04, 1e-04);
    //mystepper->SetMode(ChTimestepperHHT::POSITION);
    //mystepper->SetMode(ChTimestepperHHT::ACCELERATION); // Default
    mystepper->SetMinStepSize(1E-15);
    mystepper->SetMaxItersSuccess(4);
    mystepper->SetRequiredSuccessfulSteps(3);
    mystepper->SetStepIncreaseFactor(1.25);
    mystepper->SetStepDecreaseFactor(0.25);
    //mystepper->SetModifiedNewton(true);
    mystepper->SetJacobianUpdateMethod(ChTimestepperHHT::JacobianUpdate::NEVER);    //only at the beginning of the very first step
    // Default: JacobianUpdate::EVERY_STEP
    mystepper->SetVerbose(true);
    mystepper->SetStepControl(false);
    sys.SetTimestepper(mystepper);
    //}   
    
    
    sys.SetTimestepper(mystepper);
    
     
    
	
    //unsigned int totalEl=my_mesh->GetNumelements();
    //std::cout << "elnum " << elnum<<std::endl;
    
    
  
    if (true) { 
     

	/*
    auto f_xyz = chrono_types::make_shared<ChFunctionPosition_XYZfunctions>(); 
    f_xyz->SetFunctionX(chrono_types::make_shared<ChFunction_Ramp>(0,0.1));      
    

    auto f_sequence = chrono_types::make_shared<ChFunction_Sequence>();
    f_sequence->InsertFunct(f_xyz, 0.2, 0, true);
    
    

    auto my_motion_function = chrono_types::make_shared<ChFunction_Ramp>();
    my_motion_function->Set_ang(0.1);

    auto f_sequence = chrono_types::make_shared<ChFunction_Sequence>();
    f_sequence->InsertFunct(my_motion_function, 0.2, 0, true);

    auto my_motion_function2 = chrono_types::make_shared<ChFunction_Ramp>();
    my_motion_function2->Set_ang(-0.1);
    f_sequence->InsertFunct(my_motion_function2, 0.8, 0, true);

    */

    /*
    auto f_xyz = chrono_types::make_shared<ChFunctionPosition_XYZfunctions>();
    auto f_sequence1 = chrono_types::make_shared<ChFunction_Sequence>();
    auto my_motion_function1 = chrono_types::make_shared<ChFunction_Ramp>(0, -0.1);
    auto my_motion_function2 = chrono_types::make_shared<ChFunction_Ramp>(0, 0.1);
    f_sequence1->InsertFunct(my_motion_function1, 0.2, 0, true);
    f_sequence1->InsertFunct(my_motion_function2, 0.8, 0, true);
    f_xyz->SetFunctionX(f_sequence1);
    
    f_xyz->SetFunctionY(chrono_types::make_shared<ChFunction_Ramp>(0, -7));
    auto my_motion_function3 = chrono_types::make_shared<ChFunction_Ramp>(0, 0);
    auto my_motion_function4 = chrono_types::make_shared<ChFunction_Ramp>(0, -7);
    auto f_sequence2 = chrono_types::make_shared<ChFunction_Sequence>();
    f_sequence2->InsertFunct(my_motion_function3, 0.5, 0, true);
    f_sequence2->InsertFunct(my_motion_function4, 0.5, 0, true);
    f_xyz->SetFunctionY(f_sequence2);
    //
    auto impose_1 = chrono_types::make_shared<ChLinkMotionImposed>();
    sys.Add(impose_1);
    impose_1->Initialize(bodyB, mtruss, ChFrame<>(bodyB->GetPos()));
    impose_1->SetPositionFunction(f_xyz);
    */
    
    auto f_sequence1 = chrono_types::make_shared<ChFunctionSequence>();
    auto my_motion_function11 = chrono_types::make_shared<ChFunctionRamp>(0, 0.1);
    auto my_motion_function12 = chrono_types::make_shared<ChFunctionRamp>(0, -0.1);
    auto my_motion_function13 = chrono_types::make_shared<ChFunctionRamp>(0, 0.1);
    auto my_motion_function14 = chrono_types::make_shared<ChFunctionRamp>(0, -0.1);
    f_sequence1->InsertFunct(my_motion_function11, 0.2, 0, true);
    f_sequence1->InsertFunct(my_motion_function12, 2.2, 0, true);
    f_sequence1->InsertFunct(my_motion_function13, 2.2, 0, true);
    f_sequence1->InsertFunct(my_motion_function14, 0.2, 0, true);

    /*
    auto f_sequence2 = chrono_types::make_shared<ChFunction_Sequence>();
    auto my_motion_function21 = chrono_types::make_shared<ChFunction_Ramp>(0, 0.0816496581);
    auto my_motion_function22 = chrono_types::make_shared<ChFunction_Ramp>(0, -0.0816496581);
    auto my_motion_function23 = chrono_types::make_shared<ChFunction_Ramp>(0, 0.0816496581);
    auto my_motion_function24 = chrono_types::make_shared<ChFunction_Ramp>(0, -0.0816496581);
    f_sequence2->InsertFunct(my_motion_function21, 0.2, 0, true);
    f_sequence2->InsertFunct(my_motion_function22, 2.2, 0, true);
    f_sequence2->InsertFunct(my_motion_function23, 2.2, 0, true);
    f_sequence2->InsertFunct(my_motion_function24, 0.2, 0, true);

    auto f_sequence3 = chrono_types::make_shared<ChFunction_Sequence>();
    auto my_motion_function31 = chrono_types::make_shared<ChFunction_Ramp>(0, -0.0471404521);
    auto my_motion_function32 = chrono_types::make_shared<ChFunction_Ramp>(0, 0.0471404521);
    auto my_motion_function33 = chrono_types::make_shared<ChFunction_Ramp>(0, -0.0471404521);
    auto my_motion_function34 = chrono_types::make_shared<ChFunction_Ramp>(0, 0.0471404521);
    f_sequence3->InsertFunct(my_motion_function31, 0.2, 0, true);
    f_sequence3->InsertFunct(my_motion_function32, 2.2, 0, true);
    f_sequence3->InsertFunct(my_motion_function33, 2.2, 0, true);
    f_sequence3->InsertFunct(my_motion_function34, 0.2, 0, true);

    auto f_sequence4 = chrono_types::make_shared<ChFunction_Sequence>();
    auto my_motion_function41 = chrono_types::make_shared<ChFunction_Ramp>(0, -0.03333333335);
    auto my_motion_function42 = chrono_types::make_shared<ChFunction_Ramp>(0, 0.03333333335);
    auto my_motion_function43 = chrono_types::make_shared<ChFunction_Ramp>(0, -0.03333333335);
    auto my_motion_function44 = chrono_types::make_shared<ChFunction_Ramp>(0, 0.03333333335);
    f_sequence4->InsertFunct(my_motion_function41, 0.2, 0, true);
    f_sequence4->InsertFunct(my_motion_function42, 2.2, 0, true);
    f_sequence4->InsertFunct(my_motion_function43, 2.2, 0, true);
    f_sequence4->InsertFunct(my_motion_function44, 0.2, 0, true);

    auto f_sequence5 = chrono_types::make_shared<ChFunction_Sequence>();
    auto my_motion_function51 = chrono_types::make_shared<ChFunction_Ramp>(0, 0.09428090415);
    auto my_motion_function52 = chrono_types::make_shared<ChFunction_Ramp>(0, -0.09428090415);
    auto my_motion_function53 = chrono_types::make_shared<ChFunction_Ramp>(0, 0.09428090415);
    auto my_motion_function54 = chrono_types::make_shared<ChFunction_Ramp>(0, -0.09428090415);
    f_sequence5->InsertFunct(my_motion_function51, 0.2, 0, true);
    f_sequence5->InsertFunct(my_motion_function52, 2.2, 0, true);
    f_sequence5->InsertFunct(my_motion_function53, 2.2, 0, true);
    f_sequence5->InsertFunct(my_motion_function54, 0.2, 0, true);


    auto f_sequence6 = chrono_types::make_shared<ChFunction_Sequence>();
    auto my_motion_function61 = chrono_types::make_shared<ChFunction_Ramp>(0, -0.0816496581);
    auto my_motion_function62 = chrono_types::make_shared<ChFunction_Ramp>(0, 0.0816496581);
    auto my_motion_function63 = chrono_types::make_shared<ChFunction_Ramp>(0, -0.0816496581);
    auto my_motion_function64 = chrono_types::make_shared<ChFunction_Ramp>(0, 0.0816496581);
    f_sequence6->InsertFunct(my_motion_function61, 0.2, 0, true);
    f_sequence6->InsertFunct(my_motion_function62, 2.2, 0, true);
    f_sequence6->InsertFunct(my_motion_function63, 2.2, 0, true);
    f_sequence6->InsertFunct(my_motion_function64, 0.2, 0, true);
    */
    


    auto f_xyz1 = chrono_types::make_shared<ChFunctionPositionXYZFunctions>();
    f_xyz1->SetFunctionX(f_sequence1);
    //f_xyz1->SetFunctionY(f_sequence3);
    //f_xyz1->SetFunctionZ(f_sequence4);

    /*
    * auto f_xyz2 = chrono_types::make_shared<ChFunctionPosition_XYZfunctions>();
    f_xyz2->SetFunctionX(f_sequence2);
    f_xyz2->SetFunctionY(f_sequence3);
    f_xyz2->SetFunctionZ(f_sequence4);

    auto f_xyz3 = chrono_types::make_shared<ChFunctionPosition_XYZfunctions>();
    f_xyz3->SetFunctionY(f_sequence5);
    f_xyz3->SetFunctionZ(f_sequence4);

    auto f_xyz4 = chrono_types::make_shared<ChFunctionPosition_XYZfunctions>();
    f_xyz4->SetFunctionZ(f_sequence1);
    */
    /*
    * auto f_sequence0 = chrono_types::make_shared<ChFunctionSequence>();
    auto my_motion_function0 = chrono_types::make_shared<ChFunctionRamp>(0, 0);
    f_sequence0->InsertFunct(my_motion_function0, 5.0, 0, true);
    auto f_xyz0 = chrono_types::make_shared<ChFunctionPositionXYZfunctions>();
    f_xyz0->SetFunctionX(f_sequence0);
    */
    
    std::cout << "f_sequence1 "<< "\n";
    ChQuaternion<> q(1,0,0,0);
    //q.Q_from_AngZ(atan2(1,2));
    std::vector< std::shared_ptr<ChLinkMotionImposed> > const_list;
    
    auto impose_1 = chrono_types::make_shared<ChLinkMotionImposed>();
    impose_1->Initialize(loaded_nodes[0], mtruss, ChFrame<>(loaded_nodes[0]->GetPos(), q));
    //impose_1->SetPositionFunction(f_xyz0);
    impose_1->SetConstrainedCoords(true, true, true, true, true, true);
    sys.Add(impose_1);
    const_list.push_back(impose_1);
    
    auto impose_2 = chrono_types::make_shared<ChLinkMotionImposed>();
    impose_2->Initialize(loaded_nodes[1], mtruss, ChFrame<>(loaded_nodes[1]->GetPos(), q));
    impose_2->SetPositionFunction(f_xyz1);
    impose_2->SetConstrainedCoords(true, true, true, true, true, true);
    sys.Add(impose_2);
    const_list.push_back(impose_2);

    auto impose_3 = chrono_types::make_shared<ChLinkMotionImposed>();
    impose_3->Initialize(loaded_nodes[2], mtruss, ChFrame<>(loaded_nodes[2]->GetPos(), q));
    //impose_3->SetPositionFunction(f_xyz0);
    impose_3->SetConstrainedCoords(true, true, true, true, true, true);
    sys.Add(impose_3);
    const_list.push_back(impose_3);

    auto impose_4 = chrono_types::make_shared<ChLinkMotionImposed>();
    impose_4->Initialize(loaded_nodes[3], mtruss, ChFrame<>(loaded_nodes[3]->GetPos(), q));
    //impose_4->SetPositionFunction(f_xyz0);
    impose_4->SetConstrainedCoords(true, true, true, true, true, true);
    sys.Add(impose_4);
    const_list.push_back(impose_4);
    
    
	
    double timestep = 1.0E-3; 
     // reduce timestep for better convergence when using Euler implicit linearized
    int stepnum=0;
    
    while (sys.GetChTime() < 4.8) {
		
        //vis->BeginScene();
        //vis->Render();  
        //vis->EndScene();          
        std::cout << "\n/////////////////////////////////////////////////////////////////////////////\nStepnum " << stepnum << "\n";
        sys.DoStepDynamics(timestep); 
        //timestep = 0.0005;
        
        if (stepnum % 1 == 0) {
       
        //if(stepnum%1==0) {
    	//std::string mesh_filename="deneme"+std::to_string(stepnum)+".vtk";
    	//std::string vtk_filename="Vtkdeneme"+std::to_string(stepnum)+".vtk";
    	//WriteMesh(my_mesh, mesh_filename);
    	//WriteFrame(my_mesh, mesh_filename, vtk_filename);
        
        std::cout << "t=\t" << sys.GetChTime() << "\n";
        histfile << "t=\t" << sys.GetChTime() << "\n";
        for (int i=0; i<loaded_nodes.size(); i++){
            ChVector3d node_disp=loaded_nodes[i]->GetPos()-loaded_nodes[i]->GetX0().GetPos();
            ChVector3d node_force=const_list[i]->GetReaction1().force;
            ChVector3d node_moment = const_list[i]->GetReaction1().torque;
        std::cout << "node-"<< i << " disp\t" << node_disp.x()<<"\t"<< node_disp.y()<<"\t" << node_disp.z() << "\n";

        auto nodal_q = loaded_nodes[i]->Frame().GetRot();
        ChVector3d delta_rot_dir;
        double delta_rot_angle;
        nodal_q.GetAngleAxis(delta_rot_angle, delta_rot_dir);
        if (delta_rot_angle > CH_PI)
            delta_rot_angle -= CH_PI/2;  // no 0..360 range, use -180..+180	
        ChVector3d nodal_rot = delta_rot_angle * delta_rot_dir.eigen();
        std::cout << "node-" << i << " rot\t" << nodal_rot[0] << "\t" << nodal_rot[1] << "\t" << nodal_rot[2] << "\n" ;
        std::cout << "react\t" << node_force.x() << "\t" << node_force.y() << "\t" << node_force.z() << "\t" << node_moment.x() << "\t" << node_moment.y() << "\t" << node_moment.z() << "  \n";
        //GetLog() << " t=" << sys.GetChTime() << "  loaded_node rot.x()=" << nodal_rot[0] << "  \n";
        //GetLog() << " t=" << sys.GetChTime() << "  loaded_node rot.y()=" << nodal_rot[1] << "  \n";
        //GetLog() << " t=" << sys.GetChTime() << "  loaded_node rot.z()=" << nodal_rot[2] << "  \n";

        
        /*
        for (int i = 0; i < my_mesh->GetNelements(); ++i) {
            auto elem = std::dynamic_pointer_cast<ChElementLDPM>(my_mesh->GetElement(i));

            //std::cout << " V" << V << std::endl;
            ChMatrixNM<double, 24, 24> M;
            M.setZero();
            elem->ComputeMmatrixGlobal(M);
            ChMatrixDynamic<> SM = elem->GetStiffnessMatrix();

            std::cout << " MassMatrix\t" << M << std::endl;
            std::cout << " StiffnessMatrix\t" << SM << std::endl;

        }
        */
        

        histfile << "node-" << i << " disp\t" << node_disp.x() << "\t" << node_disp.y() << "\t" << node_disp.z() << "\t" << nodal_rot[0] << "\t" << nodal_rot[1] << "\t" << nodal_rot[2] << "\t" << "react\t" << node_force.x() << "\t" << node_force.y() << "\t" << node_force.z() << "\t" << node_moment.x() << "\t" << node_moment.y() << "\t" << node_moment.z() << "\n";
	}
    for (int i = 0; i < my_mesh->GetNumElements(); ++i) {
                auto elem = std::dynamic_pointer_cast<ChElementLDPM>(my_mesh->GetElement(i));
                for (auto facet : elem->GetSection()) {
                    auto statev = facet->Get_StateVar();
                    std::cout << " strain\t" << statev.segment(0,3).transpose() <<" stress\t" << statev.segment(3,3).transpose() << "\n";
                    histfile <<  statev.segment(3,3).transpose() << " " << statev.segment(0,3).transpose() << " ";

                }
            }
        
    std::cout << "\n";
	histfile << "\n";
       // }
        
    }        
        stepnum++;
    }
    
     histfile.close();
    
	/*
    while (vis->Run()) {
        vis->BeginScene();
        vis->Render();
        vis->EndScene();
        //sys.DoStaticLinear();
        //sys.DoStepDynamics(0.01);           
    }
    */
    
    
    
  }  
	
	std::cout<<"\nProgram completed succesfully!!!!!!!!!!!!\n\n";
    return 0;
}
        
