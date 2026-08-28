
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

						  
						  
						   
						  
										 


//#include "chrono_ldpm/ChApiLDPM.h"
#include "chrono_ldpm/ChElementLDPM.h"
#include "chrono_ldpm/ChBuilderLDPM.h"
#include "chrono_ldpm/ChMaterialVECT.h"
#include "chrono_ldpm/ChSectionLDPM.h"

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

#include <chrono/fea/ChMeshExporter.h>

#include "chrono/fea/ChContactSurfaceMesh.h"
#include "chrono/fea/ChContactSurfaceNodeCloud.h"

#include "chrono/physics/ChLinkMotorLinearPosition.h"
#include "chrono/physics/ChLinkMotorRotationSpeed.h"


#include "chrono/functions/ChFunctionPositionXYZFunctions.h"
#include "chrono/physics/ChLinkMotionImposed.h"

#include <filesystem>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <typeinfo>
#include <map>
#include <array>
#include <set>


using namespace chrono;
//using namespace chrono::geometry;
using namespace chrono::fea;
//using namespace chrono::irrlicht;
using namespace chrono::ldpm;
//using namespace irr;

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
            /*} else if (auto elementBm = std::dynamic_pointer_cast<ChElementCSL>(mesh->GetElement(iele)))  {
                mynodes.resize(2);
                out_stream << "2 ";
                int nodeOrder[] = {0, 1};
                mynodes[0] = elementBm->GetNode(nodeOrder[0])->GetIndex();
                mynodes[1] = elementBm->GetNode(nodeOrder[1])->GetIndex();
                BeamElemNodes.push_back(mynodes);
                for (int myNodeN = 0; myNodeN < mynodes.size(); myNodeN++) {
                    auto nodeA = (elementBm->GetNode(nodeOrder[myNodeN]));
                    std::vector<std::shared_ptr<ChNodeFEAbase>>::iterator it;
                    it = find(myvector.begin(), myvector.end(), nodeA);
                    if (it == myvector.end()) {
                        // name not in vector
                    } else {
                        auto index = std::distance(myvector.begin(), it);
                        out_stream << (unsigned int)index << " ";
                    }
                }
                out_stream << "\n";*/
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

    std::vector<std::vector<int>> FacetNodes;

    /*
    std::vector<std::shared_ptr<ChNodeFEAbase>> myvector;
    myvector.resize(mesh->GetNumNodes());

    for (unsigned int i = 0; i < mesh->GetNumNodes(); i++) {
        myvector[i] = std::dynamic_pointer_cast<ChNodeFEAbase>(mesh->GetNode(i));
    }

    */
    
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

        FacetNodes.push_back(mynodes);

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
            // auto pC = vertices[0]->GetPos();
            // auto pA = vertices[1]->GetPos();    
            // auto pB = vertices[2]->GetPos();
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



    /*
    out_stream << "\nCELL_DATA " << numCell*12 << "\n";
    out_stream << "VECTORS Crack float\n";
    //out_stream << "LOOKUP_TABLE default\n";

    for (unsigned int iele = 0; iele < mesh->GetNumElements(); iele++) {

        auto elem = std::dynamic_pointer_cast<ChElementLDPM>(mesh->GetElement(iele));

        for (auto facet : elem->GetSection()) {
            auto statev = facet->Get_StateVar();

            out_stream << statev(0) << " " << statev(1) << " " << statev(2) << "\n";
        }

    }
    */


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




ChVector3d calculate_Force(std::vector< std::shared_ptr<ChLinkMateGeneric> > const_list ){
    unsigned int icons=0;   
    ChVector3d Force(0.,0.,0.); 
    for (auto constraint:const_list) { 
        	auto fn=constraint->GetReaction1(); 
        	//std::cout<<"fn "<<fn<<"\n";   
        	Force+=fn.force;
    }
    
    return Force;

}


class UserContactMaterial : public ChContactContainer::AddContactCallback {
public:
     virtual void OnAddContact(const ChCollisionInfo& contactinfo, ChContactMaterialComposite* const material) override {
													   
        // Downcast to appropriate composite material type
         auto mat = static_cast<ChContactMaterialCompositeNSC* const>(material);
        ChVector3d relvel = contactinfo.vpB - contactinfo.vpA;
        auto normal_dir = contactinfo.vN;
        double relvel_n_mag = relvel.Dot(normal_dir);
        ChVector3d relvel_n = relvel_n_mag * normal_dir;
        ChVector3d relvel_t = relvel - relvel_n;
        double relvel_t_mag = relvel_t.Length();
        double dT = sys->GetStep();
        double delta_t = relvel_t_mag * dT;
        //
        double mu = mud + (mus - mud) * s0 / (s0 + delta_t);
        mat->static_friction = mu;
        mat->sliding_friction = mud;
    }

    void SetSystem(ChSystem& mysystem) { sys = &mysystem; }
    double s0 = 1.3; //mm
    double mus = 0.13; //static friction
    double mud = 0.015; //dynamic friction
    ChSystem* sys;
};

// Compute triangular facet area using cross product
double TriangleArea(const ChVector3d& a, const ChVector3d& b, const ChVector3d& c) {
    return 0.5 * ((b - a).Cross(c - a)).Length();
}

struct BoundaryFace {
    std::shared_ptr<ChNodeFEAxyzrot> n1;
    std::shared_ptr<ChNodeFEAxyzrot> n2;
    std::shared_ptr<ChNodeFEAxyzrot> n3;
};



// Compute tetrahedron volume from its four nodes
double TetraVolume(const ChVector3d& a,
                   const ChVector3d& b,
                   const ChVector3d& c,
                   const ChVector3d& d) {
    return fabs((b - a).Dot((c - a).Cross(d - a))) / 6.0;
}



int main(int argc, char** argv) {
     SetChronoDataPath(CHRONO_DATA_DIR);
    //
    // Create a Chrono::Engine physical system
    ChSystemSMC sys;
	
    sys.SetNumThreads(24);
	// Collision system type
    //auto collision_type = collision::ChCollisionSystemType::BULLET;
    //sys.SetCollisionSystemType(collision_type);
    //sys.SetMaxPenetrationRecoverySpeed(0.01);  						
    //
    //    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    ///
    /// Give the path and file name of LDPM data
    ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    std::string current_dir(argv[0]);
    
    current_dir="";  
    //  
    std::string LDPM_data_path=current_dir+"LDPMgeo000Cylinder000/";
    std::string LDPM_GeoName="LDPMgeo000";
    //
    std::string out_dir = "out/";

if (!std::filesystem::exists(out_dir)) {
    std::filesystem::create_directory(out_dir);
}
    	
    std::string history_filename="hist.dat"; 
    std::string num_iter_filename="num_iter.dat";
    std::ofstream histfile;
    histfile.open(out_dir+history_filename, std::ios::out);
    std::ofstream num_iter_file;
    num_iter_file.open(out_dir+num_iter_filename, std::ios::out);
    //
    //	
    // Create ground:   
    // 
    auto mtruss = chrono_types::make_shared<ChBody>();
    mtruss->SetFixed(true);  
    //mtruss->SetCollide(false);	
    sys.Add(mtruss);
    //    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    ///
    /// Define loading plates
    ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    //
 

    //auto mysurfmaterial = chrono_types::make_shared<ChMaterialSurfaceNSC>();    
    //mysurfmaterial->SetFriction(0.5f);   

    ///
    /*
        auto bottom_plate = chrono_types::make_shared<ChBodyEasyBox>(140,140,10,7.8E-9,true,true,mysurfmaterial);
    bottom_plate->SetPos(ChVector<>(50,50,-5));
    //bottom_plate->SetCollide(true);
    bottom_plate->GetVisualShape(0)->SetTexture("/chrono-concrete/data/textures/blue.png");
    //bottom_plate->SetBodyFixed(true);
    sys.Add(bottom_plate);    
    
	auto constr_bot_plate=chrono_types::make_shared<ChLinkMateGeneric>(true, true, true, true, true, true);
    constr_bot_plate->Initialize(bottom_plate, mtruss, bottom_plate->GetFrame_COG_to_abs()); 
    sys.Add(constr_bot_plate); 
    ///    
    ///    
    auto top_plate = chrono_types::make_shared<ChBodyEasyBox>(140, 140, 10,7.8E-9,true,true,mysurfmaterial);
    top_plate->SetPos(ChVector<>(50,50,105)); 
    top_plate->GetVisualShape(0)->SetTexture("/chrono-concrete/data/textures/blue.png");  
    sys.Add(top_plate); 
    */


    /// 
    //auto constr_top_plate = chrono_types::make_shared<ChLinkMateGeneric>(true, true, false, true, true, true);
    //constr_top_plate->Initialize(top_plate, mtruss, top_plate->GetFrame_COG_to_abs());
    //sys.Add(constr_top_plate);

    //auto constr_top_plate = chrono_types::make_shared<ChLinkLockLock>();
    //constr_top_plate->Initialize(top_plate, mtruss, ChCoordsys<>(top_plate->GetPos()));
    //sys.Add(constr_top_plate);


    //    
    //auto constr_top_plate = chrono_types::make_shared<ChLinkLockLock>();
    //constr_top_plate->Initialize(top_plate, mtruss, ChCoordsys<>(top_plate->GetPos()));
    //sys.Add(constr_top_plate);

	/*
    auto constr_top_plate = chrono_types::make_shared<ChLinkLockPrismatic>();    
    constr_top_plate->Initialize(top_plate, mtruss, ChCoordsys<>( ChCoordsys<>(top_plate->GetPos(), Q_from_AngAxis(CH_C_PI_2, ChVector<>(1,0,0))))); 
    sys.Add(constr_top_plate);	
    */
    /*
    auto constr_top_plate = chrono_types::make_shared<ChLinkLockLock>();    
    constr_top_plate->Initialize(top_plate, mtruss, ChCoordsys<>(ChVector<>(0, 0, 0))); 
    sys.Add(constr_top_plate);	
    */
    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    ///
    /// Read LDPM Freecad outputs and insert into mesh object
    ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    //
    // Create a mesh object
    //
    auto my_mesh = chrono_types::make_shared<ChMesh>();
    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    ///
    /// Read LDPM Freecad outputs and insert into mesh object
    ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    //
    // Create  vectorial material for LDPM
    //
    auto my_mat = chrono_types::make_shared<ChMaterialVECT>();
    my_mat->Set_density(2.338E-9);          // Density
    my_mat->Set_E0(60273);                // Normal modulus
    my_mat->Set_alpha(0.25);            // Alpha
    my_mat->Set_sigmat(3.44);            // Tensile strength
    my_mat->Set_sigmas(8.944);               // shear strength
    my_mat->Set_nt(0.4);                     // Softening exponent
    my_mat->Set_lt(500);                     // Tensile characteristic length
    my_mat->Set_Ed(60273);                  // Ed
    my_mat->Set_sigmac0(150);              // Compressive yielding strength
    my_mat->Set_Hc0(24109);                // 𝐻𝑐0 = Initial hardening modulus ratio * Normal modulus
    my_mat->Set_Hc1(6027.3);               // 𝐻𝑐1= Final hardening modulus ratio * Normal modulus
    my_mat->Set_beta(0);                 // Volumetric deviatoric coupling
    my_mat->Set_kc0(4);                 // Transitional strain ratio
    my_mat->Set_kc1(1);                 // Deviatoric strain threshold ratio
    my_mat->Set_kc2(5);                  // Deviatoric damage parameter
    my_mat->Set_kc3(0.1);               // Volumetric strain parameter
    my_mat->Set_mu0(0.4);               // Initial friction
    my_mat->Set_muinf(0);               // Asymptotic friction
    my_mat->Set_sigmaN0(600);            // Transitional stress
    //
    //read_LDPM_info(my_mesh, nodesFilename, elemFilename, facetFilename, tetsFilename, verticesFilename);
    ChBuilderLDPM builder;
    builder.read_LDPM_info(my_mesh, my_mat, LDPM_data_path, LDPM_GeoName);
    //	
    sys.Add(my_mesh);   
  
    double xmin = 1e20, xmax = -1e20;
    double ymin = 1e20, ymax = -1e20;
    double zmin = 1e20, zmax = -1e20;

for (unsigned int i = 0; i < my_mesh->GetNumNodes(); i++) {
    auto node = std::dynamic_pointer_cast<ChNodeFEAxyzrot>(my_mesh->GetNode(i));
    auto p = node->GetX0().GetPos();

    xmin = std::min(xmin, p.x());
    xmax = std::max(xmax, p.x());

    ymin = std::min(ymin, p.y());
    ymax = std::max(ymax, p.y());

    zmin = std::min(zmin, p.z());
    zmax = std::max(zmax, p.z());
}

   double center_x = 0.5 * (xmin + xmax);
   double center_y = 0.5 * (ymin + ymax);
   double center_z = 0.5 * (zmin + zmax);

   double D = std::max(xmax - xmin, ymax - ymin);
   double R = D / 2.0;
   double H = zmax - zmin;

   double A_specimen = 3.14159265359 * R * R;
   double A_lateral = 2.0 * 3.14159265359 * R * H;

   double tol_z = 0.001 * H;
   double tol_r = 0.02 * R;  

   std::cout << "Detected geometry:\n";
   std::cout << "D = " << D << "\n";
   std::cout << "R = " << R << "\n";
   std::cout << "H = " << H << "\n";
   std::cout << "zmin = " << zmin << "\n";
   std::cout << "zmax = " << zmax << "\n";


        //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    ///
    /// Pick the center nodes on top and bottom
    ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    //
    //loaded_node_1 = hnode2;    
    //loaded_node_1->SetForce(F_node_1);
    // Apply constrain on a group of nodes

    std::shared_ptr<ChNodeFEAxyzrot> RP1;
    //RP1->SetX0(ChFrame<>(50, 50, 200));
    std::shared_ptr<ChNodeFEAxyzrot> RP2;
    //RP2->SetX0(ChFrame<>(50, 50, 0));

    for (unsigned int i = 0; i < my_mesh->GetNumNodes(); i++) {
        auto node = std::dynamic_pointer_cast<ChNodeFEAxyzrot>(my_mesh->GetNode(i));
        auto p = node->GetPos();
        if (abs(p.z() - zmin) < 0.1 &&
    abs(p.x() - center_x) < 2 &&
    abs(p.y() - center_y) < 2) {
        RP2 = node;
                    }
    }

    for (unsigned int i = 0; i < my_mesh->GetNumNodes(); i++) {
        auto node = std::dynamic_pointer_cast<ChNodeFEAxyzrot>(my_mesh->GetNode(i));
        auto p = node->GetPos();
        if (abs(p.z() - zmax) < 0.1 &&
    abs(p.x() - center_x) < 2 &&
    abs(p.y() - center_y) < 2) {
        RP1 = node;
                           }
    }

        //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    ///
    /// Select the nodes on the top surface of the concrete cube
    ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    //
    //loaded_node_1 = hnode2;    
    //loaded_node_1->SetForce(F_node_1);
    // Apply constrain on a group of nodes
    std::vector< std::shared_ptr<ChNodeFEAxyzrot> > top_nodes;	    
    for (unsigned int i = 0; i < my_mesh->GetNumNodes(); i++) {
        auto node = std::dynamic_pointer_cast<ChNodeFEAxyzrot>(my_mesh->GetNode(i)); 
        auto p=node->GetPos();
        if (abs(p.z() - zmax) < tol_z && node.get() != RP1.get()) {
        top_nodes.push_back(node);
        }        	
        	//Stop_surface->AddNode(node,0.025);
        	//top_surf->AddNode(node);
        	//node->SetForce(F_node_1);
        	//node->SetPos(U_node);
        }       
    
    //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    ///
    /// Select the nodes on the bottom surface of the concrete cube
    ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                       
    // Apply constraint on a group of nodes
    //	
    std::vector< std::shared_ptr<ChNodeFEAxyzrot> > bottom_nodes;
    //std::vector< std::shared_ptr<ChLinkMateGeneric> > const_list;
    //auto constr_bc=chrono_types::make_shared<ChLinkMateGeneric>(true, true, true, true, true, true); 
    unsigned int icons=0;
    for (unsigned int i = 0; i < my_mesh->GetNumNodes(); i++) {
        auto node = std::dynamic_pointer_cast<ChNodeFEAxyzrot>(my_mesh->GetNode(i)); 
        auto p=node->GetPos();
        if (abs(p.z() - zmin) < tol_z && node.get() != RP2.get()) {
        bottom_nodes.push_back(node);
        }       
    }
	

    std::vector< std::shared_ptr<ChNodeFEAxyzrot> > sur_nodes;
    //std::vector< std::shared_ptr<ChLinkMateGeneric> > const_list;
    //auto constr_bc=chrono_types::make_shared<ChLinkMateGeneric>(true, true, true, true, true, true); 
    //unsigned int icons = 0;

    
    for (unsigned int i = 0; i < my_mesh->GetNumNodes(); i++) {
        auto node = std::dynamic_pointer_cast<ChNodeFEAxyzrot>(my_mesh->GetNode(i));
        auto p = node->GetPos();

        double dx = p.x() - center_x;
        double dy = p.y() - center_y;
        double r = sqrt(dx * dx + dy * dy);

        if (r > R - tol_r) {
            sur_nodes.push_back(node);
        }
    }




    //std::cout << "RP1" << RP1->GetX0() << std::endl;
    //std::cout << "RP2" << RP2->GetX0() << std::endl;
    std::cout << "TOP nodes number:" << top_nodes.size() << std::endl;
    std::cout << "BOT nodes number:" << bottom_nodes.size() << std::endl;
    std::cout << "Surround nodes number:" << sur_nodes.size() << std::endl;

    auto constr_bot = chrono_types::make_shared<ChLinkMateGeneric>(true, true, true, true, true, true);
    constr_bot->Initialize(RP2, mtruss, RP2->Frame());
    sys.Add(constr_bot);

/*
    std::vector< std::shared_ptr<ChLinkMateGeneric> > constr_sur_list;
    //unsigned int icons = 0;
    for (auto node : sur_nodes) {
        auto constr_sur = chrono_types::make_shared<ChLinkMateGeneric>(true, true, false, false, false, false);
        constr_sur->Initialize(node, mtruss, node->Frame());
        //constr_sur->Initialize(node, mtruss, false, node->Frame(), node->Frame());
        sys.Add(constr_sur);
        constr_sur_list.push_back(constr_sur);
    }
*/
    ///
    ///
    /// Tie bottom nodes to bottom plate 
    ///  
    /// 

    
        std::vector< std::shared_ptr<ChLinkMateGeneric> > const_bot;
    for (auto node: bottom_nodes) {
        auto constr_tie=chrono_types::make_shared<ChLinkMateGeneric>(false, false, true, true, true, false);
        constr_tie->Initialize(node, RP2, false, node->Frame(), node->Frame());
        const_bot.push_back(constr_tie);
        sys.Add(constr_tie);

    }
    
    // Create a rigid top loading plate reference body
auto top_plate = chrono_types::make_shared<ChBody>();
top_plate->SetMass(1e-6);
top_plate->SetInertiaXX(ChVector3d(1e-6, 1e-6, 1e-6));
top_plate->SetPos(ChVector3d(center_x, center_y, zmax));   // at top surface
sys.Add(top_plate);


//////////////////////////////////////////////////////////
// External vertical force applied to the top plate
//////////////////////////////////////////////////////////
auto top_plate_force = chrono_types::make_shared<ChForce>();
top_plate->AddForce(top_plate_force);

top_plate_force->SetMode(ChForce::FORCE);
top_plate_force->SetDir(ChVector3d(0, 0, -1));
top_plate_force->SetVrelpoint(ChVector3d(0, 0, 0));  



// Top plate boundary condition:
auto constr_top_plate = chrono_types::make_shared<ChLinkMateGeneric>(true, true, false, true, true, true);
constr_top_plate->Initialize(top_plate, mtruss, top_plate->GetFrameCOMToAbs());
sys.Add(constr_top_plate);


//////////////////////////////////////
/// Tie top nodes to top plate 
////////////////////////////////////// 
    
    std::vector< std::shared_ptr<ChLinkMateGeneric> > const_top;
    for (auto node: top_nodes) {
        auto constr_tie=chrono_types::make_shared<ChLinkMateGeneric>(false, false, true, true, true, false);
        constr_tie->Initialize(node, top_plate, false, node->Frame(), node->Frame());
        const_top.push_back(constr_tie);
        sys.Add(constr_tie);
    }
    
// Top RP boundary condition:
// Fix X, Y and all rotations; keep Z free
/*auto constr_top_RP = chrono_types::make_shared<ChLinkMateGeneric>(true, true, false, true, true, true);
constr_top_RP->Initialize(RP1, mtruss, RP1->Frame());
sys.Add(constr_top_RP); */ 
    




 
    /// traced node
    /// 
    ///
    /*
    auto traced_node= chrono_types::make_shared<ChNodeFEAxyzrot>();
    for (unsigned int i = 0; i < my_mesh->GetNumNodes(); i++) {
        auto node = std::dynamic_pointer_cast<ChNodeFEAxyzrot>(my_mesh->GetNode(i)); 
        double cy=node->GetPos().y();
        double cz=node->GetPos().z();
        double cx=node->GetPos().x();   
        if (cy>49. && abs(cx-25)<2 && abs(cz-25)<2 ){            		
        	traced_node=node;
        	break;
        }       
    }
    std::cout << "traced_node->GetIndex() : "<<traced_node->GetIndex()<<std::endl;
    */

    
    
	/*
    auto mcontactsurf = chrono_types::make_shared<ChContactSurfaceMesh>(mysurfmaterial);
    my_mesh->AddContactSurface(mcontactsurf);
    mcontactsurf->AddFacesFromBoundary(1);    
    */
    
    /*
        auto mcontactcloud = chrono_types::make_shared<ChContactSurfaceNodeCloud>(mysurfmaterial);
    my_mesh->AddContactSurface(mcontactcloud);
    mcontactcloud->AddAllNodes(0.0);
    */

    


    
    /*
   ChVector <> F_node_1(0, 1000, 0); 
    for (auto node: top_nodes) {
    	//auto link_1 = chrono_types::make_shared<ChLinkLockLock>();
    	//sys.Add(link_1);
    	//link_1->Initialize(node, mtruss);
		//node->SetForce(F_node_1);
    	auto impose_1 = chrono_types::make_shared<ChLinkMotionImposed>();    	
    	impose_1->Initialize(node, mtruss, ChFrame<>(node->GetPos()));
    	impose_1->SetPositionFunction(f_xyz);
		sys.Add(impose_1);
    	
    } 
    */   
     
    // We do not want gravity effect on FEA elements in this demo
    my_mesh->SetAutomaticGravity(false);
    //sys.Set_G_acc(ChVector<>(0, 0, 0));
    //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    ///
    /// Create a visualization system
    ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    //
    /*
    auto mvisualizebeamA = chrono_types::make_shared<ChVisualShapeFEA>(my_mesh);
    //mvisualizebeamA->SetFEMdataType(ChVisualShapeFEA::DataType::ELEM_BEAM_MZ);   
    mvisualizebeamA->SetFEMdataType(ChVisualShapeFEA::DataType::NODE_DISP_NORM); 
    mvisualizebeamA->SetColorscaleMinMax(-5., 0.);
    mvisualizebeamA->SetSmoothFaces(true);
    mvisualizebeamA->SetWireframe(false);
    my_mesh->AddVisualShapeFEA(mvisualizebeamA);


    auto mvisualizebeamC = chrono_types::make_shared<ChVisualShapeFEA>(my_mesh);
    mvisualizebeamC->SetFEMglyphType(ChVisualShapeFEA::GlyphType::NODE_DOT_POS);
    mvisualizebeamC->SetFEMdataType(ChVisualShapeFEA::DataType::NONE);
    mvisualizebeamC->SetSymbolsThickness(1);
    mvisualizebeamC->SetSymbolsScale(1);
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
    vis->AddCamera(ChVector3d(100.0, 100.0, 200.));
    vis->AttachSystem(&sys);
    */
    


    
    

	
    
    	/*
    auto mvisualizemeshcoll = chrono_types::make_shared<ChVisualShapeFEA>(my_mesh);
    mvisualizemeshcoll->SetFEMdataType(ChVisualShapeFEA::DataType::NODE_DISP_NORM);
    mvisualizemeshcoll->SetWireframe(true);
    mvisualizemeshcoll->SetDefaultMeshColor(ChColor(1, 0.5, 0));
    my_mesh->AddVisualShapeFEA(mvisualizemeshcoll);
    	*/
	
	

    //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    ///
    /// Create a Chrono solver and set solver settings
    ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    
    ///
    ///
    /// Select solver type
    ///
    ///
    // Use the ADMM solver: it has the capability of handling both FEA and NSC!
    /*
    sys.SetSolverForceTolerance(1e-10);
    //auto solver = chrono_types::make_shared<ChSolverADMM>(chrono_types::make_shared<ChSolverSparseLU>());
    auto solver = chrono_types::make_shared<ChSolverADMM>(chrono_types::make_shared<ChSolverPardisoMKL>());
    solver->EnableWarmStart(true);
    solver->SetMaxIterations(50);
    solver->SetToleranceDual(1e-6);
    solver->SetTolerancePrimal(1e-6);
    solver->SetRho(1.);
    solver->SetStepAdjustPolicy(ChSolverADMM::AdmmStepType::BALANCED_UNSCALED);
    solver->SetVerbose(false);
    sys.SetSolver(solver);

    sys.SetSolverForceTolerance(1e-10);
    */


	/*
    //auto solver = chrono_types::make_shared<ChSolverSparseQR>();   
    auto solver = chrono_types::make_shared<ChSolverSparseLU>();
    sys.SetSolver(solver);
    solver->UseSparsityPatternLearner(true);
    solver->LockSparsityPattern(true);
    solver->SetVerbose(false);	
    */

	
    
	
    auto solver = chrono_types::make_shared<ChSolverPardisoMKL>();    
    sys.SetSolver(solver); 
    solver->UseSparsityPatternLearner(true);
    solver->LockSparsityPattern(true);   
    //solver->SetVerbose(true);
	//sys.Update();
    
       
    /*	
    auto solver = chrono_types::make_shared<ChSolverMINRES>();
    sys.SetSolver(solver);
    solver->SetMaxIterations(400);
    solver->SetTolerance(1e-12);
    solver->EnableDiagonalPreconditioner(true);
    solver->EnableWarmStart(true);  // Enable for better convergence when using Euler implicit linearized
    solver->SetVerbose(false);
    sys.SetSolverForceTolerance(1e-14);
	*/  
	
    ///
    ///
    /// Select time stepper
    ///
    ///
	
        //sys.SetTimestepperType(ChTimestepper::Type::HHT);   
        //auto mystepper = std::dynamic_pointer_cast<ChTimestepperHHT>(sys.GetTimestepper());
auto mystepper = chrono_types::make_shared<ChTimestepperHHT>(&sys);
mystepper->SetAlpha(-0.05);
mystepper->SetMaxIters(50);
mystepper->SetAbsTolerances(1e-03, 1e-03);
mystepper->SetMinStepSize(1E-15);
mystepper->SetMaxItersSuccess(4);
mystepper->SetRequiredSuccessfulSteps(3);
mystepper->SetStepIncreaseFactor(1.25);
mystepper->SetStepDecreaseFactor(0.25);
mystepper->SetJacobianUpdateMethod(ChTimestepperHHT::JacobianUpdate::NEVER);
mystepper->SetVerbose(false);
mystepper->SetStepControl(false);
sys.SetTimestepper(mystepper);
          //}
    

    
    
    

    

    
  
    

    	
    	
    //auto mystepper=chrono_types::make_shared<ChTimestepperNewmark>(&sys);
    //mystepper->SetGammaBeta(0.5, 0.25);  // Newmark as const accel. method    
    //mystepper->SetGammaBeta(0.5, 1 / 6);  // Newmark as linear accel. method    
    //mystepper->SetGammaBeta(1.0, 0.25);  // Newmark with max numerical damping
    //auto mystepper=chrono_types::make_shared<ChTimestepperEulerImplicit>(&sys);
    //auto mystepper=chrono_types::make_shared<ChTimestepperEulerImplicitLinearized>(&sys);
    //auto mystepper=chrono_types::make_shared<ChTimestepperEulerImplicitProjected>(&sys);
    //auto mystepper=chrono_types::make_shared<ChTimestepperTrapezoidal>(&sys);
    //auto mystepper=chrono_types::make_shared<ChTimestepperTrapezoidalLinearized>(&sys);
    //auto mystepper=chrono_types::make_shared<ChTimestepperTrapezoidalLinearized2>(&sys);
    //auto mystepper=chrono_types::make_shared<ChTimestepperLeapfrog>(&sys);
    //auto mystepper=chrono_types::make_shared<ChTimestepperHeun>(&sys);
    //auto mystepper=chrono_types::make_shared<ChTimestepperRungeKuttaExpl>(&sys);
    //auto mystepper=chrono_types::make_shared<ChTimestepperEulerSemiImplicit>(&sys);
    //auto mystepper=chrono_types::make_shared<ChTimestepperEulerExplIIorder>(&sys);
    //auto mystepper=chrono_types::make_shared<ChTimestepperEulerExpl>(&sys);    
    //sys.SetTimestepper(mystepper);
    ///
    ///
    
    //sys.SetTimestepperType(ChTimestepper::Type::EULER_IMPLICIT_LINEARIZED);
    //sys.SetTimestepperType(ChTimestepper::Type::EULER_IMPLICIT_PROJECTED);
    //sys.SetTimestepperType(ChTimestepper::Type::EULER_IMPLICIT);
    //sys.SetTimestepperType(ChTimestepper::Type::NEWMARK);
    //sys.SetTimestepperType(ChTimestepper::Type::TRAPEZOIDAL);
    //sys.SetTimestepperType(ChTimestepper::Type::RUNGEKUTTA45);
    //sys.SetTimestepperType(ChTimestepper::Type::EULER_EXPLICIT);
    //auto mystepper = std::dynamic_pointer_cast<ChTimestepperEulerExplIIorder>(sys.GetTimestepper());
    //sys.SetTimestepper(mystepper);
    //sys.DoStaticLinear();    
    //sys.DoStaticNonlinear(50);             
    
        
        
        
    //std::string history_filename="hist.dat";   
    //std::ofstream histfile(history_filename);  
       
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    ///
    /// Analysis
    ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    
    /*
    for(int i=0; i<my_mesh->GetNumElements(); i++){    
	    	auto elem = std::dynamic_pointer_cast<ChElementLDPM>(my_mesh->GetElement(i));
	    	elem->ComputeStiffnessMatrix();
	    	//
	    	auto nI=std::dynamic_pointer_cast<ChNodeFEAxyz>(elem->GetNodeN(0));
	    	auto nJ=std::dynamic_pointer_cast<ChNodeFEAxyz>(elem->GetNodeN(1));
	    	auto nK=std::dynamic_pointer_cast<ChNodeFEAxyz>(elem->GetNodeN(2));
	    	auto nL=std::dynamic_pointer_cast<ChNodeFEAxyz>(elem->GetNodeN(3));
	    	std::cout<< i <<". element:  ";
	    	std::cout<<nI->Frame().GetPos().x()<<"\t"<<nI->Frame().GetPos().y()<<"\t"<<nI->Frame().GetPos().z()<<"\t";
	    	std::cout<<nJ->Frame().GetPos().x()<<"\t"<<nJ->Frame().GetPos().y()<<"\t"<<nJ->Frame().GetPos().z()<<"\t";
	    	std::cout<<nK->Frame().GetPos().x()<<"\t"<<nK->Frame().GetPos().y()<<"\t"<<nK->Frame().GetPos().z()<<"\t";
	    	std::cout<<nL->Frame().GetPos().x()<<"\t"<<nL->Frame().GetPos().y()<<"\t"<<nL->Frame().GetPos().z()<<std::endl;
	    	//
	    	
    		auto Km=elem->GetStiffnessMatrix();
    		std::cout<<"Km\n"<< Km <<"\n\n";
    		//
    		
    		 ChMatrixDynamic<> Mloc(24, 24);
        	Mloc.setZero();
    		elem->ComputeMmatrixGlobal(Mloc);
    		std::cout<<"mLOC\n"<< Mloc <<"\n\n";
    		
    		//
    		ChQuaternion<> qa0 = nJ->GetRot();
    		ChVector<> rotator(VNULL);  rotator[0] = 1e-5;
                ChQuaternion<> mdeltarotL;  mdeltarotL.Q_from_Rotv(rotator); // rot.in local basis - as in system wide vectors
                ChQuaternion<> qaD = qa0 * mdeltarotL;
    		nJ->SetRot(qaD);
    		ChVectorDynamic<> displ(elem->GetNdofs());
    		elem->GetStateBlock(displ);
    		std::cout<< "displ:\n"<< displ<<std::endl;
    		//std::cout<< -Km*displ<<std::endl;
    		//
    		ChVectorDynamic<> Fi(24);
    		elem->ComputeInternalForces(Fi);
    		std::cout<<"Fi\n"<< Fi+Km*displ <<"\n\n";
    		
    		exit(9);
    
    	}
        //exit(9);
    */


    		
    
    
    
    if (true) {
	std::cout << "LDPM RESULTS (DYNAMIC ANALYSIS) \n\n";
	///
        /// Displacement controlled 
        ///
        
        /*
	auto motion = chrono_types::make_shared<ChFunction_Ramp>();
    	motion->Set_ang(10.0);      
    	constr_top_plate->SetMotion_Y(motion);  
	*/
	
/*
    auto motor1 = chrono_types::make_shared<ChLinkMotorLinearPosition>();

    //Connect the guide and the slider and add the motor to the system:
    motor1->Initialize(RP1,               // body A (slave)
                       mtruss,                // body B (master)
                       ChFrame<>(RP1->GetPos(), QUNIT)  // motor frame, in abs. coords
    );    

    auto my_motion_function1 = chrono_types::make_shared<ChFunctionPoly>();
	my_motion_function1->SetCoefficients(std::vector<double>{0.0, 0.0, -25000.});

    auto my_motion_function2 = chrono_types::make_shared<ChFunctionRamp>(0, -50);

    auto f_sequence1 = chrono_types::make_shared<ChFunctionSequence>();
    f_sequence1->InsertFunct(my_motion_function1, 0.001, 1, true);
    f_sequence1->InsertFunct(my_motion_function2, 5.0, 1, true);

    motor1->SetMotionFunction(f_sequence1);
    sys.Add(motor1);
*/


    
    /*
        auto my_motion_function1 = chrono_types::make_shared<ChFunction_Poly>();
    my_motion_function1->Set_coeff(0.0, 0);
    my_motion_function1->Set_coeff(0.0, 1);
    my_motion_function1->Set_coeff(-1000, 2);
    my_motion_function1->Set_order(2);
    auto my_motion_function2 = chrono_types::make_shared<ChFunction_Ramp>(0, -2.0);

    auto f_sequence1 = chrono_types::make_shared<ChFunction_Sequence>();
    f_sequence1->InsertFunct(my_motion_function1, 0.001, 1, true);
    f_sequence1->InsertFunct(my_motion_function2, 1.0, 1, true);

    constr_top_plate->SetMotion_Z(f_sequence1);

    ChVector<double> initial_pos = top_plate->GetPos();
    */

    
	
	///
        /// Load controlled 
        ///
        //auto mod = chrono_types::make_shared<ChFunction_Ramp>();
    	//mod->Set_ang(1000000.0);
   	
   	/*
	// Actuate second slider using a body force
	auto frcY = chrono_types::make_shared<ChForce>();
	frcY->SetF_y(mod);
	top_plate->AddForce(frcY); 	
	*/

    
    //Add following lines before while loop (before beginning of dynamic analysis):
    //auto cmaterial = chrono_types::make_shared<UserContactMaterial>();
    //cmaterial->SetSystem(sys);
    //sys.GetContactContainer()->RegisterAddContactCallback(cmaterial);

	double timestep = 1E-4; 
	int stepnum=0;
    
    
    // Hydrostatic test - vertical pressure loading 
    

    double Pmax = 800;                // MPa = N/mm^2
    double T0 = 0.8;                  // loading time sec



    double P_applied = 0.0;
    double F_top = 0.0;


// Collect true boundary faces on the lateral cylindrical surface
std::vector<BoundaryFace> lateral_faces;

std::map<std::array<int, 3>, BoundaryFace> face_map;
std::map<std::array<int, 3>, int> face_count;

// Four triangular faces of a tetrahedron
std::array<std::array<int, 3>, 4> tet_faces = {{
    {{0, 1, 2}},
    {{0, 1, 3}},
    {{0, 2, 3}},
    {{1, 2, 3}}
}};

// Count how many times each triangular face appears
for (int i = 0; i < my_mesh->GetNumElements(); ++i) {

    auto elem =
    std::dynamic_pointer_cast<ChElementLDPM>(
    my_mesh->GetElement(i));

    for (auto face : tet_faces) {

        auto n1 = elem->GetTetrahedronNode(face[0]);
        auto n2 = elem->GetTetrahedronNode(face[1]);
        auto n3 = elem->GetTetrahedronNode(face[2]);

        std::array<int, 3> key = {
            n1->GetIndex(),
            n2->GetIndex(),
            n3->GetIndex()
        };

        std::sort(key.begin(), key.end());

        face_count[key]++;

        if (face_count[key] == 1) {
            face_map[key] = {n1, n2, n3};
        }
    }
}

// Keep only faces that appear once and lie on the lateral surface
for (auto item : face_count) {

    if (item.second != 1)
        continue;

    BoundaryFace bf = face_map[item.first];

    ChVector3d p1 = bf.n1->GetPos();
    ChVector3d p2 = bf.n2->GetPos();
    ChVector3d p3 = bf.n3->GetPos();

    ChVector3d pc = (p1 + p2 + p3) / 3.0;

    double dx = pc.x() - center_x;
    double dy = pc.y() - center_y;
    double rc = sqrt(dx * dx + dy * dy);

    if (fabs(rc - R) < tol_r &&
    pc.z() > zmin + tol_z &&
    pc.z() < zmax - tol_z) {

        lateral_faces.push_back(bf);
    }
}

std::cout << "Lateral boundary faces number: "
          << lateral_faces.size()
          << std::endl;

// Store initial volume for each LDPM tetrahedron
std::vector<double> tet_initial_volume;

for (int i = 0; i < my_mesh->GetNumElements(); ++i) {

    auto elem =
    std::dynamic_pointer_cast<ChElementLDPM>(
    my_mesh->GetElement(i));

    auto n1 = elem->GetTetrahedronNode(0);
    auto n2 = elem->GetTetrahedronNode(1);
    auto n3 = elem->GetTetrahedronNode(2);
    auto n4 = elem->GetTetrahedronNode(3);

    double V0 = TetraVolume(
        n1->GetX0().GetPos(),
        n2->GetX0().GetPos(),
        n3->GetX0().GetPos(),
        n4->GetX0().GetPos()
    );

    tet_initial_volume.push_back(V0);
}


    double u = 0;
    double F = 0;
    double Wext = 0;

    //std::vector<int> N_iter;



 //while (vis->Run() & sys.GetChTime() <= 0.2) {
	while (sys.GetChTime() <= T0) {
			

// Linear pressure ramp
    P_applied = Pmax * sys.GetChTime() / T0;

    if (P_applied > Pmax)
    P_applied = Pmax;

// Total vertical force = pressure * specimen area
         F_top = P_applied * A_specimen;

// Apply compressive force downward on top plate
    top_plate_force->SetMforce(F_top);    




		//vis->BeginScene();
		//vis->Render();  
		//vis->EndScene(); 
	// ==========================================

    
    // ==========================================
// ==========================================
// Lateral pressure using true boundary faces
// ==========================================

std::map<std::shared_ptr<ChNodeFEAxyzrot>, ChVector3d>
lateral_force_map;

// Initialize force map for lateral nodes
for (auto node : sur_nodes) {
    lateral_force_map[node] = ChVector3d(0, 0, 0);
}

// Loop only over pre-detected lateral boundary faces
for (auto bf : lateral_faces) {

    ChVector3d p1 = bf.n1->GetPos();
    ChVector3d p2 = bf.n2->GetPos();
    ChVector3d p3 = bf.n3->GetPos();

    double A_face = TriangleArea(p1, p2, p3);
    double F_face = P_applied * A_face;
    double F_node = F_face / 3.0;

    std::vector<std::shared_ptr<ChNodeFEAxyzrot>>
        face_nodes = { bf.n1, bf.n2, bf.n3 };

    for (auto node : face_nodes) {

        ChVector3d p = node->GetX0().GetPos();

        double dx = p.x() - center_x;
        double dy = p.y() - center_y;
        double r = sqrt(dx * dx + dy * dy);

        if (r > 1e-12) {
            double rx = dx / r;
            double ry = dy / r;

            lateral_force_map[node] +=
                ChVector3d(-F_node * rx, -F_node * ry, 0.0);
        }
    }
}

// Apply accumulated nodal forces
for (auto item : lateral_force_map) {
    if (item.first != RP1 && item.first != RP2) {
        item.first->SetForce(item.second);
    }
}

		
		sys.DoStepDynamics(timestep);  
		
        stepnum++;

/*
        double du = motor1->GetMotorPos() - u;
        Wext = Wext + abs(du * (motor1->GetMotorForce() + F) / 2);
        u = motor1->GetMotorPos();
        F = motor1->GetMotorForce();
*/

        F = F_top;
        int n_iter = mystepper->GetNumStepIterations();
        std::cout << "n_iter= " << n_iter << std::endl;
        num_iter_file << "\tt=\t" << sys.GetChTime()<< "\tstep=\t" << stepnum << "\tn_iter=\t" << n_iter << "\t\n";
        num_iter_file.flush();
        //N_iter.push_back(n_iter);

		
		if(stepnum%25==0) {

            double Wint = 0;
            for (int i = 0; i < my_mesh->GetNumElements(); ++i) {
                auto elem = std::dynamic_pointer_cast<ChElementLDPM>(my_mesh->GetElement(i));

                for (auto facet : elem->GetSection()) {
                    auto statev = facet->Get_StateVar();
                    Wint = Wint + statev(10);
                }
            }

            double Ek = 0;
            for (int i = 0; i < my_mesh->GetNumElements(); ++i) {
                auto elem = std::dynamic_pointer_cast<ChElementLDPM>(my_mesh->GetElement(i));

                ChVectorN<double, 24> V;
                V.setZero();
                for (int j = 0; j < 4; j++) {
                    auto node = elem->GetTetrahedronNode(j);
                    V((j + 1) * 6 - 6) = node->GetPosDt().x();
                    V((j + 1) * 6 - 5) = node->GetPosDt().y();
                    V((j + 1) * 6 - 4) = node->GetPosDt().z();
                    V((j + 1) * 6 - 3) = node->GetAngVelLocal().x();
                    V((j + 1) * 6 - 2) = node->GetAngVelLocal().y();
                    V((j + 1) * 6 - 1) = node->GetAngVelLocal().z();
                }


               
                //std::cout << " V" << V << std::endl;
                ChMatrixNM<double, 24, 24> M;
                M.setZero();
                elem->ComputeMmatrixGlobal(M);

                double Ekp = 0.5 * V.transpose() * M * V;
                Ek = Ek + Ekp;
            }

double min_vol_strain = 1e20;
double max_vol_strain = -1e20;


        // Compute volume-weighted average volumetric strain
double total_volume = 0.0;
double weighted_ev_sum = 0.0;

for (int i = 0; i < my_mesh->GetNumElements(); ++i) {
    auto elem = std::dynamic_pointer_cast<ChElementLDPM>(my_mesh->GetElement(i));

    auto n1 = elem->GetTetrahedronNode(0);
    auto n2 = elem->GetTetrahedronNode(1);
    auto n3 = elem->GetTetrahedronNode(2);
    auto n4 = elem->GetTetrahedronNode(3);

    double V0 = tet_initial_volume[i];

    double V = TetraVolume(
        n1->GetPos(),
        n2->GetPos(),
        n3->GetPos(),
        n4->GetPos()
    );

    double ev = ((V - V0) / V0)/3;
        if (ev < min_vol_strain)
        min_vol_strain = ev;

        if (ev > max_vol_strain)
        max_vol_strain = ev;


       weighted_ev_sum += ev * V0;
       total_volume += V0;
}

    double avg_vol_strain = weighted_ev_sum / total_volume;    

	    	std::string mesh_filename=out_dir+"deneme"+std::to_string(stepnum)+".vtk";
	    	std::string vtk_filename=out_dir+"Vtkdeneme"+std::to_string(stepnum)+".vtk";
	    	WriteMesh(my_mesh, mesh_filename);
	    	WriteFrame(my_mesh, mesh_filename, vtk_filename);

            std::string mesh_filename1 = out_dir + "crack" + std::to_string(stepnum) + ".vtk";
            std::string vtk_filename1 = out_dir + "Vtkcrack" + std::to_string(stepnum) + ".vtk";
            WriteMesh1(my_mesh, mesh_filename1);
            WriteFrame1(my_mesh, mesh_filename1, vtk_filename1);
		

        
        double top_disp = top_plate->GetPos().z() - zmax;
        double support_Rz = constr_bot->GetReaction1().force.z();
        double sigma_bottom = support_Rz / A_specimen;

std::cout << " t=\t" << sys.GetChTime()
          << "\ttop_disp_z=\t" << top_disp
          << "\tP_applied=\t" << P_applied
          << "\tAppliedForce=\t" << F_top
          << "\tSupport_Rz=\t" << support_Rz
          << "\tsigma_bottom=\t" << sigma_bottom
          << "\tavg_vol_strain=\t" << avg_vol_strain
          << "\tinternal_work\t" << Wint
          << "\texternal_work\t" << Wext
          << "\tkinetic energy\t" << Ek 
          << "\tmin_ev=\t" << min_vol_strain
          << "\tmax_ev=\t" << max_vol_strain
          << "\t\n";

histfile << " t=\t" << sys.GetChTime()
         << "\ttop_disp_z=\t" << top_disp
         << "\tP_applied=\t" << P_applied
         << "\tAppliedForce=\t" << F_top
         << "\tSupport_Rz=\t" << support_Rz
         << "\tsigma_bottom=\t" << sigma_bottom
         << "\tavg_vol_strain=\t" << avg_vol_strain
         << "\tinternal_work\t" << Wint
         << "\texternal_work\t" << Wext
         << "\tkinetic energy\t" << Ek 
         << "\tmin_ev=\t" << min_vol_strain
         << "\tmax_ev=\t" << max_vol_strain
         << "\t\n";
		histfile.flush();

		}

	    }

	    histfile.close();
        num_iter_file.close();
		
	   };
   	
    return 0;
}
        
