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
// Authors: Jibril B. Coulibaly
// =============================================================================

#include <cmath>
#include <initializer_list>
#include <map>
#include <memory>
#include "chrono/core/ChFrame.h"
#include "chrono/core/ChMatrix33.h"
#include "chrono/core/ChQuaternion.h"
#include "chrono/core/ChTypes.h"
#include "chrono/core/ChVector3.h"
#include "chrono/fea/ChElementBase.h"
#include "chrono/fea/ChMesh.h"
#include "chrono/fea/ChNodeFEAxyzrot.h"
#include "chrono/physics/ChSystemSMC.h"
#include "chrono_ldpm/ChSectionLDPM.h"
#include "gtest/gtest.h"

#include "chrono_ldpm/ChElementLDPM.h"

#include <chrono>
#include <utility>
#include <vector>

using Moment = std::chrono::high_resolution_clock::time_point;
using FloatSecs = std::chrono::duration<double>;

inline Moment now()
{
    return std::chrono::high_resolution_clock::now();
}

using namespace chrono;
using namespace fea;
using namespace ldpm;

// TODO: Write a Fixture to avoid repeating the code for setting up the test
// TODO: Branch and duplicate test with / without largeDeflection (when implemented) and coupled forces


TEST(LDPMFacetTest, internal_forces) {
    //
    ChSystemSMC sys;
    auto my_mesh = chrono_types::make_shared<ChMesh>();
    my_mesh->SetAutomaticGravity(false);
    sys.Add(my_mesh);

    auto my_LDPM_tet = chrono_types::make_shared<ChElementLDPM>();
    my_mesh->AddElement(my_LDPM_tet);

    auto node1 = chrono_types::make_shared<fea::ChNodeFEAxyzrot>(ChFrame(ChVector3<>(0,0,0), ChQuaternion<>(1, 0, 0, 0)));
    auto node2 = chrono_types::make_shared<fea::ChNodeFEAxyzrot>(ChFrame(ChVector3<>(1,0,0), ChQuaternion<>(1, 0, 0, 0)));
    auto node3 = chrono_types::make_shared<fea::ChNodeFEAxyzrot>(ChFrame(ChVector3<>(0,1,0), ChQuaternion<>(1, 0, 0, 0)));
    auto node4 = chrono_types::make_shared<fea::ChNodeFEAxyzrot>(ChFrame(ChVector3<>(0,0,1), ChQuaternion<>(1, 0, 0, 0)));
    std::vector<std::shared_ptr<ChNodeFEAxyzrot>> nodes {node1, node2, node3, node4};
    for (auto node : nodes) {
        my_mesh->AddNode(node);
    }
    my_LDPM_tet->SetNodes(node1, node2, node3, node4);

    // Assume all grain radii equal to zero to simplify construction
    ChVector3<> edge12 = 0.5 * (node1->GetPos() + node2->GetPos());
    ChVector3<> edge13 = 0.5 * (node1->GetPos() + node3->GetPos());
    ChVector3<> edge14 = 0.5 * (node1->GetPos() + node4->GetPos());
    ChVector3<> edge23 = 0.5 * (node2->GetPos() + node3->GetPos());
    ChVector3<> edge24 = 0.5 * (node2->GetPos() + node4->GetPos());
    ChVector3<> edge34 = 0.5 * (node3->GetPos() + node4->GetPos());
    ChVector3<> face1 = (0.5 * (edge34 + node2->GetPos()) + 0.5 * (edge23 + node4->GetPos()) + 0.5 * (edge24 + node3->GetPos())) / 3;
    ChVector3<> face2 = (0.5 * (edge34 + node1->GetPos()) + 0.5 * (edge14 + node3->GetPos()) + 0.5 * (edge13 + node4->GetPos())) / 3;
    ChVector3<> face3 = (0.5 * (edge24 + node1->GetPos()) + 0.5 * (edge14 + node2->GetPos()) + 0.5 * (edge12 + node4->GetPos())) / 3;
    ChVector3<> face4 = (0.5 * (edge12 + node3->GetPos()) + 0.5 * (edge23 + node1->GetPos()) + 0.5 * (edge13 + node2->GetPos())) / 3;
    ChVector3<> tet = (0.5 * (face1 + node1->GetPos()) + 0.5 * (face2 + node2->GetPos()) + 0.5 * (face3 + node3->GetPos()) + 0.5 * (face4 + node4->GetPos())) / 4;

    auto my_mat = chrono_types::make_shared<ChMaterialVECT>();
    double E0 = 43748;
    double alpha = 0.25;
    // Assume elasticity (stress calculation is not the responsibility of this test)
    // by selecting very larg stress limits
    double sigmat = 1e10;
    double sigmas = 1e10;
    double sigmac0 = 1e10;
    double sigmaN0 = 1e10;
    my_mat->Set_density(2.370E-9);
    my_mat->Set_E0(E0);
    my_mat->Set_alpha(alpha);
    my_mat->Set_sigmat(sigmat);
    my_mat->Set_sigmas(sigmas);
    my_mat->Set_nt(0.2);
    my_mat->Set_lt(250);
    my_mat->Set_Ed(43748);
    my_mat->Set_sigmac0(sigmac0);
    my_mat->Set_Hc0(17499);
    my_mat->Set_Hc1(4374.8);
    my_mat->Set_beta(0);
    my_mat->Set_kc0(2);
    my_mat->Set_kc1(1);
    my_mat->Set_kc2(5);
    my_mat->Set_kc3(0.1);
    my_mat->Set_mu0(0.2);
    my_mat->Set_muinf(0);
    my_mat->Set_sigmaN0(sigmaN0);

    std::vector<std::vector<ChVector3<>>> vertices_map{{node2->GetPos()-node1->GetPos(), edge12, face3, face4},
                                                       {node3->GetPos()-node1->GetPos(), edge13, face2, face4},
                                                       {node4->GetPos()-node1->GetPos(), edge14, face2, face3},
                                                       {node3->GetPos()-node2->GetPos(), edge23, face1, face4},
                                                       {node4->GetPos()-node2->GetPos(), edge24, face1, face3},
                                                       {node4->GetPos()-node3->GetPos(), edge34, face1, face2}};
    for (int i=0 ; i < 6 ; i++) {
        auto node_to_node = vertices_map[i][0];
        double length = node_to_node.Length();
        auto edge = vertices_map[i][1];
        for (int f : {2,3}) {
            auto face = vertices_map[i][f];
            ChVector3<>center = (tet + edge + face) / 3;
            ChVector3<> normal = (edge - tet).Cross(face - tet);;
            double area = 0.5 * normal.Dot(node_to_node) / length; // Projected area

            // The node-to-node direction is used as the facet normal "n_k" (as opposed to `normal` , which is "n_k0" in the paper)
            ChMatrix33<> facetFrame;
            facetFrame.SetFromAxisX(node_to_node, face - tet); // choice of Y-axis does not matter
            // LDPM expects the facet frame matrix to be stored as the transpose of the rotation matrix:
            // - n on the first row
            // - m on the second row
            // - l on the second row
            // TODO JBC: I think this should eventually be changed as it might be confusing
            // that this is different from all the other "frames" defined as ChMatrix33 in Chrono
            auto section = chrono_types::make_shared<ChSectionLDPM>(my_mat, area, center, facetFrame.transpose());
            section->Set_Length(length);
            my_LDPM_tet->AddFacetI(section);

            // Add vertices for initial mass calculation
            // TODO JBC: not sure I am doing this correctly (order in particular) but the whole design is pretty stupid if I am being honest
            auto nonsense_unnecessary_xyzrot_1 = chrono_types::make_shared<ChNodeFEAxyzrot>(ChFrame<>(tet, QUNIT));
            auto nonsense_unnecessary_xyzrot_2 = chrono_types::make_shared<ChNodeFEAxyzrot>(ChFrame<>(edge, QUNIT));
            auto nonsense_unnecessary_xyzrot_3 = chrono_types::make_shared<ChNodeFEAxyzrot>(ChFrame<>(face, QUNIT));
            my_LDPM_tet->AddVertNodeVec(std::vector<std::shared_ptr<fea::ChNodeFEAxyzrot>> {nonsense_unnecessary_xyzrot_1,nonsense_unnecessary_xyzrot_2,nonsense_unnecessary_xyzrot_3});
        }
    }


    // Small deflection can be used since the void ChElementLDPM::ComputeStrain function only really performs the
    // computation of the strain according to https://doi.org/10.1016/j.cemconcomp.2011.02.011
    my_LDPM_tet->SetLargeDeflection(false);

    // Performing some kind of step / analysis is required for chrono to run the appropriate setup
    // Functions setupInitial() etc, are all private functions and inaccessible from here, otherwise we use them to make a minimal setup
    sys.DoStepDynamics(0);

    // Displace and rotate node 4 to cause imposed local strain
    // Strain only on facets 1-4, 2-4, 3-4
    ChVector3d dispB(0.0, 0.0, 0.01);
    node4->SetPos(node4->GetPos() + dispB);
    ChQuaterniond qB_ini = node4->GetRot();
    double drotB_angle = 1e-2;
    ChVector3d drotB_axis(0.46, 1.5, -0.5);
    ChVector3d drotB = drotB_angle * drotB_axis.GetNormalized();
    ChQuaterniond dqB;
    dqB.SetFromRotVec(drotB);
    ChQuaterniond qB = dqB * qB_ini;
    node4->SetRot(qB);
    // TODO JBC: Current rotational DOFs increment calculation inside ComputeInternalForces() is additive in terms of rotation vector
    // Multiplicative composition of quaternion/rotation matrices, which is exact, does not give the same result:
    // If the increment of rotation is given by the quaternion dq, then q2 = dq * q1. The increment rotation vector is RotVec(dq)
    // However, RotVec(q2) - RotVec(q1) != RotVec(dq). For small rotations dq, this is a fair enough approximation.
    // But this is wrong enough for testing! The error is in the order of magnitude of the angle, i.e., 1% error for rotB_angle = 1e-2 rad.
    // To pass the test, and actually exercise the code the way it is currently written (even though that way is wrong !!!) we are using additive rotation vector decomposition
    // TODO JBC: The commented line below is the mathematically correct multiplicative decomposition
    // ChVector3d rot_incr = drotB;
    // TODO JBC: The line below is for the small rotation assumption used in the code.
    // Once/if we debug and refactor this aspect of the code, this test should fail and we should use the commented line above instead
    ChVector3d rot_incr = qB.GetRotVec() - qB_ini.GetRotVec();

    // Chrono LDPM calculation of the internal forces
    ChVectorDynamic<> Fi(24);
    my_LDPM_tet->ComputeInternalForces(Fi);

    // Analytical calculation of the internal forces
    // Equation (12) (13) of https://doi.org/10.1016/j.cemconcomp.2011.02.011
    // Assume elasticity (stress calculation is not the responsibility of this test)
    // To be performed and aggregated on each of the 12 facets
    // Facet node numbering ChElementLDPM::facetNodeNums << 0, 1, 0, 1, 0, 2, 0, 2, 0, 3, 0, 3, 1, 2, 1, 2, 1, 3, 1, 3, 2, 3, 2, 3;
    std::vector<int> nodeIind = {0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2};
    std::vector<int> nodeJind = {1, 1, 2, 2, 3, 3, 2, 2, 3, 3, 3, 3};
    ChVectorN<double, 24> Fi_analytical;
    Fi_analytical.setZero();

    // Reset position and orientation for analytical calculation
    node4->SetPos(node4->GetPos() - dispB);
    node4->SetRot(QUNIT);

    // Strain only on facets 1-4, 2-4, 3-4
    for (int i : {4, 5, 8, 9, 10, 11}) {
        auto section = my_LDPM_tet->GetSection()[i];
        auto center = section->Get_center();
        auto facetFrame = section->Get_facetFrame();
        auto area = section->Get_area();
        auto length = section->Get_Length();
        
        ChVectorDynamic<> Fi_facet(12);
        ChVector3<> posNodeA = nodes[nodeIind[i]]->GetPos();
        ChVector3<> posNodeB = nodes[nodeJind[i]]->GetPos();
        

        ChVector3<> imposed_local_strain = facetFrame * (dispB + rot_incr.Cross(center - posNodeB)) / length; // This is not very general, but sufficient if we only move node 4
        ChVector3<> stress_local = imposed_local_strain * E0 * ChVector3d(1.0, alpha, alpha);
        ChVector3<> xc_xi = center - posNodeA;
        ChVector3<> xc_xj = center - posNodeB;
        ChMatrixNM<double,3,6> Ai, Aj;
        Ai << 1.0, 0.0, 0.0, 0.0      ,  xc_xi[2], -xc_xi[1],
              0.0, 1.0, 0.0, -xc_xi[2],  0.0     ,  xc_xi[0],
              0.0, 0.0, 1.0,  xc_xi[1], -xc_xi[0],  0.0     ;
        Aj << 1.0, 0.0, 0.0, 0.0      ,  xc_xj[2], -xc_xj[1],
              0.0, 1.0, 0.0, -xc_xj[2],  0.0     ,  xc_xj[0],
              0.0, 0.0, 1.0,  xc_xj[1], -xc_xj[0],  0.0     ;
        ChVectorN<double, 6> Bn_i_tr, Bm_i_tr, Bl_i_tr;
        ChVectorN<double, 6> Bn_j_tr, Bm_j_tr, Bl_j_tr;
        ChMatrix33<> facetFrame_tr = facetFrame.transpose();
        Bn_i_tr = Ai.transpose() * facetFrame_tr.GetAxisX().eigen() / length;
        Bm_i_tr = Ai.transpose() * facetFrame_tr.GetAxisY().eigen() / length;
        Bl_i_tr = Ai.transpose() * facetFrame_tr.GetAxisZ().eigen() / length;
        Bn_j_tr = Aj.transpose() * facetFrame_tr.GetAxisX().eigen() / length;
        Bm_j_tr = Aj.transpose() * facetFrame_tr.GetAxisY().eigen() / length;
        Bl_j_tr = Aj.transpose() * facetFrame_tr.GetAxisZ().eigen() / length;
        Fi_facet.segment(0,6) = - length * area * (stress_local[0] * Bn_i_tr + stress_local[1] * Bm_i_tr + stress_local[2] * Bl_i_tr);
        Fi_facet.segment(6,6) =   length * area * (stress_local[0] * Bn_j_tr + stress_local[1] * Bm_j_tr + stress_local[2] * Bl_j_tr);

        Fi_analytical.segment(nodeIind[i]*6, 6) += -Fi_facet.segment(0,6);
        Fi_analytical.segment(nodeJind[i]*6, 6) += -Fi_facet.segment(6,6);
    }

    double tol = 1e-10;
    for (int i = 0; i < 24; i++) {
        ASSERT_NEAR(Fi(i), Fi_analytical(i), tol);
    }




}

TEST(LDPM_misc_tests, multiple_volume_calc) {
    // ChSystemSMC sys;
    
    // auto my_mesh = chrono_types::make_shared<ChMesh>();
    // my_mesh->SetAutomaticGravity(false);
    // sys.Add(my_mesh);

    // ChElementLDPM my_LDPM_tet;
    // auto nodeA = chrono_types::make_shared<fea::ChNodeFEAxyzrot>(ChFrame(ChVector3<>(0,0,0), ChQuaternion<>(1, 0, 0, 0)));
    // auto nodeB = chrono_types::make_shared<fea::ChNodeFEAxyzrot>(ChFrame(ChVector3<>(1,0,0), ChQuaternion<>(1, 0, 0, 0)));
    // auto nodeC = chrono_types::make_shared<fea::ChNodeFEAxyzrot>(ChFrame(ChVector3<>(0,1,0), ChQuaternion<>(1, 0, 0, 0)));
    // auto nodeD = chrono_types::make_shared<fea::ChNodeFEAxyzrot>(ChFrame(ChVector3<>(0,0,1), ChQuaternion<>(1, 0, 0, 0)));
    
    // my_LDPM_tet.SetNodes(nodeA, nodeB, nodeC, nodeD);


    // std::cout<<"ComputeVolume()" <<my_LDPM_tet.ComputeVolume()<<std::endl;
    

    // for (int i = 0 ; i < 100 ; i++) {
    //     double rx = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
    //     double ry = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
    //     double rz = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
    //     nodeD->SetPos(ChVector3<>(rx, ry, rz));
    //     std::cout<<my_LDPM_tet.ComputeVolume() << " ";
    // }

    // std::cout<< " PROFILING"<<std::endl;




    // double dummy = 0;
    // auto startTime = now();
    // for (int i = 0 ; i < 100000 ; i++) {
    //     nodeB->SetPos(nodeB->GetPos() + ChVector3<>(1e-4, 0, 0));
    //     dummy += my_LDPM_tet.ComputeVolume();
    // }
    // auto elapsed = now() - startTime;
    // auto seconds = std::chrono::duration_cast<FloatSecs>(elapsed);
    // std::cout << "ComputeVolume()" << seconds.count() << " dummy = " << dummy << std::endl;


}




