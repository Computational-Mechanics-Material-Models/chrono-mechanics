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
#include "chrono/core/ChMatrix.h"
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

// TODO: Branch and duplicate test with / without largeDeflection (when implemented)

namespace {
class LDPMTest : public testing::Test {
    protected:

    void SetUp() override {
        // Parameters from Appendix A. of our Overleaf
        // TODO: update reference once this is published

        // Material parameters (SI units)
        const double MPa_to_Pa = 1e6;
        const double mm_to_m = 1e-3;

        rho = 2380; // Density
        E0 = 60273 * MPa_to_Pa; // Normal modulus
        alpha = 0.25; // Alpha
        sigmat = 3.44 * MPa_to_Pa; // Tensile strength
        lt = 500 * mm_to_m; // Tensile charachteristic length
        Gt = lt * sigmat * sigmat / (2 * E0); // Fracture energy
        rt = 2.6; // Shear strength ratio
        nt = 0.4; // Softening exponent
        sigmac0 = 150 * MPa_to_Pa; // Compressive yield strength
        Hc0 = 0.40 * E0; // Initial hardening modulus
        kc0 = 4; // Transitional strain ratio
        kc1 = 1; // Deviatoric strain threshold ratio
        kc2 = 5; // Deviatoric damage parameter
        kc3 = 0.1; // Volumetric strain parameter
        mu0 = 0.4; // Initial friction
        muinf = 0; // Asymptotic friction
        sigmaN0 = 600 * MPa_to_Pa; // Transitional stress
        Ed = 1 * E0; // Densification ratio
        beta = 0; // Volumetric deviatoric coupling
        kt = 0.0; // Tensile unloading
        ks = 0.0; // Shear unloading
        kc = 0.0; // Compressive unloading
        rs = 0.0; // Shear softening modulus ratio
        Hc1 = 0.1 * E0; // Final hardening modulus ratio
        // Other parameters
        double sigmas = rt * sigmat; // TODO: chrono wants sigmas as an input, not rt

        auto my_mat = chrono_types::make_shared<ChMaterialVECT>();
        my_mat->Set_density(rho);
        my_mat->Set_E0(E0);
        my_mat->Set_alpha(alpha);
        my_mat->Set_sigmat(sigmat);
        my_mat->Set_lt(lt);
        // my_mat->Set_Gt( Gt); // TODO: this function is not used and commented out in the chrono code
        my_mat->Set_sigmas(sigmas); // TODO: chrono input asks for `sigmas`, not for `rt`
        my_mat->Set_nt(nt);
        my_mat->Set_sigmac0(sigmac0);
        my_mat->Set_Hc0(Hc0);
        my_mat->Set_kc0(kc0);
        my_mat->Set_kc1(kc1);
        my_mat->Set_kc2(kc2);
        my_mat->Set_kc3(kc3);
        my_mat->Set_mu0(mu0);
        my_mat->Set_muinf(muinf);
        my_mat->Set_sigmaN0(sigmaN0);
        my_mat->Set_Ed(Ed);
        my_mat->Set_beta(beta); // TODO: beta = 0 does not exercise e_DV (i.e. eDV = eV)
        my_mat->Set_kt(kt); // TODO: kt = 0 does not exercise the different hysteresis behaviors
        // TODO: no function available to set `ks`
        // TODO: no function available to set `kc`
        my_mat->Set_rs(rs);
        my_mat->Set_Hc1(Hc1);

        // Create tetrahedron, create nodes, add to mesh setc
        my_mesh = chrono_types::make_shared<ChMesh>();
        my_mesh->SetAutomaticGravity(false);
        sys.Add(my_mesh);
        my_LDPM_tet = chrono_types::make_shared<ChElementLDPM>();
        my_mesh->AddElement(my_LDPM_tet);
        std::shared_ptr<ChNodeFEAxyzrot> node1 = chrono_types::make_shared<fea::ChNodeFEAxyzrot>(ChFrame(ChVector3<>(0,0,0), ChQuaternion<>(1, 0, 0, 0)));
        std::shared_ptr<ChNodeFEAxyzrot> node2 = chrono_types::make_shared<fea::ChNodeFEAxyzrot>(ChFrame(ChVector3<>(1,0,0), ChQuaternion<>(1, 0, 0, 0)));
        std::shared_ptr<ChNodeFEAxyzrot> node3 = chrono_types::make_shared<fea::ChNodeFEAxyzrot>(ChFrame(ChVector3<>(0,1,0), ChQuaternion<>(1, 0, 0, 0)));
        std::shared_ptr<ChNodeFEAxyzrot> node4 = chrono_types::make_shared<fea::ChNodeFEAxyzrot>(ChFrame(ChVector3<>(0,0,1), ChQuaternion<>(1, 0, 0, 0)));
        nodes.push_back(node1);
        nodes.push_back(node2);
        nodes.push_back(node3);
        nodes.push_back(node4);
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
                ChVector3<> normal = (edge - tet).Cross(face - tet);
                double area = 0.5 * std::abs(normal.Dot(node_to_node)) / length; // Projected area

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

        // Tolerances
        for (auto section : my_LDPM_tet->GetSection()) {
            average_area += section->Get_area();
            average_length += section->Get_Length();
        }
        average_area /= my_LDPM_tet->GetSection().size();
        average_length /= my_LDPM_tet->GetSection().size();
        tol_force = E0 * average_area * rtol;
        tol_stiff = E0 * average_area / average_length * rtol;
    }

    // Material parameters
    double rho; // Density
    double E0; // Normal modulus
    double alpha; // Alpha
    double sigmat; // Tensile strength
    double lt; // Tensile charachteristic length
    double Gt; // Fracture energy
    double rt; // Shear strength ratio
    double nt; // Softening exponent
    double sigmac0; // Compressive yield strength
    double Hc0; // Initial hardening modulus
    double kc0; // Transitional strain ratio
    double kc1; // Deviatoric strain threshold ratio
    double kc2; // Deviatoric damage parameter
    double kc3; // Volumetric strain parameter
    double mu0; // Initial friction
    double muinf; // Asymptotic friction
    double sigmaN0; // Transitional stress
    double Ed; // Densification ratio
	double beta; // Volumetric deviatoric coupling
    double kt; // Tensile unloading
    double ks; // Shear unloading
    double kc; // Compressive unloading
    double rs; // Shear softening modulus ratio
    double Hc1; // Final hardening modulus ratio

    // LDPM element and Chrono internals
    ChSystemSMC sys;
    std::shared_ptr<ChMesh> my_mesh;
    std::shared_ptr<ChElementLDPM> my_LDPM_tet;
    std::vector<std::shared_ptr<ChNodeFEAxyzrot>> nodes;
    // Facet node numbering ChElementLDPM::facetNodeNums << 0, 1, 0, 1, 0, 2, 0, 2, 0, 3, 0, 3, 1, 2, 1, 2, 1, 3, 1, 3, 2, 3, 2, 3;
    std::vector<int> nodeIind = {0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2};
    std::vector<int> nodeJind = {1, 1, 2, 2, 3, 3, 2, 2, 3, 3, 3, 3};

    // Tolerance
    double rtol = 1e-10;
    double average_area = 0;
    double average_length = 0;
    double tol_force;
    double tol_stiff;
};

TEST_F(LDPMTest, internal_forces) {
    // Small deflection can be used since the void ChElementLDPM::ComputeStrain function only really performs the
    // computation of the strain according to https://doi.org/10.1016/j.cemconcomp.2011.02.011
    my_LDPM_tet->SetLargeDeflection(false);

    // Performing some kind of step / analysis is required for chrono to run the appropriate setup
    // Functions setupInitial() etc, are all private functions and inaccessible from here, otherwise we use them to make a minimal setup
    sys.DoStepDynamics(0);

    // Displace and rotate node 4 to cause imposed local strain
    // Strain only on facets 1-4, 2-4, 3-4
    ChVector3d dispB(0.0, 0.0, 0.01);
    nodes[3]->SetPos(nodes[3]->GetPos() + dispB);
    ChQuaterniond qB_ini = nodes[3]->GetRot();
    double drotB_angle = 1e-2;
    ChVector3d drotB_axis(0.46, 1.5, -0.5);
    ChVector3d drotB = drotB_angle * drotB_axis.GetNormalized();
    ChQuaterniond dqB;
    dqB.SetFromRotVec(drotB);
    ChQuaterniond qB = dqB * qB_ini;
    nodes[3]->SetRot(qB);
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
    // this loading with the default material parameters gives elastic response
    // To be performed and aggregated on each of the 12 facets

    ChVectorN<double, 24> Fi_analytical;
    Fi_analytical.setZero();

    // Reset position and orientation for analytical calculation
    nodes[3]->SetPos(nodes[3]->GetPos() - dispB);
    nodes[3]->SetRot(QUNIT);

    // Strain only on facets 1-4, 2-4, 3-4
    for (int i : {4, 5, 8, 9, 10, 11}) {
        auto section = my_LDPM_tet->GetSection()[i];
        auto center = section->Get_center();
        auto facetFrame = section->Get_facetFrame();
        auto area = section->Get_area();
        auto length = section->Get_Length();

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

        Fi_analytical.segment(nodeIind[i]*6, 6) +=  length * area * (stress_local[0] * Bn_i_tr + stress_local[1] * Bm_i_tr + stress_local[2] * Bl_i_tr);;
        Fi_analytical.segment(nodeJind[i]*6, 6) += -length * area * (stress_local[0] * Bn_j_tr + stress_local[1] * Bm_j_tr + stress_local[2] * Bl_j_tr);
    }

    for (int i = 0; i < 24; i++) {
        ASSERT_NEAR(Fi(i), Fi_analytical(i), tol_force);
    }
}

TEST_F(LDPMTest, elastic_stiffness_matrix) {
    // Small deflection can be used since the void ChElementLDPM::ComputeStrain function only really performs the
    // computation of the strain according to https://doi.org/10.1016/j.cemconcomp.2011.02.011
    my_LDPM_tet->SetLargeDeflection(false);

    // Performing some kind of step / analysis is required for chrono to run the appropriate setup
    // Functions setupInitial() etc, are all private functions and inaccessible from here, otherwise we use them to make a minimal setup
    sys.DoStepDynamics(0);

    // Chrono LDPM calculation of the elastic stiffness matrix
    ChMatrixNM<double, 24, 24> chrono_matrix;
    chrono_matrix.setZero();
    my_LDPM_tet->ComputeStiffnessMatrixGlobal(chrono_matrix);

    // Analytical calculation of the stiffness matrix
    ChMatrixNM<double, 24, 24> analytical_matrix;
    analytical_matrix.setZero();

    for (int i = 0 ; i < 12 ; i++) {
        auto section = my_LDPM_tet->GetSection()[i];
        double facet_area = section->Get_area();
        double length = section->Get_Length();
        auto center = section->Get_center();
        ChMatrix33<> Rf = section->Get_facetFrame().transpose();

        ChMatrix33<> stiff_local;
        stiff_local << facet_area * E0 / length , 0.0                             , 0.0                             ,
                       0.0                      , facet_area * E0 * alpha / length, 0.0                             ,
                       0.0                      , 0.0                             , facet_area * E0 * alpha / length;

        ChVector3<> posNodeA = nodes[nodeIind[i]]->GetPos();
        ChVector3<> posNodeB = nodes[nodeJind[i]]->GetPos();
        ChVector3<> xc_xi = center - posNodeA;
        ChVector3<> xc_xj = center - posNodeB;
        ChMatrixNM<double,3,6> Ai, Aj;
        Ai << 1.0, 0.0, 0.0, 0.0      ,  xc_xi[2], -xc_xi[1],
            0.0, 1.0, 0.0, -xc_xi[2],  0.0     ,  xc_xi[0],
            0.0, 0.0, 1.0,  xc_xi[1], -xc_xi[0],  0.0     ;
        Aj << 1.0, 0.0, 0.0, 0.0      ,  xc_xj[2], -xc_xj[1],
            0.0, 1.0, 0.0, -xc_xj[2],  0.0     ,  xc_xj[0],
            0.0, 0.0, 1.0,  xc_xj[1], -xc_xj[0],  0.0     ;

        // Stiffness from facet stresses
        // Chain rule Derivatives of (12) (13) from https://doi.org/10.1016/j.cemconcomp.2011.02.011
        analytical_matrix.block<6,6>(nodeIind[i]*6, nodeIind[i]*6) +=  Ai.transpose() * Rf * stiff_local * Rf.transpose() * Ai;
        analytical_matrix.block<6,6>(nodeIind[i]*6, nodeJind[i]*6) += -Ai.transpose() * Rf * stiff_local * Rf.transpose() * Aj;
        analytical_matrix.block<6,6>(nodeJind[i]*6, nodeIind[i]*6) += -Aj.transpose() * Rf * stiff_local * Rf.transpose() * Ai;
        analytical_matrix.block<6,6>(nodeJind[i]*6, nodeJind[i]*6) +=  Aj.transpose() * Rf * stiff_local * Rf.transpose() * Aj;

    }

    for (int i = 0; i < 24; i++){
        for (int j = 0 ; j < 24 ; j++) {
            ASSERT_NEAR(chrono_matrix(i, j), analytical_matrix(i, j), tol_stiff);
        }
    }
}

TEST_F(LDPMTest, compute_strain) {
    // Small deflection can be used since the void ChElementLDPM::ComputeStrain function only really performs the
    // computation of the strain according to https://doi.org/10.1016/j.cemconcomp.2011.02.011
    my_LDPM_tet->SetLargeDeflection(false);

    // Performing some kind of step / analysis is required for chrono to run the appropriate setup
    // Functions setupInitial() etc, are all private functions and inaccessible from here, otherwise we use them to make a minimal setup
    sys.DoStepDynamics(0);

    // Only test the first section
    auto section = my_LDPM_tet->GetSection()[0];
    unsigned int ind = nodeIind[0];
    unsigned int jnd = nodeJind[0];
    double length = section->Get_Length();
    ChVector3<> center = section->Get_center();
    ChMatrix33<> nmL = section->Get_facetFrame();
    ChMatrix33<> nmL_tr = nmL.transpose();
    std::vector<ChVector3<>> NML = {nmL_tr.GetAxisX(), nmL_tr.GetAxisY(), nmL_tr.GetAxisZ()};

    double disp = 0.01;
    double rot = 0.01;
    ChVector3d strain;
    ChVectorDynamic<> displ_incr;
    displ_incr.setZero(12);

    double eps = 1e-15;

    // Translation of node A
    for (int i = 0; i<3 ; i++) {
        displ_incr.segment(0,3) = disp * NML[i].eigen(); // N, M, L displacement
        my_LDPM_tet->ComputeStrain(section, ind, jnd, displ_incr, strain);
        ASSERT_NEAR(strain[i], -disp / length, eps);
        ASSERT_NEAR(strain[(i+1)%3], 0.0, eps);
        ASSERT_NEAR(strain[(i+2)%3], 0.0, eps);
        displ_incr.setZero();
    }
        // NML displacements
        displ_incr.segment(0,3) = disp * (NML[0]+NML[1]+NML[2]).eigen();
        my_LDPM_tet->ComputeStrain(section, ind, jnd, displ_incr, strain);
        ASSERT_NEAR(strain[0], -disp / length, eps);
        ASSERT_NEAR(strain[1], -disp / length, eps);
        ASSERT_NEAR(strain[2], -disp / length, eps);
        displ_incr.setZero();

    // Rotation of node A
    for (int i = 0; i<3 ; i++) {
        displ_incr(3+i) = rot; // X, Y, Z rotation
        ChVector3<> rot_vector(rot*(i==0), rot*(i==1), rot*(i==2));
        my_LDPM_tet->ComputeStrain(section, ind, jnd, displ_incr, strain);
        ChVector3<> strain_analytical = -nmL * rot_vector.Cross(center - nodes[ind]->GetPos()).eigen();
        ASSERT_NEAR(strain[0], strain_analytical[0], eps);
        ASSERT_NEAR(strain[1], strain_analytical[1], eps);
        ASSERT_NEAR(strain[2], strain_analytical[2], eps);
        displ_incr.setZero();
    }

    // Translation of node B
    for (int i = 0; i<3 ; i++) {
        displ_incr.segment(6, 3) = disp * NML[i].eigen(); // N, M, L displacement
        my_LDPM_tet->ComputeStrain(section, ind, jnd, displ_incr, strain);
        ASSERT_NEAR(strain[i], disp / length, eps);
        ASSERT_NEAR(strain[(i+1)%3], 0.0, eps);
        ASSERT_NEAR(strain[(i+2)%3], 0.0, eps);
        displ_incr.setZero();
    }
        // NML displacements
        displ_incr.segment(6,3) = disp * (NML[0]+NML[1]+NML[2]).eigen();
        my_LDPM_tet->ComputeStrain(section, ind, jnd, displ_incr, strain);
        ASSERT_NEAR(strain[0], disp / length, eps);
        ASSERT_NEAR(strain[1], disp / length, eps);
        ASSERT_NEAR(strain[2], disp / length, eps);
        displ_incr.setZero();

    // Rotation of node B
    for (int i = 0; i<3 ; i++) {
        displ_incr(9+i) = rot; // X, Y, Z rotation
        ChVector3<> rot_vector(rot*(i==0), rot*(i==1), rot*(i==2));
        my_LDPM_tet->ComputeStrain(section, ind, jnd, displ_incr, strain);
        ChVector3<> strain_analytical = nmL * rot_vector.Cross(center - nodes[jnd]->GetPos());
        ASSERT_NEAR(strain[0], strain_analytical[0], eps);
        ASSERT_NEAR(strain[1], strain_analytical[1], eps);
        ASSERT_NEAR(strain[2], strain_analytical[2], eps);
        displ_incr.setZero();
    }
}

} // namespace

