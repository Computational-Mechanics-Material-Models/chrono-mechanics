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
// Authors: Erol Lale, Jibril B. Coulibaly
// =============================================================================
// Class for LDPM elements:  
//
//  i)   Internal forces
//  ii)  Stiffness matrix
//  iii) Mass matrix  
//  iv)  Body forces
//
// Formulation of the LDPM element can be found in: https://doi.org/10.1016/j.cemconcomp.2011.02.011
// =============================================================================

#include "chrono_ldpm/ChElementLDPM.h"
#include <memory>
#include "chrono/core/ChMatrix.h"
#include "chrono/core/ChVector3.h"
#include "chrono/utils/ChConstants.h"
#include "chrono_ldpm/ChSectionLDPM.h"


namespace chrono {
namespace ldpm {

bool ChElementLDPM::LargeDeflection = false;
ChElementLDPM::MassMatrixType ChElementLDPM::Mmatrix_type = ChElementLDPM::MassMatrixType::CONSISTENT;

ChElementLDPM::ChElementLDPM() : V0(0) {
    nodes.resize(4);
}

ChElementLDPM::~ChElementLDPM() {}

void ChElementLDPM::SetNodes(std::shared_ptr<fea::ChNodeFEAxyzrot> nodeA,
                                     std::shared_ptr<fea::ChNodeFEAxyzrot> nodeB,
                                     std::shared_ptr<fea::ChNodeFEAxyzrot> nodeC,
                                     std::shared_ptr<fea::ChNodeFEAxyzrot> nodeD) {
    nodes[0] = nodeA;
    nodes[1] = nodeB;
    nodes[2] = nodeC;
    nodes[3] = nodeD;
    std::vector<ChVariables*> mvars;
    mvars.push_back(&nodes[0]->Variables());
    mvars.push_back(&nodes[1]->Variables());
    mvars.push_back(&nodes[2]->Variables());
    mvars.push_back(&nodes[3]->Variables());    
    Kmatr.SetVariables(mvars);    
    
}

void ChElementLDPM::Update() {	
    // parent class update:
    ChElementGeneric::Update();
    // always keep updated the rotation matrix A:
    // TODO JBC: Because we inherit from a corotational element, we might consider using
    //           the `disable_corotational` in addition to/instead of the `LargeDeflection` flag we defined here
    if (ChElementLDPM::LargeDeflection) {
        this->UpdateRotation();

        // update projection matrix and compute eigenstrains
        if (macro_strain) {
            for( auto facet:this->GetSection()) {
                facet->ComputeProjectionMatrix();
                facet->ComputeEigenStrain(this->macro_strain);
            }
        } else { // JBC: This might not be necessary since the value of 0.0 in SetupInitial() should remain if macro_strain is nullptr
            for( auto facet:this->GetSection()) {
                facet->Set_nonMechanicStrain(ChVector3d(0.0, 0.0, 0.0));
            }
        }
    }
    // TODO JBC: Erol's original code updates the length of the sections in updateRotation()
}

void ChElementLDPM::ShapeFunctions(ShapeVector& N, double r, double s, double t) {
    N(0) = 1.0 - r - s - t;
    N(1) = r;
    N(2) = s;
    N(3) = t;
}


void ChElementLDPM::GetStateBlock(ChVectorDynamic<>& mD) {
    mD.setZero(this->GetNumCoordsPosLevel());

    // Node displacement in global frame
    auto getDisplacement = [&] (std::shared_ptr<fea::ChNodeFEAxyzrot> node) {
        return node->Frame().GetPos() - node->GetX0().GetPos();
    };

    // Node rotation in global frame
    auto getRotation = [&] (std::shared_ptr<fea::ChNodeFEAxyzrot> node) {
        ChQuaternion<> global_dq = node->Frame().GetRot() * node->GetX0().GetRot().GetConjugate();
        ChVector3d delta_rot_dir;
        double delta_rot_angle;
        global_dq.GetAngleAxis(delta_rot_angle, delta_rot_dir);
        // TODO: We may want to use GetRotVec directly, but the angle is not within [-PI .. PI]
        // TODO: Consider changing GetRotVec() to return within [-PI .. PI]
        return delta_rot_angle * delta_rot_dir;
    };

    //
    // First node
    //
    mD.segment(0, 3) = getDisplacement(nodes[0]).eigen();
    mD.segment(3, 3) = getRotation(nodes[0]).eigen();
    //
    // Second node
    //
    mD.segment(6, 3) = getDisplacement(nodes[1]).eigen();
    mD.segment(9, 3) = getRotation(nodes[1]).eigen();
    //
    // third node
    //
    mD.segment(12, 3) = getDisplacement(nodes[2]).eigen();
    mD.segment(15, 3) = getRotation(nodes[2]).eigen();
    //
    // fourth node
    //
    mD.segment(18, 3) = getDisplacement(nodes[3]).eigen();
    mD.segment(21, 3) = getRotation(nodes[3]).eigen();
}




void ChElementLDPM::GetLatticeStateBlock(unsigned int& ind, unsigned int& jnd, ChVectorDynamic<>& mD) {
    mD.setZero(12);
	mD.segment(0, 6)=DUn_1.segment(ind*6,6);
	mD.segment(6, 6)=DUn_1.segment(jnd*6,6);
}



void ChElementLDPM::GetLatticeField_dt(unsigned int& ind, unsigned int& jnd, ChVectorDynamic<>& mD) {
    mD.setZero(12);    
    //
    //
    // First node
    //
    // displacement dofs:
    mD.segment(0, 3) = (nodes[ind]->Frame().GetPosDt()).eigen();
    // rotational dofs:   	
    mD.segment(3, 3) = (nodes[ind]->Frame().GetAngVelParent()).eigen();
    //
    // Second node
    //
    // displacement dofs:
    mD.segment(6, 3) = (nodes[jnd]->Frame().GetPosDt()).eigen();
    // rotational dofs:   	
    mD.segment(9, 3) = (nodes[jnd]->Frame().GetAngVelParent()).eigen();
        
}



void ChElementLDPM::ComputeStrain(std::shared_ptr<ChSectionLDPM> section, unsigned int& ind, unsigned int& jnd, ChVectorDynamic<>& displ, ChVector3d& mstrain) {
    ChVector3d ui = displ.segment(0,3);
    ChVector3d ri = displ.segment(3,3);
    ChVector3d uj = displ.segment(6,3);
    ChVector3d rj = displ.segment(9,3);
    ChVector3d xc_xi;
    ChVector3d xc_xj;

    double linv = 1.0 / section->Get_Length();
    // For small deflection, use the initial position of the nodes and facet centers
    // to determine the displacement induced by the node rotation
    if (!ChElementLDPM::LargeDeflection) {
        xc_xi = section->Get_center() - this->nodes[ind]->GetX0().GetPos();
        xc_xj = section->Get_center() - this->nodes[jnd]->GetX0().GetPos();
    }
    else { // TODO: Find out how to evolve position of the center for large displacement
            //       Erol's implementation computed the initial A matrix from Get_center and GetX0 as well, so code below does not change the behavior
        xc_xi = section->Get_center() - this->nodes[ind]->GetX0().GetPos(); // TEMPORARY
        xc_xj = section->Get_center() - this->nodes[jnd]->GetX0().GetPos(); // TEMPORARY
    }

    ChMatrix33<double> nmL = section->Get_facetFrame();
    if(ChElementLDPM::LargeDeflection) { // TODO JBC: If the facet orientation were made to coincide with the quaternion, we could save a lot and move this to Update() instead of rotating the initial frame at every step
        ChQuaternion<> q_delta = section->Get_abs_rot() *  section->Get_ref_rot().GetConjugate();
        for (int id=0; id<3; id++){
            nmL.block<1,3>(id, 0) = q_delta.Rotate(nmL.block<1,3>(id, 0)).eigen();
        }
    }

    // Strain
    ChVector3d strain_increment = (uj + rj.Cross(xc_xj) - (ui + ri.Cross(xc_xi))) * linv;
    mstrain = nmL * strain_increment.eigen();
}


double ChElementLDPM::ComputeVolume() {
    return ComputeTetVol(nodes[0]->GetPos(), nodes[1]->GetPos(), nodes[2]->GetPos(), nodes[3]->GetPos());
}

void ChElementLDPM::ComputeStiffnessMatrixGlobal(ChMatrixRef Km) {
    // Assume Km is 24*24 and zero. Not checked here
    for (int iface = 0 ; iface < my_section.size() ; iface++) {
        std::shared_ptr<ChSectionLDPM> facet = my_section[iface];
        unsigned int ind = facetNodeNums(iface, 0);
        unsigned int jnd = facetNodeNums(iface, 1);

        double normal_stiff = facet->Get_material()->Get_E0() * facet->Get_area() / facet->Get_Length();
        double tangent_stiff = normal_stiff * facet->Get_material()->Get_alpha();
        Eigen::DiagonalMatrix<double, 3> K_diag(normal_stiff, tangent_stiff, tangent_stiff);

        ChMatrix33<double> nmL = facet->Get_facetFrame();
        if(ChElementLDPM::LargeDeflection) { // TODO JBC: If the facet orientation were made to coincide with the quaternion, we could save a lot and move this to Update() instead of rotating the initial frame at every step
            ChQuaternion<> q_delta = facet->Get_abs_rot() * facet->Get_ref_rot().GetConjugate();
            for (int id=0; id<3; id++){
                nmL.block<1,3>(id, 0) = q_delta.Rotate(nmL.block<1,3>(id, 0)).eigen();
            }
        }
        ChMatrix33<double> nmL_tr = nmL.transpose();
        ChMatrix33<> K_material = nmL_tr * K_diag * nmL;

        // For small deflection, use the initial position of the nodes and facet centers
        // to determine the displacement induced by the node rotation
        ChVector3d xc_xi;
        ChVector3d xc_xj;
        if (!ChElementLDPM::LargeDeflection) {
            xc_xi = facet->Get_center() - this->nodes[ind]->GetX0().GetPos();
            xc_xj = facet->Get_center() - this->nodes[jnd]->GetX0().GetPos();
        }
        else { // TODO: Find out how to evolve position of the center for large displacement
               //       Erol's implementation computed the initial A matrix from Get_center and GetX0 as well, so code below does not change the behavior
            xc_xi = facet->Get_center() - this->nodes[ind]->GetX0().GetPos(); // TEMPORARY
            xc_xj = facet->Get_center() - this->nodes[jnd]->GetX0().GetPos(); // TEMPORARY
        }

        ChMatrix33<> Ai_cross; // 3x3 right sub-block of Ai matrix for skew-symmetric cross-product vector
        ChMatrix33<> Aj_cross; // 3x3 right sub-block of Aj matrix for skew-symmetric cross-product vector
        Ai_cross << 0.0      ,  xc_xi[2], -xc_xi[1],
                    -xc_xi[2],  0.0     ,  xc_xi[0],
                     xc_xi[1], -xc_xi[0],  0.0     ;
        Aj_cross << 0.0      ,  xc_xj[2], -xc_xj[1],
                    -xc_xj[2],  0.0     ,  xc_xj[0],
                     xc_xj[1], -xc_xj[0],  0.0     ;

        ChMatrix33<> K_material_Id_Ai = K_material * Ai_cross;
        ChMatrix33<> K_material_Ai_Id = K_material_Id_Ai.transpose();
        ChMatrix33<> K_material_Ai_Ai = K_material_Ai_Id * Ai_cross;
        ChMatrix33<> K_material_Id_Aj = K_material * Aj_cross;
        ChMatrix33<> K_material_Aj_Id = K_material_Id_Aj.transpose();
        ChMatrix33<> K_material_Aj_Aj = K_material_Aj_Id * Aj_cross;
        ChMatrix33<> K_material_Ai_Aj = K_material_Ai_Id * Aj_cross;
        ChMatrix33<> K_material_Aj_Ai = K_material_Ai_Aj.transpose();

        // dFi / dQi
        Km.block<3,3>(ind*6  , ind*6)   +=  K_material;
        Km.block<3,3>(ind*6  , ind*6+3) +=  K_material_Id_Ai;
        Km.block<3,3>(ind*6+3, ind*6)   +=  K_material_Ai_Id;
        Km.block<3,3>(ind*6+3, ind*6+3) +=  K_material_Ai_Ai;
        // dFi / dQj
        Km.block<3,3>(ind*6  , jnd*6)   += -K_material;
        Km.block<3,3>(ind*6  , jnd*6+3) += -K_material_Id_Aj;
        Km.block<3,3>(ind*6+3, jnd*6)   += -K_material_Ai_Id;
        Km.block<3,3>(ind*6+3, jnd*6+3) += -K_material_Ai_Aj;
        // dFj / dQi
        Km.block<3,3>(jnd*6  , ind*6)   += -K_material;
        Km.block<3,3>(jnd*6  , ind*6+3) += -K_material_Id_Ai;
        Km.block<3,3>(jnd*6+3, ind*6)   += -K_material_Aj_Id;
        Km.block<3,3>(jnd*6+3, ind*6+3) += -K_material_Aj_Ai;
        // dFj / dQj
        Km.block<3,3>(jnd*6  , jnd*6)   +=  K_material;
        Km.block<3,3>(jnd*6  , jnd*6+3) +=  K_material_Id_Aj;
        Km.block<3,3>(jnd*6+3, jnd*6)   +=  K_material_Aj_Id;
        Km.block<3,3>(jnd*6+3, jnd*6+3) +=  K_material_Aj_Aj;
    }
}

void ChElementLDPM::SetupInitial(ChSystem* system) {
    //
    //this->mysystem=system;
    // set initial orientation of facets
    unsigned int iface=0;    
    for( auto facet:this->GetSection()){ 
    	unsigned int ind=facetNodeNums(iface,0);
    	unsigned int jnd=facetNodeNums(iface,1);
	    
        // TODO JBC: the quaternion should be defined at the tetrahedron level, this would save a lot of compute and memory!
	    ChMatrix33<> A0;
	    ChVector3d mXele = nodes[jnd]->GetX0().GetPos() - nodes[ind]->GetX0().GetPos();
	    ChVector3d myele =
		(nodes[ind]->GetX0().GetRotMat().GetAxisY() + nodes[jnd]->GetX0().GetRotMat().GetAxisY()).GetNormalized();
	    A0.SetFromAxisX(mXele, myele); // If the Y-axes are aligned with the node-to-node direction, , or cancel out, Gram-Schmidt will fail and use an arbitrary direction
	    facet->Set_ref_rot( A0.GetQuaternion() );
	    ///
		///		
		double length = (this->GetTetrahedronNode(ind)->GetX0().GetPos() - this->GetTetrahedronNode(jnd)->GetX0().GetPos()).Length();		
		facet->Set_Length(length);
	    ///
	    ///
	    auto statev= facet->Get_StateVar();
	    statev.resize(16); // TODO JBC: this is sized to 16 but the highest it is ever accessed in (15) inside ChMaterialVECT
	    statev.setZero();
	    facet->Set_StateVar(statev); 
		///
		///
		if(this->macro_strain){
			facet->ComputeProjectionMatrix();			
			facet->ComputeEigenStrain(this->macro_strain);
		} else {
            facet->Set_nonMechanicStrain(ChVector3d(0.0, 0.0, 0.0));
        }
	    iface++;
     }

	ChVectorDynamic<> displ(24);
    this->GetStateBlock(displ);
	Un_1=displ;

    this->UpdateRotation(); // TODO JBC: not sure why this is needed here

    this->SetInitialVolume(ComputeVolume());
}

void ChElementLDPM::UpdateRotation() {
    //
    //    
    unsigned int iface=0;    
    for( auto facet:this->GetSection()){
    	    unsigned int ind=facetNodeNums(iface,0);
    	    unsigned int jnd=facetNodeNums(iface,1);   

        ChMatrix33<> A0(facet->Get_ref_rot());
        ChMatrix33<> Aabs;
        ChQuaternion<> q_lattice_abs_rot;
        if (!LargeDeflection) {
        Aabs = A0;
        q_lattice_abs_rot =facet->Get_ref_rot();
        } else {
            ChVector3d mXele = nodes[jnd]->Frame().GetPos() - nodes[ind]->Frame().GetPos();
            ChVector3d myele = (nodes[ind]->Frame().GetRotMat().GetAxisY() + nodes[jnd]->Frame().GetRotMat().GetAxisY()).GetNormalized();
            // TODO JBC: This looks fragile and relies on assumptions and current behavior of other parts of the code
            //           - If the Y-axes are aligned with the node-to-node direction, or cancel out, Gram-Schmidt will fail and use an arbitrary direction
            //              This arbitrary direction is repeatable, that is the only thing that guarantee this works, otherwise, the initial quaternion determined in SetupInitial() might be different
            //              Any change to the fallback behavior of ChMatrix33::SetFromAxisX() will break this code !
            //           - We could compute the rigid body motion of the tegrahedron, rather than compute and store an individual quaternion for all 12 facet to facets.
            //              This is what is done in `ChElementTetraCorot_4::UpdateRotation()` using polar decompositoin, which seems better suited
            //              The difficulty is that the N-axis might not be aligned with the node-to-node direction because of the deformation, albeit small
            Aabs.SetFromAxisX(mXele, myele);
            q_lattice_abs_rot = Aabs.GetQuaternion();

            double length = (this->GetTetrahedronNode(ind)->GetX0().GetPos() - this->GetTetrahedronNode(jnd)->GetX0().GetPos()).Length(); // TODO JBC: this has nothing to do with rotation	and could be taken to Update()
            facet->Set_Length(length);
        }
        facet->Set_abs_rot(q_lattice_abs_rot);
        iface++;
    }
}

void ChElementLDPM::ComputeKRMmatricesGlobal(ChMatrixRef H, double Kfactor, double Rfactor, double Mfactor) {
    assert((H.rows() == 24) && (H.cols() == 24));    
    // For K stiffness matrix and R damping matrix:
    double RayleighDampingK=this->GetFacetI(0)->Get_material()->GetRayleighDampingK();
    double RayleighDampingM=this->GetFacetI(0)->Get_material()->GetRayleighDampingM();

    if (Kfactor || Rfactor) {
        // LDPM currently uses the initial elastic stiffness matrix
        ChMatrixNM<double, 24, 24> Km;
        Km.setZero();
        ComputeStiffnessMatrixGlobal(Km);
        double mkfactor = Kfactor + Rfactor * RayleighDampingK;
        H.block(0, 0, 24, 24) = mkfactor*Km;
    }

    if (Mfactor || Rfactor) {      	    
        ChMatrixDynamic<> Mloc(24, 24); 
	this->ComputeMmatrixGlobal(Mloc);	
	double amfactor = Mfactor + Rfactor * RayleighDampingM;
	H.block(0, 0, 24, 24)+= amfactor*Mloc;
    }

    //***TO DO*** better per-node lumping, or 12x12 consistent mass matrix.
}

void ChElementLDPM::ComputeInternalForces(ChVectorDynamic<>& Fi) {
    assert(Fi.size() == 24);
    //assert(section);

    ChVectorDynamic<> displ(24);
    this->GetStateBlock(displ);
    DUn_1=displ-Un_1;
    Un_1=displ;

    Fi.setZero();
    //
    // compute volumetric strain
    double V = this->ComputeVolume();	
    double epsV = (V - V0) / V0 * CH_1_3; // This 1/3 factor is correct! Original paper (https://doi.org/10.1016/j.cemconcomp.2011.02.011) without the 1/3 is a typo

    for (int iface = 0 ; iface < my_section.size() ; iface++) {
        std::shared_ptr<ChSectionLDPM> facet = my_section[iface];
        unsigned int ind = facetNodeNums(iface, 0);
        unsigned int jnd = facetNodeNums(iface, 1);

        double area = facet->Get_area();
        double length = facet->Get_Length();
                                        
        ChMatrix33<double> nmL=facet->Get_facetFrame();
        if(ChElementLDPM::LargeDeflection) { // TODO JBC: If the facet orientation were made to coincide with the quaternion, we could save a lot and move this to Update() instead of rotating the initial frame at every step
            ChQuaternion<> q_delta = facet->Get_abs_rot() *  facet->Get_ref_rot().GetConjugate();
            for (int id=0; id<3; id++){
                nmL.block<1,3>(id, 0) = q_delta.Rotate(nmL.block<1,3>(id, 0)).eigen();
            }
        }
        ChMatrix33<double> nmL_tr = nmL.transpose();
        //
        // Get Stress values at facet center
        //
        ChVector3d mstress;
        ChVector3d dmstrain;
        ChVectorDynamic<> statev;
        ChVectorDynamic<> displ_facet_nodes(12);
        this->GetLatticeStateBlock(ind, jnd, displ_facet_nodes);	
        this->ComputeStrain(facet, ind, jnd, displ_facet_nodes, dmstrain);
        //std::cout<<"strain_INC: "<<dmstrain(0)<<"\t"<<dmstrain(1)<<"\t"<<dmstrain(2)<<"\t";
        statev=facet->Get_StateVar();	

        // TODO JBC: not sure why 2 lines below were kept commented. If not used, get rid of them
        //if (this->macro_strain)
        //	facet->ComputeEigenStrain(this->macro_strain);
        ChVector3d nonMechanicalStrain=facet->Get_nonMechanicStrain();
        facet->Get_material()->ComputeStress( dmstrain, nonMechanicalStrain, length,  epsV, statev, mstress, area);
        facet->Set_StateVar(statev);	

        ChVector3d force = area * (nmL_tr * mstress);
        ChVector3d xc_xi;
        ChVector3d xc_xj;
        // For small deflection, use the initial position of the nodes and facet centers
        // to determine the displacement induced by the node rotation
        if (!ChElementLDPM::LargeDeflection) {
            xc_xi = facet->Get_center() - this->nodes[ind]->GetX0().GetPos();
            xc_xj = facet->Get_center() - this->nodes[jnd]->GetX0().GetPos();
        }
        else { // TODO: Find out how to evolve position of the center for large displacement
               //       Erol's implementation computed the initial A matrix from Get_center and GetX0 as well, so code below does not change the behavior
            xc_xi = facet->Get_center() - this->nodes[ind]->GetX0().GetPos(); // TEMPORARY
            xc_xj = facet->Get_center() - this->nodes[jnd]->GetX0().GetPos(); // TEMPORARY
        }

        Fi.segment(ind*6+0,3) +=  force.eigen();
        Fi.segment(ind*6+3,3) +=  xc_xi.Cross(force).eigen();
        Fi.segment(jnd*6+0,3) += -force.eigen();
        Fi.segment(jnd*6+3,3) += -xc_xj.Cross(force).eigen();
    }
}


void ChElementLDPM::ComputeNodalMass() {
    // All sections materials have the same density, use first section
    // Assign same mass to all nodes (I don't think this is used by FEA so this approximation is ok)
    double nodal_mass = V0 * my_section[0]->Get_material()->Get_density() * 0.25;
    nodes[0]->m_TotalMass += nodal_mass;
    nodes[1]->m_TotalMass += nodal_mass;
    nodes[2]->m_TotalMass += nodal_mass;
    nodes[3]->m_TotalMass += nodal_mass;
}

void ChElementLDPM::LoadableGetStateBlockPosLevel(int block_offset, ChState& mD) {
    mD.segment(block_offset + 0, 3) = nodes[0]->GetPos().eigen();
    mD.segment(block_offset + 3, 4) = nodes[0]->GetRot().eigen();
    //
    mD.segment(block_offset + 7, 3) = nodes[1]->GetPos().eigen();
    mD.segment(block_offset + 10, 4) = nodes[1]->GetRot().eigen();
    //
    mD.segment(block_offset + 14, 3) = nodes[2]->GetPos().eigen();
    mD.segment(block_offset + 17, 4) = nodes[2]->GetRot().eigen();
    //
    mD.segment(block_offset + 21, 3) = nodes[3]->GetPos().eigen();
    mD.segment(block_offset + 24, 4) = nodes[3]->GetRot().eigen();
}

void ChElementLDPM::LoadableGetStateBlockVelLevel(int block_offset, ChStateDelta& mD) {
    mD.segment(block_offset + 0, 3) = nodes[0]->GetPosDt().eigen();
    mD.segment(block_offset + 3, 3) = nodes[0]->GetAngVelLocal().eigen();
    //
    mD.segment(block_offset + 6, 3) = nodes[1]->GetPosDt().eigen();
    mD.segment(block_offset + 9, 3) = nodes[1]->GetAngVelLocal().eigen();
    //
    mD.segment(block_offset + 12, 3) = nodes[2]->GetPosDt().eigen();
    mD.segment(block_offset + 15, 3) = nodes[2]->GetAngVelLocal().eigen();
    //
    mD.segment(block_offset + 18, 3) = nodes[3]->GetPosDt().eigen();
    mD.segment(block_offset + 21, 3) = nodes[3]->GetAngVelLocal().eigen();
}



void ChElementLDPM::LoadableStateIncrement(const unsigned int off_x,
                                                   ChState& x_new,
                                                   const ChState& x,
                                                   const unsigned int off_v,
                                                   const ChStateDelta& Dv) {
    for (int i = 0; i < 4; ++i) {
        nodes[i]->NodeIntStateIncrement(off_x + i * 7, x_new, x, off_v + i * 6, Dv);        
    }
}

void ChElementLDPM::LoadableGetVariables(std::vector<ChVariables*>& mvars) {
    for (int i = 0; i < nodes.size(); ++i)
        mvars.push_back(&this->nodes[i]->Variables());
}

void ChElementLDPM::ComputeNF(const double U,
                                      const double V,
                                      const double W,
                                      ChVectorDynamic<>& Qi,
                                      double& detJ,
                                      const ChVectorDynamic<>& F,
                                      ChVectorDynamic<>* state_x,
                                      ChVectorDynamic<>* state_w) {
    // evaluate shape functions (in compressed vector), btw. not dependent on state
    // note: U,V,W in 0..1 range, thanks to IsTetrahedronIntegrationNeeded() {return true;}
    ShapeVector N;
    this->ShapeFunctions(N, U, V, W);

    detJ = 6 * V0;

    Qi(0) = N(0) * F(0);
    Qi(1) = N(0) * F(1);
    Qi(2) = N(0) * F(2);
    Qi(3) = N(1) * F(0);
    Qi(4) = N(1) * F(1);
    Qi(5) = N(1) * F(2);
    Qi(6) = N(2) * F(0);
    Qi(7) = N(2) * F(1);
    Qi(8) = N(2) * F(2);
    Qi(9) = N(3) * F(0);
    Qi(10) = N(3) * F(1);
    Qi(11) = N(3) * F(2);
}


double ChElementLDPM::ComputeTetVol(ChVector3d p1, ChVector3d p2, ChVector3d p3, ChVector3d p4){
    double tetvol=0;
    tetvol =  p2[0] * (p3[1] * p4[2] - p4[1] * p3[2]) - p3[0] * (p2[1] * p4[2] - p4[1] * p2[2]) + p4[0] * (p2[1] * p3[2] - p3[1] * p2[2]);
    tetvol -= p3[0] * (p4[1] * p1[2] - p1[1] * p4[2]) - p4[0] * (p3[1] * p1[2] - p1[1] * p3[2]) + p1[0] * (p3[1] * p4[2] - p4[1] * p3[2]);
    tetvol += p4[0] * (p1[1] * p2[2] - p2[1] * p1[2]) - p1[0] * (p4[1] * p2[2] - p2[1] * p4[2]) + p2[0] * (p4[1] * p1[2] - p1[1] * p4[2]);
    tetvol -= p1[0] * (p2[1] * p3[2] - p3[1] * p2[2]) - p2[0] * (p1[1] * p3[2] - p3[1] * p1[2]) + p3[0] * (p1[1] * p2[2] - p2[1] * p1[2]);
    tetvol = abs(tetvol) / 6.0;
	return tetvol;
}

ChMatrixNM<double,6,6> ChElementLDPM::ComputeSubTetMassMatrix(std::shared_ptr<ChSectionLDPM> section, ChVector3d pN, ChVector3d pC, ChVector3d pA, ChVector3d pB){
    double tetmass = section->Get_material()->Get_density() * this->ComputeTetVol(pN, pC, pA, pB);
    ChVector3d pG = 0.25 * (pN + pC + pA + pB);
    ChVector3d X = {pA[0] - pN[0], pB[0] - pN[0], pC[0] - pN[0]};
    ChVector3d Y = {pA[1] - pN[1], pB[1] - pN[1], pC[1] - pN[1]};
    ChVector3d Z = {pA[2] - pN[2], pB[2] - pN[2], pC[2] - pN[2]};
    // first moment of area
    double Sx = (pG[0] - pN[0]) * tetmass;
    double Sy = (pG[1] - pN[1]) * tetmass;
    double Sz = (pG[2] - pN[2]) * tetmass;
    // Moment of inertia
    // Calculation in expanded bilinear form:
    // coef = {{0.1, 0.05, 0.05},{0.05, 0.1, 0.05},{0.05, 0.05, 0.1}};
    // I_XY = (X^T * coef * Y) * tetmass (for all "X" and "Y")
    // Simplified form below to avoid the overhead of matrix calculations and conversions
    ChVector3d coefX = {0.1 * X[0] + 0.05 * (X[1] + X[2]), 0.1 * X[1] + 0.05 * (X[0] + X[2]), 0.1 * X[2] + 0.05 * (X[0] + X[1])};
    ChVector3d coefY = {0.1 * Y[0] + 0.05 * (Y[1] + Y[2]), 0.1 * Y[1] + 0.05 * (Y[0] + Y[2]), 0.1 * Y[2] + 0.05 * (Y[0] + Y[1])};
    ChVector3d coefZ = {0.1 * Z[0] + 0.05 * (Z[1] + Z[2]), 0.1 * Z[1] + 0.05 * (Z[0] + Z[2]), 0.1 * Z[2] + 0.05 * (Z[0] + Z[1])};
    double Ixx = X.Dot(coefX) * tetmass;
    double Iyy = Y.Dot(coefY) * tetmass;
    double Izz = Z.Dot(coefZ) * tetmass;
    double Ixy = X.Dot(coefY) * tetmass;
    double Ixz = X.Dot(coefZ) * tetmass;
    double Iyz = Y.Dot(coefZ) * tetmass;
    ChMatrixNM<double,6,6> MN;
    MN.setZero();
    MN(0,0) = tetmass; MN(0,4) =  Sz; MN(0,5) = -Sy;
    MN(1,1) = tetmass; MN(1,3) = -Sz; MN(1,5) =  Sx;
    MN(2,2) = tetmass; MN(2,3) =  Sy; MN(2,4) = -Sx;
    MN(3,1) = -Sz; MN(3,2) =  Sy; MN(3,3) = Iyy+Izz; MN(3,4) =    -Ixy; MN(3,5) =    -Ixz;
    MN(4,0) =  Sz; MN(4,2) = -Sx; MN(4,3) =    -Ixy; MN(4,4) = Ixx+Izz; MN(4,5) =    -Iyz;
    MN(5,0) = -Sy; MN(5,1) =  Sx; MN(5,3) =    -Ixz; MN(5,4) =    -Iyz; MN(5,5) = Ixx+Iyy;
	return MN;
}


void ChElementLDPM::ComputeMmatrixGlobal(ChMatrixRef M) {
    M.setZero();
    if (Mmatrix_type == MassMatrixType::LUMPED) {
        // This hydrostatic LUMPED mass matrix is invariant by rigid-body rotation of the tetrahedron
        // It is very cheap to compute and is preferred over the CONSISTENT one for quasi-static analyses.

        // Translational DOFs
        // TODO JBC: same mass on each node: this could be refined based e.g., on volume fractions, but OK for now.
        double nodal_mass = V0 * my_section[0]->Get_material()->Get_density() * 0.25; // All sections materials have the same density, use first section
        M(0,0)   = M(1,1)   = M(2,2)   = nodal_mass; // Node 1
        M(6,6)   = M(7,7)   = M(8,8)   = nodal_mass; // Node 2
        M(12,12) = M(13,13) = M(14,14) = nodal_mass; // Node 3
        M(18,18) = M(19,19) = M(20,20) = nodal_mass; // Node 4

        // Rotational DOFs
        // The moment of inertia computed from the lumped mass at the nodes alone is already larger than
        // the true moment of inertia of the tetrahedron. Adding rotational lumped mass would make it worse.
        // This is a know issue for beams (see e.g., for beams https://quickfem.com/wp-content/uploads/IFEM.Ch31.pdf)
        // and a small value is taken so that the mass matrix is not singular, but not ill-defined either.
        // For Euler-Bernoulli beams, (see ChBeamSectionEuler::JzzJyy_factor) Chrono uses a value of 1/500.
        // We use the same value here.
        double factor = 1.0 / 500;
        M(3,3)   = M(4,4)   = M(5,5)   = nodal_mass * factor; // Node 1
        M(9,9)   = M(10,10) = M(11,11) = nodal_mass * factor; // Node 2
        M(15,15) = M(16,16) = M(17,17) = nodal_mass * factor; // Node 3
        M(21,21) = M(22,22) = M(23,23) = nodal_mass * factor; // Node 4

    } else if (Mmatrix_type == MassMatrixType::CONSISTENT) {
        // This CONSISTENT mass matrix is very expensive to compute and is not preferred for performance.
        // That is because the HHT timestepper currently computes the inertial term in residual using:
        //      integrable2->LoadResidual_Mv(R, Anew, -1 / (1 + alpha));  // -1/(1+alpha)*M*a_new
        // which requests recalculation of the mass matrix from every element (even for JacobianUpdate::NEVER).
        // TODO JBC: Storing the sparse, system-wide mass matrix in HHT is the way to go, but that requires infrastructure changes to be discussed with Alessandro Tasora.
        // Storing the 24x24 dense matrix at the element level would be too memory-intensive for tetrahedral meshes, where 1 vertex can be shared by more than 10 tetrahedra in 3D
        // so for the time being, this matrix is re-computed everytime Chrono asks for it...
        auto vertices=this->V_vert_nodes;
        unsigned int iface=0;
        for(auto verts: vertices){	
            unsigned int ind=facetNodeNums(iface,0);
            unsigned int jnd=facetNodeNums(iface,1);
            //
            ChVector3d pNA=this->nodes[ind]->GetX0().GetPos();
            ChVector3d pNB=this->nodes[jnd]->GetX0().GetPos();
            //
            auto pC=verts[0];
            auto pA=verts[1];
            auto pB=verts[2];		
            //
            ChMatrixNM<double,6,6> mA=this->ComputeSubTetMassMatrix(this->GetFacetI(iface), pNA, pC, pA, pB);
            ChMatrixNM<double,6,6> mB=this->ComputeSubTetMassMatrix(this->GetFacetI(iface), pNB, pC, pA, pB);
            M.block<6,6>(ind*6,ind*6)+=mA;
            M.block<6,6>(jnd*6,jnd*6)+=mB;
            //
            iface++;
        }

        if (ChElementLDPM::LargeDeflection) {
            // TODO JBC: the mass matrix must be rotated into the current spatial orientation of the tetrahedron
            //           I think instead of rotating the 24*24 matrix here, we should rotate the vertices into the current configuration before building the matrix above. TBD
        }
    }
}


void ChElementLDPM::EleIntLoadLumpedMass_Md(ChVectorDynamic<>& Md, double& error, const double c){
        ChMatrixDynamic<> Mloc(GetNumCoordsPosLevel(), GetNumCoordsPosLevel()); 
	this->ComputeMmatrixGlobal(Mloc);			
	
	ChVectorDynamic<> dMi = c * Mloc.diagonal();
    
    	error = Mloc.sum() - Mloc.diagonal().sum();	
	
    	int stride = 0;
    	for (int in = 0; in < GetNumNodes(); in++) {
		int node_dofs = GetNodeNumCoordsPosLevelActive(in);		
		if (!GetNode(in)->IsFixed()){		    
		    Md.segment(GetNode(in)->NodeGetOffsetVelLevel(), node_dofs) += dMi.segment(stride, node_dofs);
		}
		stride += GetNodeNumCoordsPosLevel(in);
    	}      	
    	         
    
}


ChMatrixNM<double, 1, 9> ChElementLDPM::ComputeMacroStressContribution(){
	ChMatrixNM<double, 1, 9> macro_stress;	  
    macro_stress.setZero();	
	for (auto facet:this->GetSection()){		
		double length=facet->Get_Length();		 	
		double area =facet->Get_area();	
		auto statev= facet->Get_StateVar();		
		macro_stress +=(facet->GetProjectionMatrix()).transpose()*statev.segment(3,3)*area*length;		
	}
	return macro_stress;
}


chrono::ChMatrixNM<int,12,2> fill(){
chrono::ChMatrixNM<int,12,2> facetNodeNums;
facetNodeNums << 0, 1, 0, 1, 0, 2, 0, 2, 0, 3, 0, 3, 1, 2, 1, 2, 1, 3, 1, 3, 2, 3, 2, 3;
return facetNodeNums;
}

chrono::ChMatrixNM<int,12,2> ChElementLDPM::facetNodeNums=fill();



}  // end namespace ldpm
}  // end namespace chrono
