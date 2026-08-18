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
// Authors: Erol Lale,  Wisdom Akpan, Ke Yu, Jibril B. Coulibaly
// =============================================================================
// Material class for LDPM and CSL elements 
//
// A description of the material parameter can be found in: https://doi.org/10.1016/j.cemconcomp.2011.02.011
// =============================================================================

#include "chrono_wood/ChWoodMaterialVECT.h"
#include "chrono/core/ChMatrix.h"
#include "chrono/core/ChVector3.h"

namespace chrono {
namespace wood {

// Construct an isotropic material.

ChWoodMaterialVECT::ChWoodMaterialVECT(double rho,  // material density
                                       double E0,    // Mesoscale Young's modulus
                                       double alpha,   // Mesoscale Poisson ratio                                       
                                       double sigmat, double sigmac, double sigmas,
                                       double nt, double lt, double rs)
    : m_rho(rho), m_E0(E0), m_alpha(alpha), m_sigmat(sigmat), m_sigmas(sigmas), m_sigmac(sigmac), m_nt(nt), m_lt(lt), m_rs(rs) {}

ChWoodMaterialVECT::ChWoodMaterialVECT() {}

ChWoodMaterialVECT::~ChWoodMaterialVECT() {}

// statec(0): normal N strain
// statec(1): shear M strain
// statec(2): shear L strain
// statec(3): normal N stress
// statec(4): shear M stress
// statec(5): shear L stress
// statec(6): maximum normal N strain
// statec(7): maximum shear T strain
// statec(8): effective strain
// statec(9): effective stress
// statec(10): internal work
// statec(11): crack opening

void ChWoodMaterialVECT::ComputeStress(ChVector3d& strain_increment, ChVector3d& curvature_increment, ChVector3d& eigenstrain, double &len, const ChVectorDynamic<>& statev_old, ChVectorDynamic<>& statev_new,  double& area, double& width, double& height, double& random_field, ChVector3d& mstress, ChVector3d& mcouple) {
	double E0=this->Get_E0();
	double alpha=this->Get_alpha();	
	//if (random_field)
	//	E0=E0*random_field;

    ChVector3d strain(statev_old.segment(0,3) + strain_increment.eigen() - eigenstrain.eigen()); // TODO JBC: this is the net strain
    // TODO JBC: LDPM stores the netstrain as a state variable. CBL only stores the strain, so we cannot compute net strain increment because the old eigenstrain is lost.
    //           We compute the total strain increment here, not the net strain increment, until this is resolved
    ChVector3d dstrain = strain_increment;
    ChVector3d curvature(statev_old.segment(12,3) + curvature_increment.eigen()); // TODO: no eigencurvature ?
    ChVector3d dcurvature = curvature_increment;

	double epsQ, epsQN;
	double epsT;
	double h2=height*height;		
	double w2=width*width;


	double sigmat = this->Get_sigmat();
	double sigmac = this->Get_sigmac();
	double epsN = strain[0];
	double kappa;


    if (epsN > 0.0) {
        kappa = 1.0;
    } else if (epsN < 0.0) {
        kappa = -1.0;
    } else {
        if (sigmac != sigmat) {
            double val = sigmac / sigmat - 1.0;

            if (val > 0.0)
                kappa = 1.0;
            else if (val < 0.0)
                kappa = -1.0; 
            else
                kappa = 0.0; 
        } else {
            kappa = 1.0;
        }
    }
	

	double multiplier=this->GetCoupleMultiplier()/3.;
	epsQN = kappa * std::sqrt(strain[0] * strain[0] + multiplier*( curvature[1]*curvature[1]*w2 + curvature[2]*curvature[2]*h2));
	epsT = std::sqrt(strain[1] * strain[1] + strain[2] * strain[2]+multiplier*curvature[0]*curvature[0]*(w2+h2));
	epsQ = std::sqrt(epsQN * epsQN + alpha * epsT * epsT);  
	
	statev_new(6) = std::max(statev_old(6), epsQN);
	statev_new(7) = std::max(statev_old(7), epsT);

    double eps_max = std::sqrt(statev_new(6) * statev_new(6) + alpha * statev_new(7) * statev_new(7)); 
	
	
	if (!GetElasticAnalysisFlag()) {
	///
	// INELASTIC ANALYSIS
	//
	if (epsQ != 0) {
			
		double strsQ = StressBC(strain, random_field, len, epsQ, epsQN, epsT, statev_old, eps_max);
		mstress[0] = strsQ * strain[0] / epsQ;
		mstress[1] = alpha * strsQ * strain[1] / epsQ;
		mstress[2] = alpha * strsQ * strain[2] / epsQ;

		mcouple[0]=multiplier*alpha*strsQ*curvature[0]*(w2+h2)/epsQ;
		mcouple[1]=multiplier*strsQ*curvature[1]*w2/epsQ;
		mcouple[2]=multiplier*strsQ*curvature[2]*h2/epsQ;						

		double Wint = len * area * (((mstress[0] + statev_old(3)) / 2.0) * dstrain[0] + ((mstress[1] + statev_old(4)) / 2.0) * dstrain[1] + ((mstress[2] + statev_old(5)) / 2.0) * dstrain[2] + ((mcouple[0] + statev_old(15)) / 2.0) * dcurvature[0] + ((mcouple[1] + statev_old(16)) / 2.0) * dcurvature[1] + ((mcouple[2] + statev_old(17)) / 2.0) * dcurvature[2]);
		
		double strainN_in = strain[0] - mstress[0] / E0;
		double strainM_in = strain[1] - mstress[1] / (alpha * E0);
		double strainL_in = strain[2] - mstress[2] / (alpha * E0);

		double curvatureN_in = curvature[0] - mcouple[0] / (multiplier * alpha * E0 * (w2 + h2));
		double curvatureM_in = curvature[1] - mcouple[1] / (multiplier * E0 * w2);
		double curvatureL_in = curvature[2] - mcouple[2] / (multiplier * E0 * h2);

		double w_N = len * strainN_in;
		double w_M = len * strainM_in;
		double w_L = len * strainL_in;

		double w_kN = len * std::sqrt(multiplier * (w2 + h2)) * curvatureN_in;
		double w_kM = len * std::sqrt(multiplier * w2) * curvatureM_in;
		double w_kL = len * std::sqrt(multiplier * h2) * curvatureL_in;

		double w = std::sqrt(w_N  * w_N  +  w_M  * w_M  +  w_L  * w_L  +  w_kN * w_kN + w_kM * w_kM + w_kL * w_kL);

		double sgn_strsQ = (mstress[0] >= 0.0) ? 1.0 : (mstress[0] < 0.0) ? -1.0: 0;
		
		double dstrainN_in = dstrain[0] - (mstress[0] - statev_old(3)) / E0;
		double dstrainM_in = dstrain[1] - (mstress[1] - statev_old(4)) / (alpha * E0);
		double dstrainL_in = dstrain[2] - (mstress[2] - statev_old(5)) / (alpha * E0);

		double dcurvatureN_in = dcurvature[0] - (mcouple[0] - statev_old(15)) / (multiplier * alpha * E0 * (w2 + h2));
		double dcurvatureM_in = dcurvature[1] - (mcouple[1] - statev_old(16)) / (multiplier * E0 * w2);
		double dcurvatureL_in = dcurvature[2] - (mcouple[2] - statev_old(17)) / (multiplier * E0 * h2);
		

		double DE = len * area * (0.5 * (mstress[0]  + statev_old(3))  * dstrainN_in +  0.5 * (mstress[1]  + statev_old(4))  * dstrainM_in +  0.5 * (mstress[2]  + statev_old(5))  * dstrainL_in +  0.5 * (mcouple[0] + statev_old(15)) * dcurvatureN_in +   0.5 * (mcouple[1] + statev_old(16)) * dcurvatureM_in +   0.5 * (mcouple[2] + statev_old(17)) * dcurvatureL_in );


		statev_new(0) = strain[0] + eigenstrain[0]; // TODO JBC: it is a bug to only update those inside this code branch.
		statev_new(1) = strain[1] + eigenstrain[1]; //           If you unload to exactly zero (like in the unit tests)
		statev_new(2) = strain[2] + eigenstrain[2]; //           you go to the `else` branch and the state variables do not get updated
		statev_new(3) = mstress[0];                  //           while a loading step did happen!
		statev_new(4) = mstress[1];                  //           I am leaving this as is for now because we will refactor this in depth soon!
		statev_new(5) = mstress[2];                  //           The current fix aims to retrieve the "old" behavior
		statev_new(12) = curvature[0];
        statev_new(13) = curvature[1];
        statev_new(14) = curvature[2];
            //
        statev_new(15) = mcouple[0];
        statev_new(16) = mcouple[1];
        statev_new(17) = mcouple[2];
		double epsQN_new = kappa * std::sqrt(statev_new(0) * statev_new(0) + multiplier*( statev_new(13)*statev_new(13)*w2 + statev_new(14)*statev_new(14)*h2));
		double epsT_new = std::sqrt(statev_new(1) * statev_new(1) + statev_new(2) * statev_new(2)+multiplier*statev_new(12)*statev_new(12)*(w2+h2));
		statev_new(8) =  std::sqrt(epsQN_new*epsQN_new + alpha *epsT_new*epsT_new);
		double strsQN_new = kappa * std::sqrt(statev_new(3) * statev_new(3) + multiplier * ((statev_new(16)*statev_new(16))/w2 + (statev_new(17)*statev_new(17))/h2));
		double StrsT_new =  std::sqrt(statev_new(4)*statev_new(4) + statev_new(5)*statev_new(5) + (multiplier * statev_new(15)* statev_new(15))/(w2+h2));		
		statev_new(9) = std::sqrt(strsQN_new*strsQN_new + (StrsT_new*StrsT_new)/alpha);
        statev_new(10) = Wint + statev_old(10);
		statev_new(11) = w;
		statev_new(18) = DE + statev_old(18);

	}
	else {
		mstress.Set(0.0);
		mcouple.Set(0.0);
	}
	//std::cout << "stress: " << mstress << std::endl;
	//std::cout << statev(3) << ' ' << statev(4) << ' ' << statev(5) << ' ' << statev(0) << ' ' << statev(1) << ' ' << statev(2) << std::endl;
	
	}else{
	
	//
	// ELASTIC ANALYSIS
	//
	if (epsQ!=0) {
	    double strsQ=E0*epsQ;
		mstress[0] = strsQ * strain[0] / epsQ;
		mstress[1] = alpha * strsQ * strain[1] / epsQ;
		mstress[2] = alpha * strsQ * strain[2] / epsQ;

		mcouple[0]=multiplier*alpha*strsQ*curvature[0]*(w2+h2)/epsQ;
		mcouple[1]=multiplier*strsQ*curvature[1]*w2/epsQ;
		mcouple[2]=multiplier*strsQ*curvature[2]*h2/epsQ;						

		double Wint = len * area * (((mstress[0] + statev_old(3)) / 2.0) * dstrain[0] + ((mstress[1] + statev_old(4)) / 2.0) * dstrain[1] + ((mstress[2] + statev_old(5)) / 2.0) * dstrain[2] + ((mcouple[0] + statev_old(15)) / 2.0) * dcurvature[0] + ((mcouple[1] + statev_old(16)) / 2.0) * dcurvature[1] + ((mcouple[2] + statev_old(17)) / 2.0) * dcurvature[2]);
		
		double strainN_in = strain[0] - mstress[0] / E0;
		double strainM_in = strain[1] - mstress[1] / (alpha * E0);
		double strainL_in = strain[2] - mstress[2] / (alpha * E0);

		double curvatureN_in = curvature[0] - mcouple[0] / (multiplier * alpha * E0 * (w2 + h2));
		double curvatureM_in = curvature[1] - mcouple[1] / (multiplier * E0 * w2);
		double curvatureL_in = curvature[2] - mcouple[2] / (multiplier * E0 * h2);

		double w_N = len * strainN_in;
		double w_M = len * strainM_in;
		double w_L = len * strainL_in;

		double w_kN = len * std::sqrt(multiplier * (w2 + h2)) * curvatureN_in;
		double w_kM = len * std::sqrt(multiplier * w2) * curvatureM_in;
		double w_kL = len * std::sqrt(multiplier * h2) * curvatureL_in;

		double w = std::sqrt(w_N  * w_N  +  w_M  * w_M  +  w_L  * w_L  +  w_kN * w_kN + w_kM * w_kM + w_kL * w_kL);

		double sgn_strsQ = (mstress[0] > 0.0) ? 1.0 : (mstress[0] < 0.0 ? -1.0 : 0.0);
		
		double dstrainN_in = dstrain[0] - (mstress[0] - statev_old(3)) / E0;
		double dstrainM_in = dstrain[1] - (mstress[1] - statev_old(4)) / (alpha * E0);
		double dstrainL_in = dstrain[2] - (mstress[2] - statev_old(5)) / (alpha * E0);

		double dcurvatureN_in = dcurvature[0] - (mcouple[0] - statev_old(15)) / (multiplier * alpha * E0 * (w2 + h2));
		double dcurvatureM_in = dcurvature[1] - (mcouple[1] - statev_old(16)) / (multiplier * E0 * w2);
		double dcurvatureL_in = dcurvature[2] - (mcouple[2] - statev_old(17)) / (multiplier * E0 * h2);

		double DE = len * area * (0.5 * (mstress[0]  + statev_old(3))  * dstrainN_in +  0.5 * (mstress[1]  + statev_old(4))  * dstrainM_in +  0.5 * (mstress[2]  + statev_old(5))  * dstrainL_in +  0.5 * (mcouple[0] + statev_old(15)) * dcurvatureN_in +   0.5 * (mcouple[1] + statev_old(16)) * dcurvatureM_in +   0.5 * (mcouple[2] + statev_old(17)) * dcurvatureL_in );

		statev_new(0) = strain[0] + eigenstrain[0]; // TODO JBC: it is a bug to only update those inside this code branch.
		statev_new(1) = strain[1] + eigenstrain[1]; //           If you unload to exactly zero (like in the unit tests)
		statev_new(2) = strain[2] + eigenstrain[2]; //           you go to the `else` branch and the state variables do not get updated
		statev_new(3) = mstress[0];                  //           while a loading step did happen!
		statev_new(4) = mstress[1];                  //           I am leaving this as is for now because we will refactor this in depth soon!
		statev_new(5) = mstress[2];                  //           The current fix aims to retrieve the "old" behavior
		statev_new(12) = curvature[0];
        statev_new(13) = curvature[1];
        statev_new(14) = curvature[2];
        statev_new(15) = mcouple[0];
        statev_new(16) = mcouple[1];
        statev_new(17) = mcouple[2];
		double epsQN_new = kappa * std::sqrt(statev_new(0) * statev_new(0) + multiplier*( statev_new(13)*statev_new(13)*w2 + statev_new(14)*statev_new(14)*h2));
		double epsT_new = std::sqrt(statev_new(1) * statev_new(1) + statev_new(2) * statev_new(2)+multiplier*statev_new(12)*statev_new(12)*(w2+h2));
		statev_new(8) =  std::sqrt(epsQN_new*epsQN_new + alpha *epsT_new*epsT_new);
		double strsQN_new = sgn_strsQ * std::sqrt(statev_new(3) * statev_new(3) + multiplier * ((statev_new(16)*statev_new(16))/w2 + (statev_new(17)*statev_new(17))/h2));
		double StrsT_new =  std::sqrt(statev_new(4)*statev_new(4) + statev_new(5)*statev_new(5) + (multiplier * statev_new(15)* statev_new(15))/(w2+h2));
		statev_new(9) = std::sqrt(strsQN_new*strsQN_new + (StrsT_new*StrsT_new)/alpha);
        statev_new(10) = Wint + statev_old(10);
		statev_new(11) = w;
		statev_new(18) = DE + statev_old(18);
	}else{
		mstress.Set(0.0);
		mcouple.Set(0.0);
	}
}
}


void ChWoodMaterialVECT::ComputeStress(const std::shared_ptr<ChViscoelasticity>& visco_params,ChViscoelasticityState& visco_state,ChVector3d& strain_increment, ChVector3d& curvature_increment, ChVector3d& eigenstrain, double &len,const ChVectorDynamic<>& statev_old, ChVectorDynamic<>& statev_new,double& area, double& width, double& height, double& random_field,ChVector3d& mstress, ChVector3d& mcouple){
    if (!m_CreepAnalysis || !visco_params) {
        ComputeStress(strain_increment, curvature_increment, eigenstrain, len, statev_old, statev_new, area, width, height, random_field, mstress, mcouple);
        return;
    }

    const double E0    = this->Get_E0();
    const double alpha = this->Get_alpha();
    const double h2    = height * height;
    const double w2    = width  * width;
    const double multiplier = this->GetCoupleMultiplier() / 3.0;

    const int Nkelv = visco_params->GetNkelv();

    const int idx_DSTRESS = 19;
    const int idx_DTOLD   = idx_DSTRESS + 3;
    const int idx_GAM     = idx_DTOLD + 1;
    const int idx_Z       = idx_GAM + 3 * Nkelv;
    const int idx_CMC0    = idx_Z + 1;
    const int idx_DMC0    = idx_Z + 2;
    const int idx_CMCX    = idx_Z + 3;
    const int idx_AF      = idx_Z + 4;

    const int NSVARS = idx_AF + 1;   // 28 + 3*Nkelv

    if (statev_new.size() < NSVARS) statev_new.resize(NSVARS);

    if ((int)statev_old.size() < NSVARS) {
        mstress[0] = E0 * strain_increment[0];
        mstress[1] = E0 * alpha * strain_increment[1];
        mstress[2] = E0 * alpha * strain_increment[2];
        mcouple[0] = multiplier * alpha * E0 * curvature_increment[0] * (w2 + h2);
        mcouple[1] = multiplier          * E0 * curvature_increment[1] * w2;
        mcouple[2] = multiplier          * E0 * curvature_increment[2] * h2;
        return;
    }

    const ChVector3d stress_prev (statev_old.segment(3, 3));
    const ChVector3d dstress_prev(statev_old.segment(idx_DSTRESS, 3));
    const double dtold = statev_old(idx_DTOLD);

    std::vector<ChVector3d> gamold(Nkelv);
    for (int i = 0; i < Nkelv; ++i) {
        const int idx = idx_GAM + 3 * i;
        gamold[i] = ChVector3d(statev_old(idx), statev_old(idx + 1), statev_old(idx + 2));
    }
    const double Zold  = statev_old(idx_Z);
    const double CMC0  = statev_old(idx_CMC0);
    const double DMC0  = statev_old(idx_DMC0);
    const double CMCX  = statev_old(idx_CMCX);
    const double AFold = statev_old(idx_AF);

    std::vector<ChVector3d> gamnew;
    double CN_visco = 0.0;
    ChVector3d viscous_strain = visco_params->ComputeViscoCBLCON(gamold, dtold, stress_prev, dstress_prev, visco_state, gamnew, CN_visco);

    double Znew = Zold, CN_mech = 0.0;
    ChVector3d mechano_strain = visco_params->ComputeMechanoCBLCON(Zold, stress_prev, dstress_prev, visco_state, Znew, CN_mech);

    double CMCnew = 0.0, DMCnew = 0.0, CMCXnew = 0.0, AFnew = 1.0;
    ChVector3d shrink_strain = visco_params->ComputeShrinkCBLCON(CMC0, CMCX, AFold, DMC0, visco_state, CMCnew, DMCnew, CMCXnew, AFnew);

    ChVector3d thermal_strain = visco_params->ComputeThermalCBLCON(visco_state);

    //ChVector3d eigen_creep = shrink_strain;
     ChVector3d eigen_creep = viscous_strain + AFnew*mechano_strain + shrink_strain + thermal_strain;

    ChVector3d dstrain_mech = strain_increment - eigen_creep;

    const double A0 = visco_params->GetA0();
    const double Q1 = visco_params->GetQ1();
    const double A_total = Q1 + A0 + CN_mech + CN_visco;
    double Em = (A_total > 0.0) ? 1.0 / A_total : E0;

    ChVector3d dstress_new(Em * dstrain_mech[0], Em * alpha * dstrain_mech[1], Em * alpha * dstrain_mech[2]);
    mstress = stress_prev + dstress_new;

    ChVector3d curvature(statev_old.segment(12, 3) + curvature_increment.eigen());
    mcouple[0] = multiplier * alpha * Em * curvature[0] * (w2 + h2);
    mcouple[1] = multiplier          * Em * curvature[1] * w2;
    mcouple[2] = multiplier          * Em * curvature[2] * h2;

    ChVector3d strain(statev_old.segment(0, 3) + dstrain_mech.eigen());

    statev_new(0) = strain[0]; statev_new(1) = strain[1]; statev_new(2) = strain[2];
    statev_new(3) = mstress[0]; statev_new(4) = mstress[1]; statev_new(5) = mstress[2];
    statev_new(12) = curvature[0]; statev_new(13) = curvature[1]; statev_new(14) = curvature[2];
    statev_new(15) = mcouple[0]; statev_new(16) = mcouple[1]; statev_new(17) = mcouple[2];

    statev_new(idx_DSTRESS)     = dstress_new[0];
    statev_new(idx_DSTRESS + 1) = dstress_new[1];
    statev_new(idx_DSTRESS + 2) = dstress_new[2];
    statev_new(idx_DTOLD)       = visco_state.dt;

    for (int i = 0; i < Nkelv; ++i) {
        const int idx = idx_GAM + 3 * i;
        statev_new(idx)     = gamnew[i].x();
        statev_new(idx + 1) = gamnew[i].y();
        statev_new(idx + 2) = gamnew[i].z();
    }
    statev_new(idx_Z)    = Znew;
    statev_new(idx_CMC0) = CMCnew;
    statev_new(idx_DMC0) = DMCnew;
    statev_new(idx_CMCX) = CMCXnew;
    statev_new(idx_AF)   = AFnew;
}

double ChWoodMaterialVECT::FractureBC(ChVector3d& strain, double& random_field, double& len, double& epsQ, double& epsQN,  double& epsT, const ChVectorDynamic<>& statev_old, double eps_max) {
	////
	double E0 = this->Get_E0();
	double alpha = this->Get_alpha();
	double sigmat = this->Get_sigmat();
	double sigmas = this->Get_sigmas();
	double sigmac = this->Get_sigmac();
	double nt = this->Get_nt();   
	double lt = this->Get_lt();
	//if (random_field)
	//	sigmat=sigmat*random_field;
	//double Gt = this->Get_Gt();
	////
	//double epsQ = pow(mstrain[0] * mstrain[0] + alpha * (mstrain[1] * mstrain[1] + mstrain[2] * mstrain[2]), 0.5);
	//double epsT = pow(mstrain[1] * mstrain[1] + mstrain[2] * mstrain[2], 0.5);

	// calculate stress boudary sigma_bt for tension
	double omega, sinw, cosw, cosw2, sigma0;
	double r_st = sigmas / sigmat;
	omega = atan(epsQN / (sqrt(alpha) * epsT));

	sinw = sin(omega);
	cosw = cos(omega);
	cosw2 = cos(omega) * cos(omega);

	if (cosw2 < 10e-16) {   // epsT == 0
		//omega = CH_PI * 0.5;
		sigma0 = sigmat;
	}
	else {
		//omega = atan(mstrain[0] / (pow(alpha, 0.5) * epsT));
		sigma0 = sigmat * (-sin(omega) + sqrt((sin(omega) * sin(omega) + 4 * alpha * cosw2 / r_st / r_st))) / (2 * alpha * cosw2 / r_st / r_st);

	}

	double eps0 = sigma0 / E0;
	//double lt = 2 * E0 * Gt / sigmat / sigmat;
	double Ht = 2 * E0 / (lt / len - 1);
	double H0 = Ht * pow(2 * omega / CH_PI, nt);

	double sigma_bt = sigma0 * exp(-H0 * std::max((eps_max - eps0), 0.0) / sigma0);

	double strs_ela = E0 * (epsQ - statev_old(8)) + statev_old(9); // The elastic prediction increments from old values. TODO: take this out of this function, not its role
	double sigma_fr = std::min(std::max(strs_ela, 0.0), sigma_bt);
	return sigma_fr;
}


double ChWoodMaterialVECT::CompressBC(ChVector3d& strain, double& random_field, double& len, double& epsQ, double& epsT, double& epsQN, const ChVectorDynamic<>& statev_old) {
	//
	double E0 = this->Get_E0();
	double alpha = this->Get_alpha();
	double sigmat = this->Get_sigmat();
	double sigmas = this->Get_sigmas();
	double sigmac = this->Get_sigmac();
	double nt = this->Get_nt();   
	double lt = this->Get_lt();
	//if (random_field)
	//	sigmat=sigmat*random_field;

	// calculate stress boudary sigma_bt for tension
	double omega, sinw, cosw, cosw2, beta2, sigma0;
	
	omega = atan(-epsQN / (pow(alpha, 0.5) * epsT));
	
	beta2 = (sigmas * sigmas) / (sigmac * sigmac);

	sinw = sin(omega);
	cosw = cos(omega);
	//sinw_2 = sin(omega) * sin(omega);
	cosw2 = cosw * cosw;

	//if (omega <= 0) {  
		//omega = CH_PI * 0.5;
		//sigma0 = sigmac;
	//}
	//else {
	sigma0 = sigmac / sqrt((sinw * sinw + (alpha * cosw2) / beta2));
	//}

	double eps0 = sigma0 / E0;
	//double lt = 2 * E0 * Gt / sigmat / sigmat;
	double H0 = 0;

	//double eps_max = epsQ;
	double sigma_bt = sigma0;
	//double sigma_bt = sigma0 * exp(-H0 * std::max((eps_max - eps0), 0.0) / sigma0);
	
	

	double strs_ela = E0 * (epsQ - statev_old(8)) + statev_old(9);  // The elastic prediction increments from old values. TODO: take this out of this function, not its role
	double sigma_fr = std::min(std::max(strs_ela, 0.0), sigma_bt);
	//std::cout<<"epsQ: "<<epsQ<<" epsQ_0: "<<statev(8)<<" sigma_0: "<<statev(9)<<" strs_ela: "<<strs_ela<<"\tsigma_bt: "<<sigma_bt<<"\tsigma_fr: "<<sigma_fr<<std::endl;
	//exit(9);
	return sigma_fr;
}

double ChWoodMaterialVECT::StressBC(ChVector3d& strain, double& random_field, double& len, double& epsQ, double& epsQN,  double& epsT, const ChVectorDynamic<>& statev_old, double eps_max) {
	////
	double E0 = this->Get_E0();
	double alpha = this->Get_alpha();
	double sigmat = this->Get_sigmat();
	double sigmas = this->Get_sigmas();
	double sigmac = this->Get_sigmac();
	double nt = this->Get_nt();   
	double lt = this->Get_lt();
	//if (random_field)
	//	sigmat=sigmat*random_field;
	//double Gt = this->Get_Gt();
	//
	//double epsQ = pow(mstrain[0] * mstrain[0] + alpha * (mstrain[1] * mstrain[1] + mstrain[2] * mstrain[2]), 0.5);
	//double epsT = pow(mstrain[1] * mstrain[1] + mstrain[2] * mstrain[2], 0.5);

	// calculate stress boudary sigma_bt for tension
    double omega, sinw, cosw, cosw2, sinw2, sigma0, a, c, a_b2, H0, eps1;

	a = (sigmat + sigmac)/2.0;
	c = (sigmat - sigmac)/2.0;
	a_b2 = (0.5 * alpha * sigmac * (sigmat + sigmac))/(sigmas*sigmas);

	omega = atan(epsQN / (sqrt(alpha) * epsT));
    if (std::abs(epsQN) < 1e-14) {
        omega = 0.0;
    }
    //std::cout << "omega: " << omega << std::endl;	

	sinw = sin(omega);
	cosw = cos(omega);
	cosw2 = cos(omega) * cos(omega);
	sinw2 = sin(omega) * sin(omega);

	sigma0 = (c*sinw + std::sqrt(c*c*sinw2 -(sinw2 + a_b2*cosw2)*(c*c-a*a))) / (sinw2 + a_b2*cosw2);
	
	double eps  = epsQ;      // current effective strain
	double Ht = 2 * E0 / (lt / len - 1);
	if (omega < 0){
		H0 = 0;
		eps1 = eps;
	}else{
		H0 = Ht * std::pow(2 * omega / CH_PI, nt);
		eps1 = eps_max;
		////eps1 = eps_max;
	};

	double eps0 = sigma0 / E0; // elastic strain limit

	double sigma_bt = sigma0 * exp(-H0 * std::max((eps1 - eps0), 0.0) / sigma0);

	double strs_ela = E0 * (epsQ - statev_old(8)) + statev_old(9); // The elastic prediction increments from old values. TODO: take this out of this function, not its role
	double sigma_fr = std::min(std::max(strs_ela, 0.0), sigma_bt);
	return sigma_fr;
}

}  // end of namespace wood
}  // end of namespace chrono
