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
// Authors: Erol Lale
//			Ke Yu
// =============================================================================
// Material class for LDPM and CSL elements 
//
// A description of the material parameter can be found in: https://doi.org/10.1016/j.cemconcomp.2011.02.011
// =============================================================================

#include "chrono_ldpm/ChMaterialVECT.h"

namespace chrono {
namespace ldpm {

// Construct an isotropic material.

ChMaterialVECT::ChMaterialVECT(double rho,  // material density
                                       double E0,    // Mesoscale Young's modulus
                                       double alpha,   // Mesoscale Poisson ratio                                       
                                       double sigmat, double sigmas, double nt, double lt, 
									   double Ed, double sigmac0, double beta, double Hc0,
									   double Hc1, double kc0, double kc1, double kc2, double kc3,
									   double mu0, double muinf, double sigmaN0, double kt, bool ela_flag, 
									   double ksp, double ksn, double Ef, double t0,
									   double Gd, double beta_fiber, double sigmauf, double krup, double salpha)
    : m_rho(rho), m_E0(E0), m_alpha(alpha), m_sigmat(sigmat), m_sigmas(sigmas), m_nt(nt),
		m_lt(lt), m_Ed(Ed), m_sigmac0(sigmac0), m_beta(beta), m_Hc0(Hc0), m_Hc1(Hc1),
		m_kc0(kc0), m_kc1(kc1), m_kc2(kc2), m_kc3(kc3), m_mu0(mu0), m_muinf(muinf), 
		m_sigmaN0(sigmaN0), m_kt(kt), m_ela(ela_flag), 
        m_ksp(ksp), m_ksn(ksn), m_Ef(Ef), m_t0(t0), m_Gd(Gd), m_beta_fiber(beta_fiber),
        m_sigmauf(sigmauf), m_krup(krup), m_salpha(salpha) {
        
}


ChMaterialVECT::ChMaterialVECT(){
    /*double m_rho=2.4E-9;
    double m_E0=30000;
    double m_alpha=0.2;*/
};

/*
// Copy constructor
ChMaterialVECT::ChMaterialVECT(const ChMaterialVECT& my_material)
	: m_rho(my_material.m_rho), m_E0(my_material.m_E0),
		m_alpha(my_material.m_alpha)
{};
*/

// Destructor
ChMaterialVECT::~ChMaterialVECT()	{};

// statev(0): normal N strain
// statev(1): shear M strain
// statev(2): shear L strain
// statev(3): normal N stress
// statev(4): shear M stress
// statev(5): shear L stress
// statev(6): maximum normal N strain
// statev(7): maximum shear T strain
// statev(8): effective strain
// statev(9): effective stress
// statev(10): internal work
// statev(11): crack opening
// statev(12): normal N eigenstrain
// statev(13): shear M eigenstrain
// statev(14): shear L eigenstrain
// statev(15): dissipated energy
// statev(16): average fiber force
// statev(17): volumetric strain

void ChMaterialVECT::ComputeStress(ChVectorDynamic<>& dmstrain,
                                   double& len,
                                   double& epsV,
                                   ChVectorDynamic<>& statev,
                                   ChVectorDynamic<>& mstress,
                                   double& area,
                                   std::vector<std::vector<double>>& fiber,
                                   ChMatrix33<double>& nmL) {
    mstress.resize(3);    	
	//
	double E0=this->Get_E0();
	double alpha=this->Get_alpha();	
	// 
	ChVectorDynamic<> mstrain(3);
    ChVectorDynamic<> stress_ldpm(3);

	mstrain=statev.segment(0,3)+dmstrain;	
	//
	double epsQ = pow(mstrain(0) * mstrain(0) + alpha * (mstrain(1) * mstrain(1) + mstrain(2) * mstrain(2)), 0.5);
	double epsT = pow(mstrain(1) * mstrain(1) + mstrain(2) * mstrain(2), 0.5);
	
	if (statev(6) < mstrain(0)) {
		statev(6) = mstrain(0);
	}
	if (statev(7) < epsT) {
		statev(7) = epsT;
	}
	
	if (epsQ != 0) {
		if (mstrain(0) > 10e-16) {     // fracture behaivor
			double strsQ = FractureBC(mstrain, len, statev);
			stress_ldpm(0) = strsQ * mstrain(0) / epsQ;
			stress_ldpm(1) = alpha * strsQ * mstrain(1) / epsQ;
			stress_ldpm(2) = alpha * strsQ * mstrain(2) / epsQ;
		}
		else {
			double sigmaN = CompressBC(mstrain, epsV, statev);
			std::pair<double, double> sigmaT = ShearBC(mstrain, statev,sigmaN);
			double sigmaM = sigmaT.first;
			double sigmaL = sigmaT.second;
			//double sigmaTpre = pow(statev(4) * statev(4) + statev(5) * statev(5), 0.5);
			//double sigmaM = sigmaT * statev(4) / sigmaTpre;
			//double sigmaL = sigmaT * statev(5) / sigmaTpre;

			stress_ldpm(0) = sigmaN;
			stress_ldpm(1) = sigmaM;
			stress_ldpm(2) = sigmaL;
		}
		//double Wint = len * area * ((mstress(0)+ statev(3))/2 * dmstrain(0) + (mstress(1)+statev(4)) / 2 * dmstrain(1) +( mstress(2) + statev(5))/2 * dmstrain(2));
		const auto stress_old = statev.segment(3,3);
		
		double Wint = len * area * ((stress_ldpm + stress_old)/2).dot(dmstrain);

		ChVectorN<double, 3> dmstrain_in;
		ChVectorN<double, 3> dstress_in;
		dstress_in(0) = (stress_ldpm(0) - statev(3)) / E0;
		dstress_in(1) = (stress_ldpm(1) - statev(4)) / (E0 * alpha);
		dstress_in(2) = (stress_ldpm(2) - statev(5)) / (E0 * alpha);
		dmstrain_in = dmstrain - dstress_in;
		double DE = len * area * ((stress_ldpm + stress_old)/2).dot(dmstrain_in); // Dissipated energy
		
		//double Wint = len * area * ((stress_ldpm(0)+ statev(3))/2 * dmstrain(0) + (stress_ldpm(1)+statev(4)) / 2 * dmstrain(1) +( stress_ldpm(2) + statev(5))/2 * dmstrain(2));
        
		double w_N, w_M, w_L;

		if (mstrain(0) > 10e-16) {
            w_N = len * (mstrain(0) - stress_ldpm(0) / E0);
            //if (w_N < -10e-16) {
                //std::cout << "w_N" << w_N << "mstrain" << mstrain(0) << "stress" << stress_ldpm(0) / E0 << std::endl;
			//}
            w_M = len * (mstrain(1) - stress_ldpm(1) / (E0 * alpha));
            w_L = len * (mstrain(2) - stress_ldpm(2) / (E0 * alpha));
		}

		//double w = pow(w_N * w_N + w_M * w_M + w_L * w_L, 0.5);

		ChVector3d w(w_N, w_M, w_L);

		// for fluid to solid transition, the effective stress is multiplied by a hydration degree factor s_alpha.
		double salpha = this->Get_salpha();
        //stress_ldpm = stress_ldpm* salpha;
		//mstress = stress_ldpm;
		mstress = stress_ldpm* salpha;
        
		if (!fiber.empty() && w.Length() > 10e-16 && mstrain(0) > 10e-16) {

			if (w_N < 0) {
                w = -w;
            } 

			//std::cout << "start processing fiber "<< std::endl;

            //ChVector3d nn(nmL(0, 0), nmL(0, 1), nmL(0, 2));
            //ChVector3d mm(nmL(1, 0), nmL(1, 1), nmL(1, 2));
            //ChVector3d ll(nmL(2, 0), nmL(2, 1), nmL(2, 2));
            double P_n=0;
            double P_m = 0;
            double P_l = 0;
            double P_fiber = 0;

            for (std::vector<double>& ifiber : fiber) {
                ChVector3d Pf = FiberForce(w, statev(16), ifiber, nmL);
                P_n += Pf.x();
                P_m += Pf.y();
                P_l += Pf.z();
                P_fiber += Pf.Length();

				/*
				if (fiber.size() > 1) {
                    std::cout << "2 fibers in 1 facet" << std::endl;
				}*/

				//std::cout << "nf vector" << ifiber[0] << " " << ifiber[1] << " "<< ifiber[2] << std::endl;
				//std::cout << "Pf: " << Pf.Length() << std::endl;
                //std::cout << "w: " << w.Length() << std::endl;
                //std::cout << "n vector" << nn << std::endl;
				
			}
            
			//std::cout << "stress  " << mstress(0) << " " << mstress(1) << " " << mstress(2) << " " << std::endl;
			//std::cout << "fibers  " << P_n / area << " " << P_m / area << " " << P_l / area << " " << std::endl;


			mstress(0) = mstress(0) + P_n / area;
            mstress(1) = mstress(1) + P_m / area;
            mstress(2) = mstress(2) + P_l / area;

			statev(16) = P_fiber / fiber.size();
			//std::cout << "final  " << mstress(0) << " " << mstress(1) << " " << mstress(2) << " " << std::endl;

			//std::cout << "  " << std::endl;
		}
		
		//std::cout << "stress_ldpm  " << stress_ldpm(0) << " " << stress_ldpm(1) << " " << stress_ldpm(2) << " " << std::endl;
		//std::cout << "stress  " << mstress(0) << " " << mstress(1) << " " << mstress(2) << " " << std::endl;

		statev(0) = mstrain(0);
		statev(1) = mstrain(1);
		statev(2) = mstrain(2);
		statev(3) = stress_ldpm(0);
		statev(4) = stress_ldpm(1);
		statev(5) = stress_ldpm(2);
		statev(8) = pow(statev(0) * statev(0) + alpha * (statev(1) * statev(1) + statev(2) * statev(2)), 0.5);
		statev(9) = pow(statev(3) * statev(3) + (statev(4) * statev(4) + statev(5) * statev(5)) / alpha, 0.5);
		statev(10) = Wint + statev(10);
        statev(11) = w.Length();
		statev(15) = DE + statev(15);
		statev(17) = epsV;
	}
	else {
		mstress << 0.0, 0.0, 0.0;	
		statev(0) = mstrain(0);
		statev(1) = mstrain(1);
		statev(2) = mstrain(2);
		statev(3) = mstress(0);
		statev(4) = mstress(1);
		statev(5) = mstress(2);
		statev(8) = pow(statev(0) * statev(0) + alpha * (statev(1) * statev(1) + statev(2) * statev(2)), 0.5);
		statev(9) = pow(statev(3) * statev(3) + (statev(4) * statev(4) + statev(5) * statev(5)) / alpha, 0.5);	
	}
	
}


void ChMaterialVECT::ComputeStress(ChVectorDynamic<>& dmstrain, ChVectorDynamic<>& eigenstrain, double &len, double &epsV, ChVectorDynamic<>& statev, ChVectorDynamic<>& mstress, double& area) {
    	mstress.resize(3);    	
	//
	double E0=this->Get_E0();
	double alpha=this->Get_alpha();	
	// 
	ChVectorDynamic<> mstrain(3);
	mstrain=statev.segment(0,3)+dmstrain;
	//	
	//
	ChVectorDynamic<> netstrain=mstrain-eigenstrain;
	ChVectorDynamic<> netdmstrain=netstrain-statev.segment(12,3);	
	//std::cout<<"eigenstrain: "<<eigenstrain.x()<<"\t"<<eigenstrain.y()<<"\t"<<eigenstrain.z()<<"\t"	;
	//std::cout<<"mstrain: "<<mstrain.x()<<"\t"<<mstrain.y()<<"\t"<<mstrain.z()<<"\t"	;
	///std::cout<<"dmstrain: "<<dmstrain.x()<<"\t"<<dmstrain.y()<<"\t"<<dmstrain.z()<<"\t"	;
	//std::cout<<"netstrain: "<<netstrain.x()<<"\t"<<netstrain.y()<<"\t"<<netstrain.z()<<"\n";
	//
	double epsQ = pow(netstrain(0) * netstrain(0) + alpha * (netstrain(1) * netstrain(1) + netstrain(2) * netstrain(2)), 0.5);
	double epsT = pow(netstrain(1) * netstrain(1) + netstrain(2) * netstrain(2), 0.5);
	//std::cout<<"epsQ: "<<epsQ<<"\n";
	
	if (statev(6) < netstrain(0)) {
		statev(6) = netstrain(0);
	}
	if (statev(7) < epsT) {
		statev(7) = epsT;
	}
	
	if (epsQ != 0) {
		if (netstrain(0) > 10e-16) {     // fracture behaivor
			double strsQ = FractureBC(netstrain, len, statev);
			mstress(0) = strsQ * netstrain(0) / epsQ;
			mstress(1) = alpha * strsQ * netstrain(1) / epsQ;
			mstress(2) = alpha * strsQ * netstrain(2) / epsQ;
		}
		else {
			double sigmaN = CompressBC(netstrain, netdmstrain, epsV, statev);
			std::pair<double, double> sigmaT = ShearBC(netstrain, netdmstrain, statev,sigmaN);
			double sigmaM = sigmaT.first;
			double sigmaL = sigmaT.second;
			//double sigmaTpre = pow(statev(4) * statev(4) + statev(5) * statev(5), 0.5);
			//double sigmaM = sigmaT * statev(4) / sigmaTpre;
			//double sigmaL = sigmaT * statev(5) / sigmaTpre;
			
			mstress(0) = sigmaN;
			mstress(1) = sigmaM;
			mstress(2) = sigmaL;
			
		}
		
		const auto stress_old = statev.segment(3,3);
		
		double Wint = len * area * ((mstress + stress_old)/2).dot(netdmstrain);

		//double Wint = len * area * ((mstress(0)+ statev(3))/2 * netdmstrain(0) + (mstress(1)+statev(4)) / 2 * netdmstrain(1) +( mstress(2) + statev(5))/2 * netdmstrain(2));
		ChVectorDynamic<> netdmstrain_in(3);
		ChVectorDynamic<> dstress_in(3);
		dstress_in(0) = (mstress(0) - statev(3)) / E0;
		dstress_in(1) = (mstress(1) - statev(4)) / (E0 * alpha);
		dstress_in(2) = (mstress(2) - statev(5)) / (E0 * alpha);
		netdmstrain_in = netdmstrain - dstress_in;
		//double DE = len * area * ((mstress(0)+ statev(3))/2 * dmstrain_in(0) + (mstress(1)+statev(4)) / 2 * dmstrain_in(1) +( mstress(2) + statev(5))/2 * dmstrain_in(2));
		double DE = len * area * ((mstress + stress_old)/2).dot(netdmstrain_in); // Dissipated energy
		
		double w_N = len * (netstrain(0) - mstress(0) / E0);
		double w_M = len * (netstrain(1) - mstress(1) / (E0*alpha));
		double w_L = len * (netstrain(2) - mstress(2) / (E0 * alpha));

		double w = pow(w_N * w_N + w_M * w_M + w_L * w_L, 0.5);

		statev(0) = mstrain(0);
		statev(1) = mstrain(1);
		statev(2) = mstrain(2);
		statev(3) = mstress(0);
		statev(4) = mstress(1);
		statev(5) = mstress(2);
		statev(8) = pow(netstrain(0) * netstrain(0) + alpha * (netstrain(1) * netstrain(1) + netstrain(2) * netstrain(2)), 0.5);
		statev(9) = pow(mstress(0) * mstress(0) + (mstress(1) * mstress(1) + mstress(2) * mstress(2)) / alpha, 0.5);
		statev(10) = Wint + statev(10);
		statev(11) = w;
		statev(12) = netstrain(0);
		statev(13) = netstrain(1);
		statev(14) = netstrain(2);
		statev(15) = DE + statev(15);
		statev(17) = epsV;
	}
	else {
		mstress << 0.0, 0.0, 0.0;	
		statev(0) = mstrain(0);
		statev(1) = mstrain(1);
		statev(2) = mstrain(2);
		statev(3) = mstress(0);
		statev(4) = mstress(1);
		statev(5) = mstress(2);
		statev(8) = pow(statev(0) * statev(0) + alpha * (statev(1) * statev(1) + statev(2) * statev(2)), 0.5);
		statev(9) = pow(statev(3) * statev(3) + (statev(4) * statev(4) + statev(5) * statev(5)) / alpha, 0.5);		
	}
}

double ChMaterialVECT::FractureBC(ChVectorDynamic<>& mstrain, double& len, ChVectorDynamic<>& statev) {
	//
	double E0 = this->Get_E0();
	double alpha = this->Get_alpha();
	double sigmat = this->Get_sigmat();
	double sigmas = this->Get_sigmas();
	double nt = this->Get_nt();   
	double lt = this->Get_lt();
	double kt = this->Get_kt();
	double rs = this->Get_rs();
	//double Gt = this->Get_Gt();
	//
	double epsQ = pow(mstrain(0) * mstrain(0) + alpha * (mstrain(1) * mstrain(1) + mstrain(2) * mstrain(2)), 0.5);
	double epsT = pow(mstrain(1) * mstrain(1) + mstrain(2) * mstrain(2), 0.5);

	// calculate stress boudary sigma_bt for tension
	double omega, sinw, cosw, cosw2, sigma0;
	double r_st = sigmas / sigmat;

	omega = atan(mstrain(0) / (pow(alpha, 0.5) * epsT));
	sinw = sin(omega);
	cosw = cos(omega);
	cosw2 = cos(omega) * cos(omega);

	if (cosw2 < 10e-16) {   // epsT == 0
		//omega = CH_PI * 0.5;
		sigma0 = sigmat;
	}
	else {
		//omega = atan(mstrain(0) / (pow(alpha, 0.5) * epsT));
		sigma0 = sigmat * (-sin(omega) + pow((sin(omega) * sin(omega) + 4 * alpha * cosw2 / r_st / r_st), 0.5)) / (2 * alpha * cosw2 / r_st / r_st);

	}

	double eps0 = sigma0 / E0;
	//double lt = 2 * E0 * Gt / sigmat / sigmat;
	double Ht = 2 * E0 / (lt / len - 1);	
	double Hs = rs* E0;
	//double H0 = Ht * pow(2 * omega / CH_PI, nt);
	double H0 = Hs/alpha + ( Ht- Hs/alpha )* pow(2 * omega / CH_PI, nt);

	double eps_max = pow(statev(6) * statev(6) + alpha * statev(7) * statev(7), 0.5);
	double sigma_bt = sigma0 * exp(-H0 * std::max((eps_max - eps0), 0.0) / sigma0);

	
	double eps_tr = kt * (eps_max - sigma_bt / E0);
	if (eps_tr > epsQ) {
		sigma_bt = 0;
	}
	
	

	double strs_ela = E0 * (epsQ - statev(8)) + statev(9);
	double sigma_fr = std::min(std::max(strs_ela, 0.0), sigma_bt);
	
	bool ela_flag = this->Get_ela();
	if (ela_flag) {
		sigma_fr = strs_ela;
	}
	/*
	if (abs(sigma_fr) < sigmat/1000 && epsQ > eps0) { // to avoid numerical issue when the stress is very small, set it to zero
		sigma_fr = 0.0 * sigma_fr;
		//std::cout << "in the small stress " << std::endl;

	}
	*/
	
	return sigma_fr;
}

double ChMaterialVECT::CompressBC(ChVectorDynamic<>& mstrain, double& epsV, ChVectorDynamic<>& statev) {
	double E0 = this->Get_E0();
	double Ed = this->Get_Ed();
	double sigmac0 = this->Get_sigmac0();
	double beta = this->Get_beta();
	double Hc0 = this->Get_Hc0();
	double Hc1 = this->Get_Hc1();
	double kc0 = this->Get_kc0();
	double kc1 = this->Get_kc1();
	double kc2 = this->Get_kc2();
	double kc3 = this->Get_kc3();

	double epsD = mstrain(0) - epsV;
	double epsDV = epsV + beta * epsD;
	double epsc0 = sigmac0 / E0;
	double epsc1 = epsc0 * kc0;
	double epsv0 = epsc0 * kc3;

	//std::cout << "epsV " << epsV << std::endl;

	double r_DV;
	if (epsV <= 0) {
		r_DV = -abs(epsD) / (epsV - epsv0);
	}
	else {
		r_DV = abs(epsD) / epsv0;
	}

	double Hc = (Hc0 - Hc1) / (1 + kc2 * std::max(r_DV - kc1, 0.0)) + Hc1;
	double sigmac1 = sigmac0 + (epsc1 - epsc0) * Hc;

	double sigma_bc;
	if (-epsDV <= 0) {
		sigma_bc = sigmac0;
	}
	else if (-epsDV > 0 && -epsDV <= epsc1) {
		sigma_bc = sigmac0 + std::max((-epsDV - epsc0), 0.0) * Hc;
		//std::cout << "case2 "  << std::endl;
	}
	else {
		sigma_bc = sigmac1 * exp((-epsDV - epsc1) * Hc / sigmac1);
		//std::cout << "case3 " << std::endl;
	}

	//std::cout << "epsDV " << epsDV << std::endl;
	//std::cout << "Hc " << Hc << std::endl;
	//std::cout << "epsc1 " << epsc1 << std::endl;

	double ENc;
	if (-statev(3) < sigmac0) {
		ENc = E0;
	}
	else {
		ENc = Ed;
	}

	double sigma_ela = statev(3) + ENc * (mstrain(0) - statev(0));
	double sigma_com = std::min(std::max(sigma_ela, -sigma_bc), 0.0);
	
	bool ela_flag = this->Get_ela();
	if (ela_flag) {
		sigma_com = sigma_ela;
	}
	return sigma_com;
}

std::pair<double, double> ChMaterialVECT::ShearBC(ChVectorDynamic<>& mstrain, ChVectorDynamic<>& statev, double& sigmaN) {
	double E0 = this->Get_E0();
	double alpha = this->Get_alpha();
	double mu0 = this->Get_mu0();
	double muinf = this->Get_muinf();
	double sigmas = this->Get_sigmas();
	double sigmaN0 = this->Get_sigmaN0();

	//double sigmabs = sigmas + (mu0 - muinf) * sigmaN0 - muinf * statev(3) - (mu0 - muinf) * sigmaN0 * exp(statev(3) / sigmaN0);
	double sigmabs = sigmas + (mu0 - muinf) * sigmaN0 - muinf * sigmaN - (mu0 - muinf) * sigmaN0 * exp(sigmaN / sigmaN0);

	double ET = alpha * E0;
	double sigmaM_ela = statev(4) + (mstrain(1) - statev(1)) * ET;
	double sigmaL_ela = statev(5) + (mstrain(2) - statev(2)) * ET;
	double sigmaT_ela = pow(sigmaM_ela * sigmaM_ela + sigmaL_ela * sigmaL_ela, 0.5);

	double sigmaT = std::min(std::max(sigmaT_ela, 0.0), sigmabs);

	bool ela_flag = this->Get_ela();
	if (ela_flag) {
		sigmaT = sigmaT_ela;
	}

	double sigmaM, sigmaL;

	if (sigmaT_ela != 0) {
		sigmaM = sigmaT * sigmaM_ela / sigmaT_ela;
		sigmaL = sigmaT * sigmaL_ela / sigmaT_ela;
	}
	else {
		sigmaM = 0;
		sigmaL = 0;
	}
	
	return std::make_pair(sigmaM, sigmaL);
}



double ChMaterialVECT::CompressBC(ChVectorDynamic<>& mstrain, ChVectorDynamic<>& dmstrain, double& epsV, ChVectorDynamic<>& statev) {
	double E0 = this->Get_E0();
	double Ed = this->Get_Ed();
	double sigmac0 = this->Get_sigmac0();
	double beta = this->Get_beta();
	double Hc0 = this->Get_Hc0();
	double Hc1 = this->Get_Hc1();
	double kc0 = this->Get_kc0();
	double kc1 = this->Get_kc1();
	double kc2 = this->Get_kc2();
	double kc3 = this->Get_kc3();

	double epsD = mstrain(0) - epsV;
	double epsDV = epsV + beta * epsD;
	double epsc0 = sigmac0 / E0;
	double epsc1 = epsc0 * kc0;
	double epsv0 = epsc0 * kc3;

	//std::cout << "epsV " << epsV << std::endl;

	double r_DV;
	if (epsV <= 0) {
		r_DV = -abs(epsD) / (epsV - epsv0);
	}
	else {
		r_DV = abs(epsD) / epsv0;
	}

	double Hc = (Hc0 - Hc1) / (1 + kc2 * std::max(r_DV - kc1, 0.0)) + Hc1;
	double sigmac1 = sigmac0 + (epsc1 - epsc0) * Hc;

	double sigma_bc;
	if (-epsDV <= 0) {
		sigma_bc = sigmac0;
	}
	else if (-epsDV > 0 && -epsDV <= epsc1) {
		sigma_bc = sigmac0 + std::max((-epsDV - epsc0), 0.0) * Hc;
		//std::cout << "case2 "  << std::endl;
	}
	else {
		sigma_bc = sigmac1 * exp((-epsDV - epsc1) * Hc / sigmac1);
		//std::cout << "case3 " << std::endl;
	}

	//std::cout << "epsDV " << epsDV << std::endl;
	//std::cout << "Hc " << Hc << std::endl;
	//std::cout << "epsc1 " << epsc1 << std::endl;

	double ENc;
	if (-statev(3) < sigmac0) {
		ENc = E0;
	}
	else {
		ENc = Ed;
	}

	double sigma_ela = statev(3) + ENc * (dmstrain(0));
	double sigma_com = std::min(std::max(sigma_ela, -sigma_bc), 0.0);
	
	bool ela_flag = this->Get_ela();
	if (ela_flag) {
		sigma_com = sigma_ela;
	}

	return sigma_com;
}

std::pair<double, double> ChMaterialVECT::ShearBC(ChVectorDynamic<>& mstrain, ChVectorDynamic<>& dmstrain, ChVectorDynamic<>& statev, double& sigmaN) {
	double E0 = this->Get_E0();
	double alpha = this->Get_alpha();
	double mu0 = this->Get_mu0();
	double muinf = this->Get_muinf();
	double sigmas = this->Get_sigmas();
	double sigmaN0 = this->Get_sigmaN0();

	//double sigmabs = sigmas + (mu0 - muinf) * sigmaN0 - muinf * statev(3) - (mu0 - muinf) * sigmaN0 * exp(statev(3) / sigmaN0);
	double sigmabs = sigmas + (mu0 - muinf) * sigmaN0 - muinf * sigmaN - (mu0 - muinf) * sigmaN0 * exp(sigmaN / sigmaN0);

	double ET = alpha * E0;
	double sigmaM_ela = statev(4) + (dmstrain(1) ) * ET;
	double sigmaL_ela = statev(5) + (dmstrain(2) ) * ET;
	double sigmaT_ela = pow(sigmaM_ela * sigmaM_ela + sigmaL_ela * sigmaL_ela, 0.5);

	double sigmaT = std::min(std::max(sigmaT_ela, 0.0), sigmabs);

	bool ela_flag = this->Get_ela();
	if (ela_flag) {
		sigmaT = sigmaT_ela;
	}

	double sigmaM, sigmaL;

	if (sigmaT_ela != 0) {
		sigmaM = sigmaT * sigmaM_ela / sigmaT_ela;
		sigmaL = sigmaT * sigmaL_ela / sigmaT_ela;
	}
	else {
		sigmaM = 0;
		sigmaL = 0;
	}
	
	return std::make_pair(sigmaM, sigmaL);
}

ChVector3d ChMaterialVECT::FiberForce(ChVector3d w, double Pf0, std::vector<double>& ifiber, ChMatrix33<double>& nmL) {

	ChVector3d nf(ifiber[0], ifiber[1], ifiber[2]);
    double Ls = ifiber[3];
    double Ll = ifiber[4];
    double df = ifiber[5];
    double vmaxs = ifiber[7];
    double vmaxl = ifiber[8];
    double pmaxs = ifiber[9];
    double pmaxl = ifiber[10];
    //ChVector3d n(nmL(0, 0), nmL(0, 1), nmL(0, 2));

	double ksp = this->Get_ksp();
    double ksn = this->Get_ksn();
    double Ef = this->Get_Ef();
    double t0 = this->Get_t0();
    double Gd = this->Get_Gd();
    double beta_fiber = this->Get_beta_fiber();



	ChVector3d nn(nmL(0, 0), nmL(0, 1), nmL(0, 2));

    ChVector3d nf_local = (nmL * nf).GetNormalized();

	if (nf_local.x() < 0) {
        nf_local = nf_local * (-1);
	}

	auto [Pf,Ps,vs, vl, kf] =
        ChFiber::FRP(w, Pf0, nf_local, nn, ksp, ksn, df, Ef, t0, Gd, Ls, Ll, beta_fiber, vmaxs, vmaxl, pmaxs, pmaxl);

	if (vs > vmaxs) {
        ifiber[7] = vs;
	}
    if (vl > vmaxl) {
        ifiber[8] = vl;
    }
    if (Ps > pmaxs) {
        ifiber[9] = Ps;
	}
    if (Ps > pmaxl) {
        ifiber[10] = Ps;
    }
	
	ifiber[11] = kf;
	
	double sigmauf = this->Get_sigmauf();
    double krup = this->Get_krup();

    double phi_prime = acos(nf_local.Dot(Pf.GetNormalized()));
    double Pf_u = sigmauf * exp(-krup * phi_prime) * CH_PI * df * df/4;

	if (Pf_u <= 0) {
        std::cout << "Pf_u < 0 "<< std::endl;
	}

    if (Pf.Length() > Pf_u) {
        Pf = Pf_u * Pf.GetNormalized();
	}
    //std::cout << "nf vector" << nf << std::endl;
    //std::cout << "nf' vector" << nf_local << std::endl;
	//std::cout << "kf: " << kf << std::endl;


	return Pf;
}


}  // end of namespace ldpm
}  // end of namespace chrono
