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
//          Ke Yu
// =============================================================================
// Material class for LDPM and CSL elements 
//
// A description of the material parameter can be found in: https://doi.org/10.1016/j.cemconcomp.2011.02.011
// =============================================================================

#ifndef CHMATERIALVECT_H
#define CHMATERIALVECT_H

#include "chrono_ldpm/ChLdpmApi.h"
#include "chrono_ldpm/ChFiber.h"
#include "chrono/core/ChMatrix33.h"
#include <vector>
#include <string>
//#include "chrono/fea/ChElementBeam.h"

namespace chrono {
namespace ldpm {

/// @addtogroup fea_elements
/// @{

/// Definition of materials to be used for CSL beams and LDPM tets utilizing the lattice discrete particle model.
class ChLdpmApi ChMaterialVECT {
  public:
    /// Construct an isotropic elastic material.
    ChMaterialVECT(double rho,  		///< material density
                       double E0,    	///< Mesoscale Young's modulus
                       double alpha,   	///< Mesoscale Poisson like parameter
                       double sigmat, double sigmas, double nt, double lt, 
					   double Ed, double sigmac0, double beta, double Hc0,
					   double Hc1, double kc0, double kc1, double kc2, double kc3,
					   double mu0, double muinf, double sigmaN0, double kt, bool ela_flag, 
                       double ksp, double ksn, double Ef, double t0, double Gd,
                       double beta_fiber, double sigmauf, double krup, double salpha
    );
    
    ChMaterialVECT();

    /*
    //copy constructor:
	ChMaterialVECT(const ChMaterialVECT& my_material);
	*/

    // Destructor declared:
    ~ChMaterialVECT();

    /// Return the material density.
     double Get_density() const { return m_rho; }
     
     void Set_density(double rho) { m_rho=rho; }
	
    /// Return the Macroscale Modulus of Elasticity.
    double Get_E0() const { return m_E0; }
	  
    void Set_E0(double E0) { m_E0=E0; }
    
    /// Return the Macroscale Poisson Ratio.
    double Get_alpha() const { return m_alpha; }
    
    void Set_alpha(double alpha) { m_alpha=alpha; }
    
    /// Return the tensile strength.
    double Get_sigmat() const { return m_sigmat; }
    void Set_sigmat(double sigmat) { m_sigmat = sigmat; }

    /// Return the shear strength.
    double Get_sigmas() const { return m_sigmas; }
    void Set_sigmas(double sigmas) { m_sigmas = sigmas; }

    /// Return the softening exponent.
    double Get_nt() const { return m_nt; }
    void Set_nt(double nt) { m_nt = nt; }

    /// Return the fracture energy.
    // double Get_Gt() const { return m_Gt; }
    // void Set_Gt(double Gt) { m_Gt = Gt; }

    /// Return the tensile characteristic length.
    double Get_lt() const { return m_lt; }
    void Set_lt(double lt) { m_lt = lt; }

    /// Return the densified normal modulus.
    double Get_Ed() const { return m_Ed; }
    void Set_Ed(double Ed) { m_Ed = Ed; }

    /// Return the compressive yield strength.
    double Get_sigmac0() const { return m_sigmac0; }
    void Set_sigmac0(double sigmac0) { m_sigmac0 = sigmac0; }

    /// Return the volunmetric deviatoric coupling.
    double Get_beta() const { return m_beta; }
    void Set_beta(double beta) { m_beta = beta; }

    /// Return the initial hardening modulus.
    double Get_Hc0() const { return m_Hc0; }
    void Set_Hc0(double Hc0) { m_Hc0 = Hc0; }

    /// Return the final hardening modulus.
    double Get_Hc1() const { return m_Hc1; }
    void Set_Hc1(double Hc1) { m_Hc1 = Hc1; }

    /// Return the transitional strain ratio.
    double Get_kc0() const { return m_kc0; }
    void Set_kc0(double kc0) { m_kc0 = kc0; }

    /// Return the deviatoric strain threshold ratio.
    double Get_kc1() const { return m_kc1; }
    void Set_kc1(double kc1) { m_kc1 = kc1; }

    /// Return the deviatoric damage parameter.
    double Get_kc2() const { return m_kc2; }
    void Set_kc2(double kc2) { m_kc2 = kc2; }

    /// Return the material parameter, 0.1.
    double Get_kc3() const { return m_kc3; }
    void Set_kc3(double kc3) { m_kc3 = kc3; }

    /// Return the initial friction.
    double Get_mu0() const { return m_mu0; }
    void Set_mu0(double mu0) { m_mu0 = mu0; }

    /// Return the asymptotic friction
    double Get_muinf() const { return m_muinf; }
    void Set_muinf(double muinf) { m_muinf = muinf; }

    /// Return the transitional stress.
    double Get_sigmaN0() const { return m_sigmaN0; }
    void Set_sigmaN0(double sigmaN0) { m_sigmaN0 = sigmaN0; }

    /// Return the tensile unloading.
    double Get_kt() const { return m_kt; }
    void Set_kt(double kt) { m_kt = kt; }
	
	/// Return the Shear softening modulus ratio.
    double Get_rs() const { return m_rs; }
    void Set_rs(double rs) { m_rs = rs; }
    
    /// Return the tensile unloading.
    double Get_ela() const { return m_ela; }
    void Set_ela(bool ela_flag) { m_ela = ela_flag; }

    /// Set and Get RayleighDampingK coefficient.
    double GetRayleighDampingK() const { return RayleighDampingK; }
    void SetRayleighDampingK(double myRayleighDampingK) { RayleighDampingK = myRayleighDampingK; }
    
    /// Set and Get RayleighDampingM coefficient.
    double GetRayleighDampingM() const { return RayleighDampingM; }
    void SetRayleighDampingM(double myRayleighDampingM) { RayleighDampingM = myRayleighDampingM; }

    /// Parameters for fiber
    double Get_ksp() const { return m_ksp; }
    void Set_ksp(double ksp) { m_ksp = ksp; }

    double Get_ksn() const { return m_ksn; }
    void Set_ksn(double ksn) { m_ksn = ksn; }

    double Get_Ef() const { return m_Ef; }
    void Set_Ef(double Ef) { m_Ef = Ef; }

    double Get_t0() const { return m_t0; }
    void Set_t0(double t0) { m_t0 = t0; }

    double Get_Gd() const { return m_Gd; }
    void Set_Gd(double Gd) { m_Gd = Gd; }

    double Get_beta_fiber() const { return m_beta_fiber; }
    void Set_beta_fiber(double beta_fiber) { m_beta_fiber = beta_fiber; }

    double Get_krup() const { return m_krup; }
    void Set_krup(double krup) { m_krup = krup; }

    double Get_sigmauf() const { return m_sigmauf; }
    void Set_sigmauf(double sigmauf) { m_sigmauf = sigmauf; }

    double Get_salpha() const { return m_salpha; }
    void Set_salpha(double salpha) { m_salpha = salpha; }


    /// Compute stresses from given strains and state variables.
    void ComputeStress(ChVectorDynamic<>& mstrain, double &len, double& epsV, ChVectorDynamic<>& statev,ChVectorDynamic<>& mstress, double& area, std::vector<std::vector<double>>& fiber, ChMatrix33<double>& nmL);
    void ComputeStress(ChVectorDynamic<>& mstrain, ChVectorDynamic<>& eigenstrain, double &len, double& epsV, ChVectorDynamic<>& statev,ChVectorDynamic<>& mstress, double& area);

    double FractureBC(ChVectorDynamic<>& mstrain, double& len, ChVectorDynamic<>& statev);

    double CompressBC(ChVectorDynamic<>& mstrain, double& epsV, ChVectorDynamic<>& statev);

    std::pair<double, double> ShearBC(ChVectorDynamic<>& mstrain, ChVectorDynamic<>& statev, double& sigmaN);
    
    double CompressBC(ChVectorDynamic<>& mstrain, ChVectorDynamic<>& dmstrain, double& epsV, ChVectorDynamic<>& statev);
    std::pair<double, double> ShearBC(ChVectorDynamic<>& mstrain, ChVectorDynamic<>& dmstrain, ChVectorDynamic<>& statev, double& sigmaN);

    ChVector3d FiberForce(ChVector3d w, double Pf0, std::vector<double>& ifiber, ChMatrix33<double>& nmL);
    
  private:
    
    double m_rho;               ///< density
    double m_E0;    			///< mesoscale modulus of elasticity
    double m_alpha; 			///< mesoscale Poisson like parameter
    double m_sigmat;            ///< tensile strength
    double m_sigmas;            ///< shear strength
    double m_nt;                ///< softening exponent
    //double m_Gt;                ///< fracture energy
    double m_lt;                ///< tensile characteristic length
    double m_Ed;                ///< densified normal modulus
    double m_sigmac0;           ///< compressive yield strength
    double m_beta;              ///< volunmetric deviatoric coupling
    double m_Hc0;               ///< initial hardening modulus
    double m_Hc1;               ///< final hardening modulus
    double m_kc0;               ///< transitional strain ratio
    double m_kc1;               ///< deviatoric strain threshold ratio
    double m_kc2;               ///< deviatoric damage parameter
    double m_kc3;               ///< material parameter, 0.1
    double m_mu0;               ///< initial friction
    double m_muinf;             ///< asymptotic friction
    double m_sigmaN0;           ///< transitional stress.
    double m_kt = 0;                ///< tensile unloading
	  double m_rs = 0;              ///< Shear softening modulus ratio
    bool m_ela = false;                 ///< elastic analysis or not
    double RayleighDampingK=0;
    double RayleighDampingM=0;

    double m_ksp;
    double m_ksn;
    double m_Ef;
    double m_t0;
    double m_Gd;
    double m_beta_fiber;
    double m_sigmauf;
    double m_krup;

    double m_salpha = 1;
  //public:
  //EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

/// @} fea_elements

}  // end of namespace ldpm
}  // end of namespace chrono

#endif
