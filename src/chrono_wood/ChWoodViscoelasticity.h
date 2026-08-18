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
// Authors: Wisdom Akpan, Erol Lale
// =============================================================================

#ifndef CHWOOD_VISCOELASTICITY_H
#define CHWOOD_VISCOELASTICITY_H

#include <cmath>
#include <memory>
#include <vector>

#include "chrono_wood/ChWoodApi.h"
#include "chrono/core/ChVector3.h"
#include "chrono/core/ChMatrix33.h"
#include "chrono/core/ChMatrix.h"

namespace chrono {
namespace wood {

// 
// Environmental input interfaces  (shared, time-function based)
// 

class ChWoodApi ChRHFunction {
public:
    virtual ~ChRHFunction() = default;
    virtual double GetRH(double t) const = 0;
};

class ChWoodApi ChTempFunction {
public:
    virtual ~ChTempFunction() = default;
    virtual double GetTemp(double t) const = 0;
};

// 

struct ChWoodApi ChViscoelasticityState {

    double RH    = 0.0;  
    double Temp  = 0.0;  
    double DRH   = 0.0;  
    double DTemp = 0.0;   

    double dt    = 0.0;   
    double CTime = 0.0;   

    double RH_prev        = 0.0;
    double Temp_prev      = 0.0;
    double last_env_time  = -1.0e30;   
    bool   env_initialized = false;

    // When true: RH/Temp/DRH/DTemp already set by coupling; skip global 
    bool use_coupled_env = false;

    /// Update environment from global time-functions
    void UpdateFromGlobalFunctions(
    double t, double step,
    const std::shared_ptr<ChRHFunction>& rh_fn, double initial_RH,
    const std::shared_ptr<ChTempFunction>& temp_fn, double initial_Temp)
{
    if (use_coupled_env) return;

    const double rh   = rh_fn   ? rh_fn->GetRH(t)     : initial_RH;
    const double temp = temp_fn ? temp_fn->GetTemp(t) : initial_Temp;

    if (!env_initialized) {
        RH_prev = rh;
        Temp_prev = temp;
        RH = rh;
        Temp = temp;
        DRH = 0.0;
        DTemp = 0.0;
        last_env_time = t;
        env_initialized = true;
        return;
    }

    if (std::abs(t - last_env_time) < 1.0e-15) return;

    DRH   = (step > 0.0) ? rh   - RH_prev   : 0.0;
    DTemp = (step > 0.0) ? temp - Temp_prev : 0.0;

    RH   = RH_prev;
    Temp = Temp_prev;

    RH_prev = rh;
    Temp_prev = temp;
    last_env_time = t;
}

    /// Set environment directly from the flow–mechanics coupling layer.
    void SetFromCoupling(double T, double h, double dT, double dh) {
        Temp = T;  RH  = h;
        DTemp = dT; DRH = dh;
        use_coupled_env = true;
    }

    void ClearCoupledEnv() {
        use_coupled_env = false;
    }
};


class ChWoodApi ChViscoelasticity {
public:

    ChViscoelasticity();

    ChViscoelasticity(double mEl,   double mQ1,
                      double mCFA,  double mCFB,
                      double mtau0, double mA0,
                      double mQ3,   double mKaps, double mTaus,
                      double mSMF,  double mKappaS, double mKappaT);

    virtual ~ChViscoelasticity();

    void ComputeKelvinProps();

    void ComputeGinv_Beam(double nu_rat);

    void ComputeGinvCBLCON(double alpha, double& width, double& height);


    ChVector3d ComputeViscoBeam(
        const std::vector<ChVector3d>& gamold, double dtold,
        const ChVector3d& stress_prev, const ChVector3d& dstress,
        const ChViscoelasticityState& env_state,
        std::vector<ChVector3d>& gamnew_out, double& CN_visco_out);

    ChVector3d ComputeViscoCBLCON(
        const std::vector<ChVector3d>& gamold, double dtold,
        const ChVector3d& stress_prev, const ChVector3d& dstress,
        const ChViscoelasticityState& env_state,
        std::vector<ChVector3d>& gamnew_out, double& CN_visco_out);


    ChVector3d ComputeMechanoBeam(
        double Zold,
        const ChVector3d& stress_prev, const ChVector3d& dstress,
        const ChViscoelasticityState& env_state,
        double& Znew_out, double& CN_mech_out);

    ChVector3d ComputeMechanoCBLCON(
        double Zold,
        const ChVector3d& stress_prev, const ChVector3d& dstress,
        const ChViscoelasticityState& env_state,
        double& Znew_out, double& CN_mech_out);


    ChVector3d ComputeShrinkBeam(
        double CMC0_old, double CMCX_old, double AF_old, double DMC0_old,
        const ChViscoelasticityState& env_state,
        double& CMCnew_out, double& DMCnew_out,
        double& CMCXnew_out, double& AFnew_out);

    ChVector3d ComputeShrinkCBLCON(
        double CMC0_old, double CMCX_old, double AF_old, double DMC0_old,
        const ChViscoelasticityState& env_state,
        double& CMCnew_out, double& DMCnew_out,
        double& CMCXnew_out, double& AFnew_out);

    ChVector3d ComputeThermalBeam(const ChViscoelasticityState& env_state);
    ChVector3d ComputeThermalCBLCON(const ChViscoelasticityState& env_state);

    void   SetEl(double v)     { El     = v; }
    double GetEl()       const { return El; }

    void   SetQ1(double v)     { Q1     = v; }
    double GetQ1()       const { return Q1; }

    void   SetA0(double v)     { A0     = v; }
    double GetA0()       const { return A0; }

    void   SetCFA(double v)    { CFA    = v; }
    double GetCFA()      const { return CFA; }

    void   SetCFB(double v)    { CFB    = v; }
    double GetCFB()      const { return CFB; }

    void   SetTau0(double v)   { tau0   = v; }
    double GetTau0()     const { return tau0; }

    void   SetNkelv(int v)     { Nkelv  = v; }
    int    GetNkelv()    const { return Nkelv; }

    void   SetQ3(double v)     { Q3     = v; }
    double GetQ3()       const { return Q3; }

    void   SetKaps(double v)   { Kaps   = v; }
    double GetKaps()     const { return Kaps; }

    void   SetTaus(double v)   { Taus   = v; }
    double GetTaus()     const { return Taus; }

    void   SetSMF(double v)    { SMF    = v; }
    double GetSMF()      const { return SMF; }

    void   SetKappaS(double v) { KappaS = v; }
    double GetKappaS()   const { return KappaS; }

    void   SetKappaT(double v) { KappaT = v; }
    double GetKappaT()   const { return KappaT; }

    void SetRHFunction(const std::shared_ptr<ChRHFunction>& f)     { RH_function   = f; }
    void SetInitialRH(double rh)                                    { initial_RH    = rh; }
    const std::shared_ptr<ChRHFunction>& GetRHFunction()     const { return RH_function; }
    double GetInitialRH()                                    const { return initial_RH; }

    void SetTempFunction(const std::shared_ptr<ChTempFunction>& f) { Temp_function = f; }
    void SetInitialTemp(double T)                                   { initial_Temp  = T; }
    const std::shared_ptr<ChTempFunction>& GetTempFunction() const { return Temp_function; }
    double GetInitialTemp()                                  const { return initial_Temp; }

    const ChMatrix33<double>& GetGinv_Beam()   const { return G_inv_beam; }
    const ChMatrix33<double>& GetGinvCBLCON()  const { return G_inv_CBLCON; }

private:

    ChVector3d ComputeViscoImpl(
        const std::vector<ChVector3d>& gamold, double dtold,
        const ChVector3d& stress_prev, const ChVector3d& dstress,
        double dt,
        std::vector<ChVector3d>& gamnew_out, double& CN_visco_out);

    ChVector3d ComputeMechanoImpl(
        double Zold,
        const ChVector3d& stress_prev,
        const ChVector3d& dstress,
        const ChViscoelasticityState& env_state,
        double& Znew_out, double& CN_mech_out);

    ChVector3d ComputeShrinkImpl(
        double CMC0_old, double CMCX_old, double AF_old, double DMC0_old,
        const ChViscoelasticityState& env_state,
        double& CMCnew_out, double& DMCnew_out,
        double& CMCXnew_out, double& AFnew_out);

    ChVector3d ComputeThermalImpl(const ChViscoelasticityState& env_state);
    // 
    double El     = 0.0;
    double Q1     = 0.0;
    double CFA    = 0.0;
    double CFB    = 0.0;
    double tau0   = 0.0;
    double A0     = 0.0;
    int    Nkelv  = 0;
    double Taus   = 0.0;
    double Q3     = 0.0;
    double Kaps   = 0.0;
    double KappaS = 0.0;
    double KappaT = 0.0;
    double SMF    = 0.0;

    ChMatrixDynamic<>  KelvinChain;  
    ChMatrix33<double> G_inv_beam;
    ChMatrix33<double> G_inv_CBLCON;

    std::shared_ptr<ChRHFunction>   RH_function;
    std::shared_ptr<ChTempFunction> Temp_function;
    double initial_RH   = 0.0;
    double initial_Temp = 0.0;
};

} // end namespace wood
} // end namespace chrono

#endif  // CHWOOD_VISCOELASTICITY_H
