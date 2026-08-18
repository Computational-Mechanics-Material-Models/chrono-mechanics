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

#include "chrono_wood/ChWoodViscoelasticity.h"

#include <cstdlib>
#include <iostream>

namespace chrono {
namespace wood {

// =============================================================================
// Constructors / destructor
// =============================================================================

ChViscoelasticity::ChViscoelasticity() {}

ChViscoelasticity::ChViscoelasticity(
    double mEl,   double mQ1,
    double mCFA,  double mCFB,
    double mtau0, double mA0,
    double mQ3,   double mKaps, double mTaus,
    double mSMF,  double mKappaS, double mKappaT)
    : El(mEl), Q1(mQ1),
      CFA(mCFA), CFB(mCFB), tau0(mtau0), A0(mA0),
      Q3(mQ3), Kaps(mKaps), Taus(mTaus),
      SMF(mSMF), KappaS(mKappaS), KappaT(mKappaT) {}

ChViscoelasticity::~ChViscoelasticity() {}


void ChViscoelasticity::ComputeKelvinProps() {
    KelvinChain.resize(Nkelv, 2);

    // Retardation times: logarithmic spacing starting from tau0
    KelvinChain(0, 0) = tau0;
    for (int i = 1; i < Nkelv; ++i)
        KelvinChain(i, 0) = 10.0 * KelvinChain(i - 1, 0);

    // Retardation moduli from power-law creep compliance
    const double factor = (CFA / 2.0) * CFB * (CFB - 1.0) * (CFB - 2.0);
    for (int i = 0; i < Nkelv; ++i) {
        const double tau  = KelvinChain(i, 0);
        const double Ltau = factor * std::pow(3.0 * tau, CFB);
        KelvinChain(i, 1) = Ltau * std::log(10.0);
    }
}

void ChViscoelasticity::ComputeGinv_Beam(double nu_rat) {
    G_inv_beam.setZero();
    G_inv_beam(0, 0) = 1.0;
    G_inv_beam(1, 1) = 2.0 * (1.0 + nu_rat);
    G_inv_beam(2, 2) = 2.0 * (1.0 + nu_rat);
}

void ChViscoelasticity::ComputeGinvCBLCON(double alpha, double& width, double& height) {
    const double h2 = height * height;
    const double w2 = width  * width;
    G_inv_CBLCON.setZero();
    G_inv_CBLCON(0, 0) = (alpha / 3.0) * (h2 + w2);
    G_inv_CBLCON(1, 1) = w2 / 3.0;
    G_inv_CBLCON(2, 2) = h2 / 3.0;
}

// =============================================================================
// PRIVATE SHARED IMPLEMENTATIONS
// =============================================================================

ChVector3d ChViscoelasticity::ComputeViscoImpl(
    const std::vector<ChVector3d>& gamold,
    double dtold,
    const ChVector3d& stress_prev,
    const ChVector3d& dstress,
    double dt,
    std::vector<ChVector3d>& gamnew_out,
    double& CN_visco_out)
{
    ChVector3d strain(0, 0, 0);
    gamnew_out.assign(Nkelv, ChVector3d(0, 0, 0));
    CN_visco_out = 0.0;

   
    if (dt <= 0.0 || dtold <= 0.0) {
        gamnew_out = gamold;
        return strain;
    }

    double CN = 0.0;
    for (int i = 0; i < Nkelv; ++i) {
        const double tau = KelvinChain(i, 0);
        const double Amu = KelvinChain(i, 1);
        const double lmd = (1.0 - std::exp(-dt / tau)) / (dt / tau);
        CN += (1.0 - lmd) * Amu;
    }

    for (int i = 0; i < Nkelv; ++i) {
        const double tau     = KelvinChain(i, 0);
        const double Amu     = KelvinChain(i, 1);
        const double exp_old = std::exp(-dtold / tau);
        const double lmd_old = (1.0 - exp_old) / (dtold / tau);

        gamnew_out[i] = gamold[i] * exp_old
                      + (1.0 - exp_old)        * Amu * (stress_prev - dstress)
                      + Amu * (1.0 - lmd_old)  * dstress;
    }

    for (int i = 0; i < Nkelv; ++i) {
        const double tau = KelvinChain(i, 0);
        const double Amu = KelvinChain(i, 1);
        strain += (1.0 - std::exp(-dt / tau)) * (Amu * stress_prev - gamnew_out[i]);
    }

    CN_visco_out = CN;
    return strain;
}

// Mechano-sorptive creep 
ChVector3d ChViscoelasticity::ComputeMechanoImpl(
    double Zold,
    const ChVector3d& stress_prev,
    const ChVector3d& /*dstress*/,
    const ChViscoelasticityState& env,
    double& Znew_out,
    double& CN_mech_out)
{
    constexpr double Tol   = 1.0e-10;
    constexpr double TRef  = 273.0;
    constexpr int    MIter = 1000;

    ChVector3d strain(0, 0, 0);
    CN_mech_out = 0.0;
    Znew_out    = Zold;

    // No mechano-sorptive driving force at the very start
    if (env.CTime <= 0.0 || env.dt <= 0.0)
        return strain;

    const double Temp_m = env.Temp + 0.5 * env.DTemp;
    const double RH_m   = env.RH   + 0.5 * env.DRH;

    if (RH_m <= 0.0)
        return strain;

    const double CEnv = env.DTemp * std::log(RH_m)
                      + (Temp_m + TRef) * env.DRH / RH_m;

    double DZ = 0.0;
    double Z  = Zold;

    for (int i = 0; i < MIter; ++i) {
        const double Zm  = Zold + 0.5 * DZ;
        const double DZ0 = DZ;

        DZ = -Taus * env.dt * Zm * Zm + (Kaps / Taus) * std::abs(CEnv);

        Z = Zold + DZ;
        const double err = std::abs(DZ0 - DZ);

        if (err <= Tol) {
            if (Z < 0.0) {
                std::cerr << "[ChViscoelasticity] MPS less than zero.\n";
                std::exit(EXIT_FAILURE);
            }
            const double Zm_final = Zold + 0.5 * DZ;
            strain      = Q3 * Z * stress_prev;
            CN_mech_out = 0.5 * Q3 * Zm_final * env.dt;
            Znew_out    = Z;
            return strain;
        }
    }

    std::cerr << "[ChViscoelasticity] MPS did not converge.\n";
    std::exit(EXIT_FAILURE);
}

// Shrinkage / swelling  
ChVector3d ChViscoelasticity::ComputeShrinkImpl(
    double CMC0_old, double CMCX_old, double AF_old, double DMC0_old,
    const ChViscoelasticityState& env,
    double& CMCnew_out, double& DMCnew_out,
    double& CMCXnew_out, double& AFnew_out)
{
    const double RHM = env.RH + 0.5 * env.DRH;

    // sorption isotherm constants
    constexpr double XM = 6.350;
    constexpr double C  = 7.818;
    constexpr double K  = 0.785;

    const double CMC = (XM * C * K * RHM)
                     / ((1.0 - K * RHM) * (1.0 + (C - 1.0) * K * RHM))
                     / 100.0;

    // Initialise on first call (CTime == 0 or moisture state not yet set)
    if (env.CTime <= 0.0 || (CMC0_old == 0.0 && CMCX_old == 0.0)) {
        CMCnew_out  = CMC;
        DMCnew_out  = 0.0;
        CMCXnew_out = CMC;
        AFnew_out   = 1.0;
        return ChVector3d(0, 0, 0);
    }

    const double DMC = CMC - CMC0_old;

    double AFnew   = 1.0;
    double CMCXnew = CMCX_old;
    if (CMC > CMCX_old) {
        AFnew   = 1.0 + SMF * (CMC - CMCX_old);
        CMCXnew = CMC;
    }

    ChVector3d strain(0, 0, 0);
    strain.x() = (DMC / 0.3) * KappaS;

    CMCnew_out  = CMC;
    DMCnew_out  = DMC;
    CMCXnew_out = CMCXnew;
    AFnew_out   = AFnew;

    return strain;
}
// Thermal expansion
ChVector3d ChViscoelasticity::ComputeThermalImpl(const ChViscoelasticityState& env) {
    ChVector3d strain(0, 0, 0);
    strain.x() = KappaT * env.DTemp;
    return strain;
}

// =============================================================================
// PUBLIC API  —  Beam variants  
// =============================================================================

ChVector3d ChViscoelasticity::ComputeViscoBeam(
    const std::vector<ChVector3d>& gamold, double dtold,
    const ChVector3d& stress_prev, const ChVector3d& dstress,
    const ChViscoelasticityState& env_state,
    std::vector<ChVector3d>& gamnew_out, double& CN_visco_out)
{
    return ComputeViscoImpl(gamold, dtold, stress_prev, dstress,
                            env_state.dt, gamnew_out, CN_visco_out);
}

ChVector3d ChViscoelasticity::ComputeMechanoBeam(
    double Zold,
    const ChVector3d& stress_prev, const ChVector3d& dstress,
    const ChViscoelasticityState& env_state,
    double& Znew_out, double& CN_mech_out)
{
    return ComputeMechanoImpl(Zold, stress_prev, dstress, env_state,
                              Znew_out, CN_mech_out);
}

ChVector3d ChViscoelasticity::ComputeShrinkBeam(
    double CMC0_old, double CMCX_old, double AF_old, double DMC0_old,
    const ChViscoelasticityState& env_state,
    double& CMCnew_out, double& DMCnew_out,
    double& CMCXnew_out, double& AFnew_out)
{
    return ComputeShrinkImpl(CMC0_old, CMCX_old, AF_old, DMC0_old, env_state,
                             CMCnew_out, DMCnew_out, CMCXnew_out, AFnew_out);
}

ChVector3d ChViscoelasticity::ComputeThermalBeam(const ChViscoelasticityState& env_state) {
    return ComputeThermalImpl(env_state);
}

// =============================================================================
// PUBLIC API 
// =============================================================================

ChVector3d ChViscoelasticity::ComputeViscoCBLCON(
    const std::vector<ChVector3d>& gamold, double dtold,
    const ChVector3d& stress_prev, const ChVector3d& dstress,
    const ChViscoelasticityState& env_state,
    std::vector<ChVector3d>& gamnew_out, double& CN_visco_out)
{
    return ComputeViscoImpl(gamold, dtold, stress_prev, dstress,
                            env_state.dt, gamnew_out, CN_visco_out);
}

ChVector3d ChViscoelasticity::ComputeMechanoCBLCON(
    double Zold,
    const ChVector3d& stress_prev, const ChVector3d& dstress,
    const ChViscoelasticityState& env_state,
    double& Znew_out, double& CN_mech_out)
{
    return ComputeMechanoImpl(Zold, stress_prev, dstress, env_state,
                              Znew_out, CN_mech_out);
}

ChVector3d ChViscoelasticity::ComputeShrinkCBLCON(
    double CMC0_old, double CMCX_old, double AF_old, double DMC0_old,
    const ChViscoelasticityState& env_state,
    double& CMCnew_out, double& DMCnew_out,
    double& CMCXnew_out, double& AFnew_out)
{
    return ComputeShrinkImpl(CMC0_old, CMCX_old, AF_old, DMC0_old, env_state,
                             CMCnew_out, DMCnew_out, CMCXnew_out, AFnew_out);
}

ChVector3d ChViscoelasticity::ComputeThermalCBLCON(const ChViscoelasticityState& env_state) {
    return ComputeThermalImpl(env_state);
}

} // end namespace wood
} // end namespace chrono
