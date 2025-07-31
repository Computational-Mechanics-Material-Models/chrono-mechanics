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

#include "chrono/core/ChMatrix.h"
#include "chrono/core/ChTypes.h"
#include "chrono/core/ChVector3.h"
#include "gtest/gtest.h"

#include "chrono_ldpm/ChMaterialVECT.h"

#include <vector>
#include <fstream>

using namespace chrono;
using namespace ldpm;

// TODO: Write a Fixture to avoid repeating the code for setting up the test

TEST(LdpmMaterialVECTTest, stress_no_eigenstrain){
    bool debug_mode = true; // true prints results to file. false runs the tests
    std::ofstream res;

    // Parameters from Appendix A. of our Overleaf
    // TODO: update reference once this is published
    
    // Material parameters (SI units)
    const double MPa_to_Pa = 1e6;
    const double mm_to_m = 1e-3;

    double rho = 2380; // Density
    double E0 = 60273 * MPa_to_Pa; // Normal modulus
    double alpha = 0.25; // Alpha
    double sigmat = 3.44 * MPa_to_Pa; // Tensile strength
    double lt = 500 * mm_to_m; // Tensile charachteristic length
    double Gt = lt * sigmat * sigmat / (2 * E0); // Fracture energy
    double rt = 2.6; // Shear strength ratio
    double nt = 0.4; // Softening exponent
    double sigmac0 = 150 * MPa_to_Pa; // Compressive yield strength
    double Hc0 = 0.40 * E0; // Initial hardening modulus
    double kc0 = 4; // Transitional strain ratio
    double kc1 = 1; // Deviatoric strain threshold ratio
    double kc2 = 5; // Deviatoric damage parameter
    double kc3 = 0.1; // Volumetric strain parameter
    double mu0 = 0.4; // Initial friction
    double muinf = 0; // Asymptotic friction
    double sigmaN0 = 600 * MPa_to_Pa; // Transitional stress
    double Ed = 1 * E0; // Densification ratio
	double beta = 0; // Volumetric deviatoric coupling
    double kt = 0.0; // Tensile unloading
    double ks = 0.0; // Shear unloading
    double kc = 0.0; // Compressive unloading
    double rs = 0.0; // Shear softening modulus ratio
    double Hc1 = 0.1 * E0; // Final hardening modulus ratio
    
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
    

    // LDPM facet parameters
    double length = 12 * mm_to_m; // MUST BE SMALLER THAN `lt`
    double facet_area = (length * 0.2) * (length * 0.2);

    // Funtions for analytical calculation of the stress in facets according to Appendix A
    // our overleaf project https://www.overleaf.com/project/65e64a729ad76d75d59ec7e2 (replace with Ref once available)

    // -- Fracturing behavior (Appendix A.1.) -- //

    // Coupling variable
    auto omega = [alpha](double eN, double eT) {
        double w;
        if (eT <= 0.0) { // Effective tangent strain is zero
            if (eN > 0.0)      w = 0.5 * M_PI;
            else if (eN < 0.0) w = - 0.5 * M_PI;
            else               w = 0.0; // TO DISCUSS: not sure the value for this 0/0 undefined case matters
        }
        else                   w = std::atan(eN / (eT * std::sqrt(alpha)));
        return w;
    };

    // Elastic increment + Loading - unloading hysteresis in tension
        // A simple elastic increment is not directly implemented because when the
        // reloading strain is crossed upon loading, the stress must be zero up to
        // the reloading strain, and incrementally elastic above the reloading strain.
        // We implement this distinction, which is correct for arbitrarily large steps.
        // This does not enforce bounds (The return algorithm is responsible for that)
    auto t_elastic_tension = [E0](double t_old, double e_old, double etr_old, double e, double etr) {
        double t = 0.0;
        if (e > etr) { // Stress may be non-zero only if strain is above reloading strain
            if (e_old < etr_old) // Reloading from below the reloading strain
                t = E0 * (e - etr_old); // Elastic increment from reloading strain at zero stress
            else // Reloading or unloading from above the reloading strain
                t = t_old + E0 * (e - e_old); // Elastic increment from previous state
        } // etr increases when emax increases: impossible to *load* from (e_old > etr_old) and get to (e < etr)
        return t;
    };

    // Strength limit for effective traction - Equation (A.3)
    auto sigma_0 = [sigmas, sigmat, alpha](double w) {
        double rst = sigmas / sigmat;
        // The standard equation: s0 = sigmat * (-std::sin(w) + std::sqrt(std::sin(w)*std::sin(w) +
        //                                       4 * a * std::cos(w)*std::cos(w) / (rst * rst)) ) /
        //                                      (2 * a * std::cos(w)*std::cos(w) / (rst * rst) );
        // is indeterminate for w -> pi/2 so we rationalize the expression
        // by multiplying both the numerator and denominator by the numerator's conjugate
        return sigmat * 2 / (std::sin(w) + std::sqrt(std::sin(w)*std::sin(w) + 4 * alpha * std::cos(w)*std::cos(w) / (rst * rst)) );
    };

    // Softening modulus - Equation (A.4)
    auto H_0 = [E0, alpha, nt, lt, rs](double w, double l) {
        double Ht =  2 * E0 / (lt / l - 1.0);
        double Hs =  rs * E0 ;
        return Hs / alpha + (Ht - Hs / alpha) * std::pow(2 * w / M_PI, nt);
    };

    // Tensile stress boundary - Equation (A.2)
    auto sigma_bt = [&sigma_0, &H_0, E0](double w, double e, double emax, double l) {
        double sig0 = sigma_0(w);
        double H0 = H_0(w, l);
        double eps0 = sig0 / E0;
        return sig0 * std::exp(-H0 * std::max(emax - eps0, 0.0) / sig0);
    };



    // -- Compressive behavior (Appendix A.2.) -- //

    // Loading - unloading stiffness modulus - Equation (A.5)
        // Equation (A.5) is not directly implemented because during
        // a large increment of strain, the loading-unloading stiffness might change
        // if the stress crosses the compressive yield strength sigmac0.
        // Instead of using the instantaneous stiffness at the current value of stress,
        // we implement the entire bilinear increment, which is always correct.
        // This does not enforce bounds (The return algorithm is responsible for that)
    auto tN_elastic_compression = [E0, Ed, sigmac0](double tN_old, double eN_incr) {
        double ENc = (-tN_old <= sigmac0) ? E0 : Ed;
        double eN_incr_0 = -(sigmac0 + tN_old) / ENc;
        double Ebranch = (eN_incr < eN_incr_0) ? Ed : E0;
        return -sigmac0 + (eN_incr - eN_incr_0) * Ebranch;
    };


    // Initial hardening modulus - Equation (A.7)
    auto H_c = [Hc0, Hc1, kc1, kc2](double rdv) {
        return (Hc0 - Hc1) / (1.0 + kc2 * std::max(rdv - kc1, 0.0)) + Hc1;
    };

    // Deviatoric to volumetric strain ratio - Equation (A.8)
    auto r_DV = [Hc0, Hc1, kc3, sigmac0, E0](double ev, double ed) {
        double ev0 = kc3 * sigmac0 / E0;
        return (ev <= 0) ? - std::abs(ed) / (ev - ev0) : std::abs(ed) / ev0 ;
    };

    // Compressive stress boundary - Equation (A.6)
    auto sigma_bc = [&H_c, &r_DV, beta, sigmac0, E0, kc0](double ev, double ed) {
        double sigbc;
        double edv = ev + beta * ed;

        double ec0 = sigmac0 / E0;
        double ec1 = kc0 * ec0;

        double rdv = r_DV(ev, ed);
        double Hc = H_c(rdv);
        double sigmac1 = sigmac0 + (ec1 - ec0) * Hc;

        if (-edv <= 0.0) {
            sigbc = sigmac0;
        } else if (-edv <= ec1) {
            sigbc = sigmac0 + std::max(-edv - ec0, 0.0) * Hc;
        } else {
            sigbc = sigmac1 * std::exp((-edv - ec1) * Hc / sigmac1);
        }
        return sigbc;
    };


    // -- Frictional behavior in compression (Appendix A.3.) -- //

    // Shear strength - Equation (A.9)
    auto sigma_bs = [sigmas, mu0, muinf, sigmaN0](double tN) {
        return sigmas + (mu0 - muinf) * sigmaN0 - muinf * tN - (mu0 - muinf) * sigmaN0 * std::exp(tN / sigmaN0);
    };


    // --------- Loading paths --------- //
    int npaths = 9;
    int nsteps = 1000;
    // Strain
    std::vector<std::vector<double>> epsN_paths(npaths, std::vector<double>(nsteps));
    std::vector<std::vector<double>> epsM_paths(npaths, std::vector<double>(nsteps));
    std::vector<std::vector<double>> epsL_paths(npaths, std::vector<double>(nsteps));
    std::vector<std::vector<double>> epsV_paths(npaths, std::vector<double>(nsteps));

    double eps_tension(10 * sigmat / E0);
    double eps_compression(-10 * sigmac0 / E0);
    double eps_shearM(0.7 * (eps_tension - eps_compression));
    double eps_shearL(-0.4 * (eps_tension - eps_compression));
    double eps_compactionV(eps_compression);

    // TODO: discuss what are the most relevant loading path to exercise the material in a demanding manner
    for (int t = 0 ; t < nsteps ; t++) {
        double ramp = double(t)/(nsteps - 1);
        // Path 0 : Monotonic axial tension, no volumetric strain
        epsN_paths[0][t] = eps_tension * ramp;
        epsM_paths[0][t] = 0.0;
        epsL_paths[0][t] = 0.0;
        epsV_paths[0][t] = 0.0;

        // Path 1 : cyclic axial tension with increasing max and min strain, no volumetric strain
        double cycle = 0.5*(1 - std::cos(double(t)/nsteps * (2*M_PI) * 4));
        epsN_paths[1][t] = eps_tension * cycle * ramp + ramp*0.75*eps_tension;
        epsM_paths[1][t] = 0.0;
        epsL_paths[1][t] = 0.0;
        epsV_paths[1][t] = 0.0;

        // Path 2 : Monotonic axial compression, no volumetric strain
        epsN_paths[2][t] = eps_compression * ramp;
        epsM_paths[2][t] = 0.0;
        epsL_paths[2][t] = 0.0;
        epsV_paths[2][t] = 0.0;

        // Path 3 : Monotonic axial compression, compressive volumetric strain
        epsN_paths[3][t] = eps_compression * ramp;
        epsM_paths[3][t] = 0.0;
        epsL_paths[3][t] = 0.0;
        epsV_paths[3][t] = eps_compactionV * ramp;

        // Path 4 : cyclic axial compression with increasing max and min strain, compressive volumetric strain
        epsN_paths[4][t] = eps_compression * cycle * ramp + ramp*0.1*eps_compression;
        epsM_paths[4][t] = 0.0;
        epsL_paths[4][t] = 0.0;
        epsV_paths[4][t] = eps_compactionV * cycle * ramp;

        // Path 5 : Monotonic axial compression + shear, no volumetric strain
        epsN_paths[5][t] = eps_compression * ramp;
        epsM_paths[5][t] = eps_shearM * ramp;
        epsL_paths[5][t] = eps_shearL * ramp;
        epsV_paths[5][t] = 0.0;

        // Path 6 : Monotonic axial compression + shear, compressive volumetric strain
        epsN_paths[6][t] = eps_compression * ramp;
        epsM_paths[6][t] = eps_shearM * ramp;
        epsL_paths[6][t] = eps_shearL * ramp;
        epsV_paths[6][t] = eps_compactionV * ramp;

        // Path 7 : Monotonic axial compression + cyclic shear, no volumetric strain
        epsN_paths[7][t] = eps_compression * ramp;
        epsM_paths[7][t] = eps_shearM * cycle;
        epsL_paths[7][t] = eps_shearL * cycle;
        epsV_paths[7][t] = 0.0;

        // Path 8 : instant axial tension followed cyclic shear, no volumetric strain
        // Question: by definition of the effective stress, cycles of shear will create normal stress at constant normal strain. Is this ok ?
        epsN_paths[8][t] = eps_tension;
        epsM_paths[8][t] = eps_shearM * cycle;
        epsL_paths[8][t] = eps_shearL * cycle;
        epsV_paths[8][t] = 0.0;
    }

    // Comparison of analytical calculations and Chrono implementation for all loading paths
    for (int path = 0 ; path < npaths ; path++) {
        // Initialize history and incremental variables
        double sig = 0.0;
        double epsNmax = 0.0;
        double epsTmax = 0.0;
        double epsold = 0.0;
        double epstrold = 0.0;
        double sigN = 0.0;
        double sigM = 0.0;
        double sigL = 0.0;
        double sigNold = 0.0;
        double sigMold = 0.0;
        double sigLold = 0.0;
        double Wint = 0.0;

        // TODO: ChMaterialVECT still uses ChVectorDynamic<> for vector of known size 12, to replace by ChVectorN<double, 12>
        ChVectorDynamic<> statev; // ChVectorN<double, 12> statev;
        statev.setZero(12);

        if (debug_mode) { // Output to file to plot
            res.open ("results_path"+std::to_string(path)+".csv");
            res << "# Analytical solution from unit test LdpmMaterialVECTTest - stress_no_eigenstrain.\n";
            res << "# E0="<<E0<<", alpha="<<alpha<<", sigmat="<<sigmat<<", lt="<<lt<<", Gt="<<Gt<<", rt"<<rt
                <<", nt="<<nt<<", sigmac0="<<sigmac0<<", Hc0="<<Hc0<<", kc0="<<kc0<<", kc1="<<kc2<<", kc3="<<kc3
                <<", mu0="<<mu0<<", muinf="<<muinf<<", sigmaN0="<<sigmaN0<<", Ed="<<Ed<<", beta="<<beta<<", ke="<<kt
                <<", ks="<<ks<<", kc="<<kc<<", rs="<<rs<<", Hc1="<<Hc1<<", length="<<length<<"\n";

            res << "step, epsN, epsM, epsL, eps, epsV, epsNmax, epsTmax, sigN, sigM, sigL, sig, Wint, crack\n";
            res << "0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0\n";
        }
        for (int t = 1 ; t < nsteps ; t++) {
            double epsN = epsN_paths[path][t];
            double epsM = epsM_paths[path][t];
            double epsL = epsL_paths[path][t];
            double epsV = epsV_paths[path][t];
            double epsD = epsN - epsV;
            double epsNincr = epsN - epsN_paths[path][t-1];
            double epsMincr = epsM - epsM_paths[path][t-1];
            double epsLincr = epsL - epsL_paths[path][t-1];


            // -------- Analytical calculation -------- //

            double epsT = std::sqrt(epsM * epsM + epsL * epsL);
            double eps = std::sqrt(epsN * epsN + alpha * epsT * epsT);   
            epsNmax = std::max(epsNmax, epsN);
            epsTmax = std::max(epsTmax, epsT);

            if (epsN >= 0.0) { // Fracturing behavior

                // Elastic prediction
                double w = omega(epsN, epsT);
                double epsmax = std::sqrt(epsNmax * epsNmax + alpha * epsTmax * epsTmax);
                double sigbt = sigma_bt(w, eps, epsmax, length);
                double epstr = std::max(kt * (epsmax - sigbt / E0), 0.0);

                sig = t_elastic_tension(sig, epsold, epstrold, eps, epstr); // Manages unloading - reloading hysteresis

                // Correction (return)
                if (sig < 0.0) sig = 0.0;
                if (sig > sigbt) sig = sigbt;

                // Stress components
                if (eps > 0.0) { // TODO: is there anything safer ? MIN_DBL ?
                    sigN = sig * epsN / eps;
                    sigM = sig * alpha * epsM / eps;
                    sigL = sig * alpha * epsL / eps;
                } else {
                    sigN = 0.0;
                    sigM = 0.0;
                    sigL = 0.0;
                }

                // Update for next increment
                epsold = eps;
                epstrold = epstr;
            } else { // Compressive + Frictional behavior

                // Elastic prediction
                sigN = tN_elastic_compression(sigNold, epsNincr);
                sigM += E0 * alpha * epsMincr;
                sigL += E0 * alpha * epsLincr;

                // Normal correction (return)
                double sigbc = sigma_bc(epsV, epsD);
                if (sigN > 0.0) sigN = 0.0;
                if (sigN < -sigbc) sigN = -sigbc;

                // Tangential correction (radial return)
                double sigbs = sigma_bs(sigN);
                double sigT = std::sqrt(sigM * sigM + sigL * sigL);
                double phi = sigT - sigbs;
                if (phi > 0.0) {
                    sigM *= sigbs / sigT;
                    sigL *= sigbs / sigT;
                }
                sig = std::sqrt(sigN * sigN + (sigM * sigM + sigL * sigL) / alpha);
            }

            // Total Energy of mechanical deformation (trapezoidal integration)
            Wint += length * facet_area * (epsNincr * 0.5 * (sigNold + sigN) + epsMincr * 0.5 * (sigMold + sigM) + epsLincr * 0.5 * (sigLold + sigL));

            // Crack Opening
            double crackN = length * std::max(epsN - sigN / E0, 0.0); // TODO: Normal Crack only if positive ?
            double crackM = length * (epsM - sigM / (alpha * E0)); // Is this tangential sliding really a "crack" ?
            double crackL = length * (epsL - sigL / (alpha * E0)); // Is this tangential sliding really a "crack" ?
            double crack = std::sqrt(crackN * crackN + crackM * crackM + crackL * crackL); // TODO: should this be zero entirely if normal crack is zero? i.e., shear in compression is not a crack?

            // Update stress for next increment
            sigNold = sigN;
            sigMold = sigM;
            sigLold = sigL;

            // -------- Chrono implementation -------- //

            // TODO: chMaterialVECT still uses ChVectorDynamic<> for vector of known size 3, to replace by chVector3d
            ChVectorDynamic<> strain_increment(3); // ChVector3d strain_increment(epsNincr, epsMincr, epsMincr);
            strain_increment(0) = epsNincr; strain_increment(1) = epsMincr; strain_increment(2) = epsMincr;
            ChVectorDynamic<> stress;
            stress.setZero(3);
            my_mat->ComputeStress(strain_increment, length, epsV, statev, stress, facet_area);

            // ---------------------------------------- //
            // -------- Comparison and Testing -------- //
            // ---------------------------------------- //

            if (!debug_mode) {
                double tol = 1e-6;

                // The state variables are updatead inside my_mat->ComputeStress()
                // so they contain the current values.
                // Maybe we should test for the old values before calling my_mat->ComputeStress() ?
                ASSERT_NEAR(epsN, statev(0), tol);
                ASSERT_NEAR(epsM, statev(1), tol);
                ASSERT_NEAR(epsL, statev(2), tol);
                ASSERT_NEAR(sigN, statev(3), tol);
                ASSERT_NEAR(sigM, statev(4), tol);
                ASSERT_NEAR(sigL, statev(5), tol);
                ASSERT_NEAR(epsNmax, statev(6), tol);
                ASSERT_NEAR(epsTmax, statev(7), tol);
                ASSERT_NEAR(eps, statev(8), tol);
                ASSERT_NEAR(sig, statev(9), tol);
                ASSERT_NEAR(Wint, statev(10), tol);
                ASSERT_NEAR(crack, statev(11), tol);

                ASSERT_NEAR(sigN, stress[0], tol);
                ASSERT_NEAR(sigM, stress[1], tol);
                ASSERT_NEAR(sigL, stress[2], tol);
            }

            if (debug_mode) {// Output
                res<<t<<", "<<epsN<<", "<<epsM<<", "<<epsL<<", "<<eps<<", "<<epsV<<", "<<epsNmax<<", "<<epsTmax<<", "<<sigN<<", "<<sigM<<", "<<sigL<<", "<<sig<<", "<<Wint<<", "<<crack<<std::endl;
            }
        }

        res.close();
    }
}
