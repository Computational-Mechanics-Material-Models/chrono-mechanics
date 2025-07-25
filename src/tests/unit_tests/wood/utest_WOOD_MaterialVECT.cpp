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

#include "chrono/core/ChTypes.h"
#include "chrono/core/ChVector3.h"
#include "gtest/gtest.h"

#include "chrono_wood/ChWoodMaterialVECT.h"

#include <vector>
#include <fstream>

using namespace chrono;
using namespace wood;

// TODO: Write a Fixture to avoid repeating the code for setting up the test

TEST(WoodMaterialVECTTest, stress_no_eigenstrain){
    bool debug_mode = true; // true prints results to file. false runs the tests
    std::ofstream res;
    
    // Material parameters
    double rho = 5e-8;
    double E0 = 8000;
    double alpha = 0.2373;
    double sigmat = 30.0;
    double sigmac = 120.0;
    double sigmas = 78.0;
    double nt = 0.2;   
    double lt = 5.0;
    double rs = 0.0;
    double couple_multiplier = 0.761;
    double material_is_elastic = false;

    auto my_mat = chrono_types::make_shared<ChWoodMaterialVECT>();
    my_mat->Set_density(rho);
    my_mat->Set_E0(E0);
    my_mat->Set_alpha(alpha);
    my_mat->Set_sigmat(sigmat);
    my_mat->Set_sigmac(sigmac);
    my_mat->Set_sigmas(sigmas);
    my_mat->Set_nt(nt);
    my_mat->Set_lt(lt);
    my_mat->Set_rs(rs);
    my_mat->SetCoupleMultiplier(couple_multiplier);
    my_mat->SetElasticAnalysisFlag(material_is_elastic);

    // CBL connector parameters
    double length = 3.0;
    double facet_width = 3.165;
    double facet_height = 12.4864;
    double facet_area = facet_height * facet_width;
    ChVectorN<double, 18> statev;
    statev.setZero();
    ChVector3d strain(4.186, 1.16, -86.48);
    ChVector3d curvature(0.123, -0.864, 0.793);

    // CBL stress calculation
    ChVector3d stress;
    ChVector3d couple;
    double random_field = 1.0; // NO RANDOM FIELD
    my_mat->ComputeStress(strain,curvature, length, statev, facet_area, facet_width, facet_height, random_field, stress, couple);

    double hM = 0.5 * facet_height;
    double hL = 0.5 * facet_width;
    double rNsq = (hL * hL + hM * hM) / 3.0;
    double rMsq = (hL * hL) / 3.0;
    double rLsq = (hM * hM) / 3.0;
    double betaN(couple_multiplier), betaM(couple_multiplier), betaL(couple_multiplier);

    // Funtions for analytical calculation of the stress in connectors according to
    // our overleaf project https://www.overleaf.com/project/681aa3904a3faf7f1699077b (replace with Ref once available)

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

    // Strength decay function
    auto epsilon_1 = [](double e, double emax, double w) {
        double e1 = e;
        if (w > 0.0) e1 = emax;
        return e1;
    };

    // Elastic stress limit
    auto sigma_0 = [sigmas, sigmat, sigmac, alpha](double w) {
        double s0;
        double rst = sigmas / sigmat;
        double beta = (sigmas * sigmas) / (sigmac * sigmac);
        if (w > 0.0) {
            // The standard equation: s0 = sigmat * (-std::sin(w) + std::sqrt(std::sin(w)*std::sin(w) +
            //                                       4 * a * std::cos(w)*std::cos(w) / (rst * rst)) ) /
            //                                      (2 * a * std::cos(w)*std::cos(w) / (rst * rst) );
            // is indeterminate for w -> pi/2 so we rationalize the expression
            // by multiplying both the numerator and denominator by the numerator's conjugate
            s0 = sigmat * 2 / (std::sin(w) + std::sqrt(std::sin(w)*std::sin(w) + 4 * alpha * std::cos(w)*std::cos(w) / (rst * rst)) );
        }
        else         s0 = sigmac / std::sqrt(std::sin(w)*std::sin(w) + alpha * std::cos(w)*std::cos(w) / beta);
        return s0;
    };

    // Softening modulus
    // TODO: what is rs ?!?! Not given i the CBL code right now. Is it the same as rst ?
    auto H_0 = [E0, alpha, nt, lt, rs](double w, double l) {
        double h0 = 0.0;
        if (w > 0.0) {
            double Ht =  2 * E0 / (lt / l - 1.0);
            double Hs =  rs * E0 ;
            h0 = Hs / alpha + (Ht - Hs / alpha) * std::pow(2*w / M_PI, nt);
        }
        return h0;
    };

    // Limiting stress
    auto sigma_bs = [&omega, &epsilon_1, &sigma_0, &H_0, E0](double w, double e, double emax, double l) {
        double sig0 = sigma_0(w);
        double H0 = H_0(w, l);
        double eps1 = epsilon_1(e, emax, w);
        double eps0 = sig0 / E0;
        return sig0 * std::exp(-H0 * std::max(emax - eps0, 0.0) / sig0);
    };


    // --------- Loading paths --------- //
    int npaths = 10;
    int nsteps = 1000;
    // Strain
    std::vector<std::vector<double>> epsN_paths(npaths, std::vector<double>(nsteps));
    std::vector<std::vector<double>> epsM_paths(npaths, std::vector<double>(nsteps));
    std::vector<std::vector<double>> epsL_paths(npaths, std::vector<double>(nsteps));
    // Curvature
    std::vector<std::vector<double>> chiN_paths(npaths, std::vector<double>(nsteps));
    std::vector<std::vector<double>> chiM_paths(npaths, std::vector<double>(nsteps));
    std::vector<std::vector<double>> chiL_paths(npaths, std::vector<double>(nsteps));

    double eps_tension(2 * sigmat / E0);
    double eps_compression(-2 * sigmac / E0);
    double eps_shearM(0.1 * (eps_tension - eps_compression));
    double eps_shearL(-0.2 * (eps_tension - eps_compression));

    double chi_torsion(2 * sigmas / (alpha * E0) / std::sqrt(rNsq));
    double chi_bendingM(2 * sigmat / E0 / hL);
    double chi_bendingL(-2 * sigmat / E0 / hM);

    for (int t = 0 ; t < nsteps ; t++) { // TODO: it is a bit clunky to add new cases by manually giving the index. Just use push_back()
        double ramp = double(t)/(nsteps - 1);
        // Path 0 : pure monotonic tension
        epsN_paths[0][t] = eps_tension * ramp;
        epsM_paths[0][t] = 0.0;
        epsL_paths[0][t] = 0.0;
        chiN_paths[0][t] = 0.0;
        chiM_paths[0][t] = 0.0;
        chiL_paths[0][t] = 0.0;

        // Path 1 : pure monotonic compression
        epsN_paths[1][t] = eps_compression * ramp;
        epsM_paths[1][t] = 0.0;
        epsL_paths[1][t] = 0.0;
        chiN_paths[1][t] = 0.0;
        chiM_paths[1][t] = 0.0;
        chiL_paths[1][t] = 0.0;

        // Path 2 : monotonic compression + shear
        epsN_paths[2][t] = eps_compression * ramp;
        epsM_paths[2][t] = eps_shearM * ramp;
        epsL_paths[2][t] = eps_shearL * ramp;
        chiN_paths[2][t] = 0.0;
        chiM_paths[2][t] = 0.0;
        chiL_paths[2][t] = 0.0;

        // Path 3 : monotonic tension + shear
        epsN_paths[3][t] = eps_tension * ramp;
        epsM_paths[3][t] = eps_shearM * ramp;
        epsL_paths[3][t] = eps_shearL * ramp;
        chiN_paths[3][t] = 0.0;
        chiM_paths[3][t] = 0.0;
        chiL_paths[3][t] = 0.0;

        // Path 4: monotonic compression + bending
        epsN_paths[4][t] = eps_compression * ramp;
        epsM_paths[4][t] = 0.0;
        epsL_paths[4][t] = 0.0;
        chiN_paths[4][t] = 0.0;
        chiM_paths[4][t] = chi_bendingM * ramp;
        chiL_paths[4][t] = chi_bendingL * ramp;

        // Path 5: monotonic tension + bending
        epsN_paths[5][t] = eps_tension * ramp;
        epsM_paths[5][t] = 0.0;
        epsL_paths[5][t] = 0.0;
        chiN_paths[5][t] = 0.0;
        chiM_paths[5][t] = chi_bendingM * ramp;
        chiL_paths[5][t] = chi_bendingL * ramp;

        // Path 6: monotonic tension + torsion
        epsN_paths[6][t] = eps_tension * ramp;
        epsM_paths[6][t] = 0.0;
        epsL_paths[6][t] = 0.0;
        chiN_paths[6][t] = chi_torsion * ramp;
        chiM_paths[6][t] = 0.0;
        chiL_paths[6][t] = 0.0;

        // Path 7: monotonic tension + torsion + bending
        epsN_paths[7][t] = eps_tension * ramp;
        epsM_paths[7][t] = 0.0;
        epsL_paths[7][t] = 0.0;
        chiN_paths[7][t] = chi_torsion * ramp;
        chiM_paths[7][t] = chi_bendingM * ramp;
        chiL_paths[7][t] = chi_bendingL * ramp;

        // Path 8: monotonic tension + shear + torsion + bending
        epsN_paths[8][t] = eps_tension * ramp;
        epsM_paths[8][t] = eps_shearM * ramp;
        epsL_paths[8][t] = eps_shearL * ramp;
        chiN_paths[8][t] = chi_torsion * ramp;
        chiM_paths[8][t] = chi_bendingM * ramp;
        chiL_paths[8][t] = chi_bendingL * ramp;
    }

    // Comparison of analytical calculations and Chrono implementation for all loading paths
    for (int path = 0 ; path < npaths ; path++) {
        // Initialize history and incremental variables
        double sig = 0.0;
        double epsNeqMax = 0.0;
        double epsTeqMax = 0.0;
        double epsold = 0.0;
        double sigNold = 0.0;
        double sigMold = 0.0;
        double sigLold = 0.0;
        double muNold = 0.0;
        double muMold = 0.0;
        double muLold = 0.0;
        double Wint = 0.0;

        ChVectorN<double, 18> statev;
        statev.setZero();

        if (debug_mode) { // Output to file to plot
            res.open ("results_path"+std::to_string(path)+".csv");
            res << "# Analytical solution from unit test WOOD_MaterialVECTTest - stress_no_eigenstrain.\n";
            res << "# E0="<<E0<<", alpha="<<alpha<<", sigmat="<<sigmat<<", sigmac="<<sigmac<<", sigmas="<<sigmas<<", nt="<<nt<<", lt="<<lt<<", couple_mult="<<couple_multiplier<<"\n";
            res << "epsN, epsM, epsL, epsNeqMax, epsTeqMax, epsmax, sigN, sigM, sigL, sig, muN, muM, muL, Wint, crack\n";
            res << "0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0\n";
        }
        for (int t = 1 ; t < nsteps ; t++) {
            double epsN = epsN_paths[path][t];
            double epsM = epsM_paths[path][t];
            double epsL = epsL_paths[path][t];
            double epsNincr = epsN - epsN_paths[path][t-1];
            double epsMincr = epsM - epsM_paths[path][t-1];
            double epsLincr = epsL - epsL_paths[path][t-1];
            double chiN = chiN_paths[path][t];
            double chiM = chiM_paths[path][t];
            double chiL = chiL_paths[path][t];
            double chiNincr = chiN - chiN_paths[path][t-1];
            double chiMincr = chiM - chiM_paths[path][t-1];
            double chiLincr = chiL - chiL_paths[path][t-1];

            // -------- Analytical calculation -------- //

            double epsT = std::sqrt(epsM * epsM + epsL * epsL);
            double epsNeq = std::sqrt(epsN * epsN +  betaM * rMsq * chiM * chiM + betaL * rLsq * chiL * chiL);
            double epsTeq = std::sqrt(epsM * epsM + epsL * epsL + betaN * rNsq * chiN * chiN);
            double eps = std::sqrt(epsNeq * epsNeq + alpha * epsTeq * epsTeq);            
            epsNeqMax = std::max(epsNeqMax, epsNeq);
            epsTeqMax = std::max(epsTeqMax, epsTeq);
            double epsmax = std::sqrt(epsNeqMax * epsNeqMax + alpha * epsTeqMax * epsTeqMax);
            double w = omega(epsN, epsT);

            // Elastic prediction
            sig += E0 * (eps - epsold);

            // Correction (return)
            double sigbs = sigma_bs(w, eps, epsmax, length);
            if (sig < 0.0) sig = 0.0;
            if (sig > sigbs) sig = sigbs;

            // Stress and torque components
            double sigN(0.0), sigM(0.0), sigL(0.0);
            double muN(0.0), muM(0.0), muL(0.0);
            if (eps > 0.0) { // TODO: is there anything safer ? MIN_DBL ?
                sigN = sig * epsN / eps;
                sigM = sig * alpha * epsM / eps;
                sigL = sig * alpha * epsL / eps;
                muN = sig * betaN * rNsq * alpha * chiN / eps;
                muM = sig * betaM * rMsq * chiM / eps;
                muL = sig * betaL * rLsq * chiL / eps;
            }

            // Total Energy of mechanical deformation (trapezoidal integration)
            Wint += length * facet_area * (epsNincr * 0.5 * (sigNold + sigN) + epsMincr * 0.5 * (sigMold + sigM) + epsLincr * 0.5 * (sigLold + sigL) +
                                           chiNincr * 0.5 * (muNold + muN)  + chiMincr * 0.5 * (muMold + muM) + chiLincr * 0.5 * (muLold + muL));

            // Crack Opening
            double crackN = length * std::max(epsN - sigN / E0, 0.0); // Normal Crack if positive
            double crackM = length * (epsM - sigM / (alpha * E0)); // Is this tangential sliding really a "crack" ?
            double crackL = length * (epsL - sigL / (alpha * E0)); // Is this tangential sliding really a "crack" ?
            double crack = std::sqrt(crackN * crackN + crackM * crackM + crackL * crackL);

            // Update for next increment
            epsold = eps;
            sigNold = sigN;
            sigMold = sigM;
            sigLold = sigL;
            muNold = muN;
            muMold = muM;
            muLold = muL;

            // -------- Chrono implementation -------- //

            ChVector3d strain_increment(epsNincr, epsMincr, epsMincr);
            ChVector3d curvature_increment(chiNincr, chiMincr, chiLincr);
            ChVector3d stress, couple;
            my_mat->ComputeStress(strain_increment,curvature_increment, length, statev, facet_area, facet_width, facet_height, random_field, stress, couple);

            // ---------------------------------------- //
            // -------- Comparison and Testing -------- //
            // ---------------------------------------- //

            if (!debug_mode) {
                double tol = 1e-10;

                // The state variables are updatead inside my_mat->ComputeStress()
                // so they contain the current values.
                // Maybe we should test for the old values before calling my_mat->ComputeStress() ?
                ASSERT_NEAR(epsN, statev(0), tol);
                ASSERT_NEAR(epsM, statev(1), tol);
                ASSERT_NEAR(epsL, statev(2), tol);
                ASSERT_NEAR(sigN, statev(3), tol);
                ASSERT_NEAR(sigM, statev(4), tol);
                ASSERT_NEAR(sigL, statev(5), tol);
                ASSERT_NEAR(epsNeqMax, statev(6), tol);
                ASSERT_NEAR(epsTeqMax, statev(7), tol);
                ASSERT_NEAR(eps, statev(8), tol);
                ASSERT_NEAR(sig, statev(9), tol);
                ASSERT_NEAR(Wint, statev(10), tol);
                ASSERT_NEAR(crack, statev(11), tol);
                ASSERT_NEAR(chiN, statev(12), tol);
                ASSERT_NEAR(chiM, statev(13), tol);
                ASSERT_NEAR(chiL, statev(14), tol);
                ASSERT_NEAR(muN, statev(15), tol);
                ASSERT_NEAR(muM, statev(16), tol);
                ASSERT_NEAR(muL, statev(17), tol);

                ASSERT_NEAR(sigN, stress[0], tol);
                ASSERT_NEAR(sigM, stress[1], tol);
                ASSERT_NEAR(sigL, stress[2], tol);
                ASSERT_NEAR(muN, couple[0], tol);
                ASSERT_NEAR(muM, couple[1], tol);
                ASSERT_NEAR(muL, couple[2], tol);
            }

            if (debug_mode) {// Output
                res<<epsN<<", "<<epsM<<", "<<epsL<<", "<<epsNeqMax<<", "<<epsTeqMax<<", "<<epsmax<<", "<<sigN<<", "<<sigM<<", "<<sigL<<", "<<sig<<", "<<muN<<", "<<muM<<", "<<muL<<", "<<Wint<<", "<<crack<<std::endl;
            }
        }

        res.close();
    }
}

