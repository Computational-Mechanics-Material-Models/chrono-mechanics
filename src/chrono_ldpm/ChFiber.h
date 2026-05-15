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
// Authors: Ke Yu
//          
// =============================================================================
// Fiber class for LDPM and CSL elements 
//
// A description of the material parameter can be found in: https://doi.org/10.1016/j.cemconcomp.2011.02.011
// =============================================================================

#ifndef CHFIBER_H
#define CHFIBER_H

#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <stdexcept>
//#include <boost/math/tools/roots.hpp> 

#include "chrono_ldpm/ChLdpmApi.h"
#include "chrono/core/ChMatrix33.h"

#include <vector>
#include <string>
//#include "chrono/fea/ChElementBeam.h"

namespace chrono {
namespace ldpm {

/// @addtogroup fea_elements
/// @{

/// Definition of materials to be used for CSL beams and LDPM tets utilizing the lattice discrete particle model.
class ChLdpmApi ChFiber {
  public:
    /// Construct an isotropic elastic material. 
    ChFiber();

    /*
    //copy constructor:
	ChMaterialVECT(const ChMaterialVECT& my_material);
	*/

    // Destructor declared:
    ~ChFiber();


    /// Compute stresses from given strains and state variables.

    /*
    std::tuple<double, double, double> ChFiber::FRP(ChVector3d w, double Pf0,
                                                    ChVector3d nf,
                                                    ChVector3d n,
                                                    double ksp,
                                                    double ksn,
                                                    double df,
                                                    double Ef,
                                                    double t0,
                                                    double Gd,
                                                    double Ls,
                                                    double Ll,
                                                    double beta);

    std::tuple<double, double, double> ChFiber::Pf_Cal(ChVector3d w, double Pf0,
                                                       ChVector3d nf,
                                                       ChVector3d n,
                                                       double ksp,
                                                       double ksn,
                                                       double df,
                                                       double Ef,
                                                       double t0,
                                                       double Gd,
                                                       double Ls,
                                                       double Ll,
                                                       double beta);

    std::tuple<double, double, double> ChFiber::compatibility(ChVector3d w_prime,
                                                              double sf,
                                                              double Ef,
                                                              double df,
                                                              double t0,
                                                              double Gd,
                                                              double Ls,
                                                              double Ll,
                                                              double beta);
    */
    static std::tuple<ChVector3d, double, double, double, double> FRP(ChVector3d w,
                                                              double Pf0,
                                                              ChVector3d nf,
                                                              ChVector3d n,
                                                              double ksp,
                                                              double ksn,
                                                              double df,
                                                              double Ef,
                                                              double t0,
                                                              double Gd,
                                                              double Ls,
                                                              double Ll,
                                                              double beta,
                                                              double vmaxs,
                                                              double vmaxl,
                                                              double pmaxs,
                                                              double pmaxl);

    static std::tuple<ChVector3d, double, double, double> Pf_cal(ChVector3d w,
                                                                 double Pf0,
                                                                 ChVector3d nf,
                                                                 ChVector3d n,
                                                                 double ksp,
                                                                 double ksn,
                                                                 double df,
                                                                 double Ef,
                                                                 double t0,
                                                                 double Gd,
                                                                 double Ls,
                                                                 double Ll,
                                                                 double beta,
                                                                 double vmaxs,
                                                                 double vmaxl,
                                                                 double pmaxs,
                                                                 double pmaxl);   

    static double
    Pv(double v, double Ef, double df, double t0, double Gd, double Le, double beta, double vmax, double pmax);
	
	static double
    ComputeFiberTangent(double v, double Ef, double df, double t0, double Gd, double Le, double beta, double vmax, double pmax);

    static double error(double x,
                        double w,
                        double sf,
                        double Ef,
                        double df,
                        double t0,
                        double Gd,
                        double Ls,
                        double Ll,
                        double beta,
                        double vmaxs,
                        double vmaxl,
                        double pmaxs,
                        double pmaxl);

    static double solve_toms748(double w, double sf, double Ef, double df, double t0, double Gd, double Ls, double Ll, double beta);

    static double solve_bisect(double w, double sf, double Ef, double df, double t0, double Gd, double Ls, double Ll, double beta);

    static std::tuple<double, double, double> compatibility(double w,
                                                            double sf,
                                                            double Ef,
                                                            double df,
                                                            double t0,
                                                            double Gd,
                                                            double Ls,
                                                            double Ll,
                                                            double beta,
                                                            double vmaxs,
                                                            double vmaxl,
                                                            double pmaxs,
                                                            double pmaxl);

    static double solve_newton_raphson(double w, double sf, double Ef, double df, double t0, double Gd, double Ls, double Ll, double beta);

    
  private:
    

  //public:
  //EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

/// @} fea_elements

}  // end of namespace ldpm
}  // end of namespace chrono

#endif
