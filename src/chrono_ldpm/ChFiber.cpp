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

#include "chrono_ldpm/ChFiber.h"

namespace chrono {
namespace ldpm {

// Construct an isotropic material.

ChFiber::ChFiber()
    {
        
}



// Destructor
ChFiber::~ChFiber()	{};

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

/*
std::tuple<double, double, double> ChFiber::FRP(ChVector3d w,double Pf0,
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
                                                double beta) {
    double e = 1;

    while (abs(e) > 0.001) {
        auto [Pf, vs, vl] = Pf_cal(ChVector3d w, Pf0, nf, n, ksp, ksn, df, Ef, t0, Gd, Ls, Ll, beta);

        e = Pf - Pf0;
        Pf0 = Pf;
    }

    return td::make_tuple(Pf, vs, vl);
}



std::tuple<double, double, double> ChFiber::Pf_Cal(ChVector3d w,
     double Pf0, ChVector3d nf, ChVector3d n, double ksp, double ksn, double df, double Ef, double t0, double Gd, double
Ls, double Ll, double beta) {

    ChVector3d w(w_N, w_M, w_L);

    double theta = acos(nf ^ n);
    double Pfn = Pf * sin(theta);  // updated
    //print(theta)
    //double sf = Pfn*math.sin(theta/2)/(ksp*sigma_t*df*math.cos(theta/2)*math.cos(theta/2))
    double sf = 0;

    ChVector3d w_prime = w + 2 * sf * nf;
    if (w_prime.Length() == 0) {
        ChVector3d nf_prime = nf;
    } else {
        ChVector3d nf_prime = w_prime / (w_prime.Length());
    }

    double phi_prime = acos(nf ^ nf_prime);
    double Ps;
    double vs;
    double vl;
    [ Ps, vs, vl ] = compatibility(w_prime, sf, Ef, df, t0, Gd, Ls, Ll, beta);
    double Pf_prime = exp(ksn * phi_prime) * Ps;
    return std::make_tuple(Pf_prime, vs, vl);
}


std::tuple<double, double, double> ChFiber::compatibility(ChVector3d w_prime, double sf, double Ef, double df, double
t0, double Gd, double Ls, double Ll, double beta): { double a = 0; double b = 1; double x0; double fa = error(a,
w_prime, sf, Ef, df, t0, Gd, Ls, Ll, beta); double fb = error(b, w_prime, sf, Ef, df, t0, Gd, Ls, Ll, beta);

    if (abs(fa) < 0.000001) {
         x0= 0;
    }
    else if (abs (fb) < 0.000001) {
        x0 = 1;
    }
    else{
        while (a <= b) {
        x0 = (a + b) / 2;
        double fx0 = error(x0, w_prime, sf, Ef, df, t0, Gd, Ls, Ll, beta);
        if (abs(fx0) < 0.000001) {
            break;
        }

        if (fa * fx0 < 0) {
            b = x0;
            fb = fx0;
        } else if (fb * fx0 < 0) {
            a = x0;
            fa = fx0;
        }

        if (b - a < 0.000001) {
            std::cout << "no root for x"<< std::endl;
        }
        }
    }

    double vs = x0 * (w_prime.Length() - 2 * sf);
    double vl = (1 - x0) * (w_prime.Length() - 2 * sf);
    double Ps = Pv(vs, Ef, df, t0, Gd, Ls - sf, beta);
    return std::make_tuple(Ps, vs, vl);
}
*/

std::tuple<ChVector3d, double, double, double,double> ChFiber::FRP(ChVector3d w,
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
                                                            double pmaxl) {
    double e = 1;
    ChVector3d Pf_result;
    double Ps_result = 0.0;
    double vs_result = 0.0;
    double vl_result = 0.0;

    while (abs(e) > 0.001) {
        auto [Pf, Ps, vs, vl] =
            Pf_cal(w, Pf0, nf, n, ksp, ksn, df, Ef, t0, Gd, Ls, Ll, beta, vmaxs, vmaxl, pmaxs, pmaxl);

        double Pf_len = Pf.Length();
        e = Pf_len - Pf0;
        Pf0 = Pf_len;

        Pf_result = Pf;
        Ps_result = Ps;
        vs_result = vs;
        vl_result = vl;
    }
	
	double kfs = ComputeFiberTangent(vs_result, Ef, df, t0, Gd, Ls, beta, vmaxs, pmaxs);
    double kfl = ComputeFiberTangent(vl_result, Ef, df, t0, Gd, Ll, beta, vmaxl, pmaxl);
	double kf = (kfl + kfs)/1;
    return std::make_tuple(Pf_result, Ps_result, vs_result, vl_result,kf);
}

double ChFiber::ComputeFiberTangent(double v, double Ef, double df, double t0, double Gd, double Le, double beta, double vmax, double pmax) {
    const double v0 = 1e-6;  
    const double pi = CH_PI; 

    if (v < vmax) {
        return pmax / vmax;
    }
	double vd = 2 * t0 * Le * Le / Ef / df + sqrt(8 * Gd * Le * Le / Ef / df);
    if (v < vd) {   
        double tmp = pi * pi * Ef * df * df * df / 2.0;

        double Gd_modifier = Gd * (1 - exp(-v / v0));
        double P = sqrt(tmp * (t0 * v + Gd_modifier));
        if (P < 1e-12) return 0; 
        double dGd_dv = Gd * exp(-v / v0) / v0;
        double dPv_dv = (0.5 / P) * tmp * (t0 + dGd_dv);

        return dPv_dv;
    }
    else {
        
        double p0 = pi * Le * df * t0;

        double term1 = - (1.0 / Le) * (1.0 + beta * (v - vd) / df);
        double term2 = (beta / df) * (1.0 - (v - vd) / Le);

        double dPv_dv = p0 * (term1 + term2);

        return dPv_dv > 0 ? dPv_dv : 0.0;
    }
}

/*
ChVector3d ChFiber::FRP(ChVector3d w, double Pf0, ChVector3d nf,
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
                        double pmaxl) {
    double e = 1;
    ChVector3d Pf;

    while (abs(e) > 0.001) {
        //std::cout << "initial Pf0 " << Pf0 << std::endl;
        auto [Pf, Ps, vs, vl] =
            Pf_cal(w, Pf0, nf, n, ksp, ksn, df, Ef, t0, Gd, Ls, Ll, beta, vmaxs, vmaxl, pmaxs, pmaxl);
        double Pf_len = Pf.Length();
        //std::cout << "Pf_len " << Pf_len << std::endl;
        e = Pf_len - Pf0;
        Pf0 = Pf_len;
    }
    //std::cout << "final Pf_len " << Pf.Length() << std::endl;
    return Pf;
}
*/




std::tuple<ChVector3d, double, double, double> ChFiber::Pf_cal(ChVector3d w,
                           double Pf,
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
                           double pmaxl) {
   	
	//ChVector3d w(w_N, w_M, w_L);

    double theta = acos(nf.Dot(n));
    double Pfn = Pf * nf.Dot(n);  // updated
    //print(theta)
    //double sf = Pfn*math.sin(theta/2)/(ksp*sigma_t*df*math.cos(theta/2)*math.cos(theta/2))
    double sf = 0;

	ChVector3d w_prime = w + 2 * sf * nf;
    if (w_prime.Length() < 1e-16) {
        //ChVector3d nf_prime = nf;
        throw std::runtime_error("Error when calculating fiber bridging force vector: w_prime=0 \n");
        exit(EXIT_FAILURE);
    }
        
    ChVector3d nf_prime = w_prime / (w_prime.Length());
    
    double dot = nf.Dot(nf_prime);
    dot = std::max(-1.0, std::min(1.0, dot));
    double phi_prime = acos(dot);


    // Ps = solve_toms748(w_prime.Length(), sf, Ef, df, t0, Gd, Ls, Ll, beta);
    //double Ps = solve_bisect(w_prime.Length(), sf, Ef, df, t0, Gd, Ls, Ll, beta);
    auto [Ps, vs, vl] = compatibility(w_prime.Length(), sf, Ef, df, t0, Gd, Ls, Ll, beta, vmaxs, vmaxl, pmaxs, pmaxl);
    //double Ps = solve_newton_raphson(w_prime.Length(), sf, Ef, df, t0, Gd, Ls, Ll, beta);

    ChVector3d Pf_prime = (exp(ksn * phi_prime) * Ps) * nf_prime;
    return std::make_tuple(Pf_prime, Ps, vs, vl);
}




double ChFiber::Pv(double v, double Ef, double df, double t0, double Gd, double Le, double beta, double vmax, double pmax){
    double vd = 2 * t0 * Le * Le / Ef / df + pow(8 * Gd * Le * Le / Ef / df, 0.5);

    double P;

    if (v < vmax) {
        return pmax * pow(v / vmax, 1);
    }
    
    if (v < vd) {
        double tmp = CH_PI * CH_PI * Ef * df * df * df / 2;
		double v0 = 1e-6;
		double Gd_modifier = Gd * (1-exp(-v/v0));
        P = pow(tmp * (t0 * v + Gd_modifier), 0.5);
    }
        
    else {
        double p0 = CH_PI * Le * df * t0;
        P = p0 * (1 - (v - vd) / Le) * (1 + beta * (v - vd) / df);
    
    }
    if (P < 0) {
        std::cout << "v= " << v << " vd= " << vd << " Le= " << Le << std::endl;
    
    }
    return P;
}



double ChFiber::error(double x,
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
                      double pmaxl) {
    //double w = w_prime.Length();
    double vs = x * (w - 2 * sf);
    double vl = (1 - x) * (w - 2 * sf);

    double Ps = Pv(vs, Ef, df, t0, Gd, Ls - sf, beta, vmaxs, pmaxs);

    double Pl = Pv(vl, Ef, df, t0, Gd, Ll - sf, beta, vmaxl, pmaxl);

    double e = Ps - Pl;
    return e;
}

/*
double ChFiber::solve_toms748(double w, double sf, double Ef, double df, double t0, double Gd, double Ls, double Ll,
double beta) { auto f = [&](double x) { return error(x, w, sf, Ef, df, t0, Gd, Ls, Ll, beta); };

    double a = 0.0;
    double b = 1.0;
    auto tol = [](double l, double r) { return std::abs(r - l) < 1e-6; };
    boost::uintmax_t max_iter = 1000;

    double fa = f(a);
    double fb = f(b);

    if (fa * fb > 0) {
        throw std::runtime_error("solve_toms748: Root not bracketed in [0, 1].");
        exit(EXIT_FAILURE);
    }

    auto result = boost::math::tools::toms748_solve(f, a, b, tol, max_iter);

    double x0 = (result.first + result.second) / 2.0;
    double vs = x0 * (w - 2 * sf);
    double Ps = Pv(vs, Ef, df, t0, Gd, Ls - sf, beta);
    return Ps;
}

double ChFiber::solve_bisect(double w, double sf, double Ef, double df, double t0, double Gd, double Ls, double Ll,
double beta) { auto f = [&](double x) { return error(x, w, sf, Ef, df, t0, Gd, Ls, Ll, beta); };

    double a = 0.0;
    double b = 1.0;
    double tol = 1e-6;
    boost::uintmax_t max_iter = 1000;

    double fa = f(a);
    double fb = f(b);

    if (fa * fb > 0) {
        throw std::runtime_error("solve_bisect: Root not bracketed in [0, 1].");
        exit(EXIT_FAILURE);
    }

    auto result = boost::math::tools::bisect(
        f, a, b, [=](double l, double r) { return std::abs(r - l) <= tol; }, max_iter);

    double x0 = (result.first + result.second) / 2.0;
    double vs = x0 * (w - 2 * sf);
    double Ps = Pv(vs, Ef, df, t0, Gd, Ls - sf, beta);
    return Ps;

}
*/



std::tuple<double, double, double> ChFiber::compatibility(double w,
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
                                                          double pmaxl) {
    double a = 0;
    double b = 1;
    double x0;

    if (Ls < w) {
        b = Ls / (w - 2 * sf);
    }
    double fa = error(a, w, sf, Ef, df, t0, Gd, Ls, Ll, beta, vmaxs, vmaxl, pmaxs, pmaxl);
    double fb = error(b, w, sf, Ef, df, t0, Gd, Ls, Ll, beta, vmaxs, vmaxl, pmaxs, pmaxl);

    double tol = 1e-5;


     if (std::abs(fa) < tol) {
        double vs = a * (w - 2 * sf);
        double vl = (1 - a) * (w - 2 * sf);
        double Ps = Pv(vs, Ef, df, t0, Gd, Ls - sf, beta, vmaxs, pmaxs);
        return std::make_tuple(Ps, vs, vl);
    }

    if (std::abs(fb) < tol) {
        double vs = b * (w - 2 * sf);
        double vl = (1 - b) * (w - 2 * sf);
        double Ps = Pv(vs, Ef, df, t0, Gd, Ls - sf, beta, vmaxs, pmaxs);
        return std::make_tuple(Ps, vs, vl);
    }

    if (fa * fb > 0) {
        std::cout << "No root found because fa0 * fb0 > 0, w="<<w << std::endl;
        double vs_a = a * (w - 2 * sf);
        double vl_a = (1 - a) * (w - 2 * sf);
        double Ps_a = Pv(vs_a, Ef, df, t0, Gd, Ls - sf, beta, vmaxs, pmaxs);
        double Pl_a = Pv(vl_a, Ef, df, t0, Gd, Ll - sf, beta, vmaxl, pmaxl);
        std::cout << "vmaxs = " << vmaxs << " vmaxl = " << vmaxl << std::endl;
        std::cout << "pmaxs = " << pmaxs << " pmaxl = " << pmaxl << std::endl;
        std::cout << "vs_a = " << vs_a << " vl_a = " << vl_a << std::endl;
        std::cout << "Ps_a = " << Ps_a << " Pl_a = " << Pl_a << std::endl;
        double vs_b = b * (w - 2 * sf);
        double vl_b = (1 - b) * (w - 2 * sf);
        double Ps_b = Pv(vs_b, Ef, df, t0, Gd, Ls - sf, beta, vmaxs, pmaxs);
        double Pl_b = Pv(vl_b, Ef, df, t0, Gd, Ll - sf, beta, vmaxl, pmaxl);
        std::cout << "vs_b = " << vs_b << " vl_b = " << vl_b << std::endl;
        std::cout << "Ps_b = " << Ps_b << " Pl_b = " << Pl_b << std::endl;
        return std::make_tuple(0, Ls, (w - 2 * sf) - Ls);
        //throw std::runtime_error("No root found because fa0 * fb0 > 0");
        //std::cout<< "No root found because fa0 * fb0 > 0"<< std::endl;
    }

    double fx0;
    while ((b - a) > 1e-14) {
        x0 = (a + b) / 2.0;
        fx0 = error(x0, w, sf, Ef, df, t0, Gd, Ls, Ll, beta, vmaxs, vmaxl, pmaxs, pmaxl);

        if (abs(fx0) < tol) {
            break;  
        }
        if (fa * fx0 < 0) { 
            b = x0;
            fb = fx0;
        } 
        else if (fb * fx0 < 0) {
            a = x0;
            fa = fx0;
        } else {
            b = x0;
            fb = fx0;
            std::cout << " No root found because fa * fb > 0 " << std::endl;
            std::cout << " a= " << a << " b= "<<b << std::endl;
            double vs_a = a * (w - 2 * sf);
            double vl_a = (1 - a) * (w - 2 * sf);
            double Ps_a = Pv(vs_a, Ef, df, t0, Gd, Ls - sf, beta, vmaxs, pmaxs);
            double Pl_a = Pv(vl_a, Ef, df, t0, Gd, Ll - sf, beta, vmaxl, pmaxl);
            std::cout << "Ps_a = " << Ps_a << " Pl_a = " << Pl_a << std::endl;
            double vs_b = b * (w - 2 * sf);
            double vl_b = (1 - b) * (w - 2 * sf);
            double Ps_b = Pv(vs_b, Ef, df, t0, Gd, Ls - sf, beta, vmaxs, pmaxs);
            double Pl_b = Pv(vl_b, Ef, df, t0, Gd, Ll - sf, beta, vmaxl, pmaxl);
            std::cout << "Ps_b = " << Ps_b << " Pl_b = " << Pl_b << std::endl;
            //throw std::runtime_error("No root found because fa * fb > 0");
        }

    }

    if (abs(fx0) > tol) {
        std::cout << " No root found for x within tolerance. " << std::endl;
        std::cout << " a= " << a << " b= "<< b
                  << " error= " << abs(fx0) << std::endl;
        throw std::runtime_error("No root found for x within tolerance.");
        exit(EXIT_FAILURE);
    }


    double vs = x0 * (w - 2 * sf);
    double vl = (1-x0) * (w - 2 * sf);
    //double vl = (1 - x0) * (w - 2 * sf);
    double Ps = Pv(vs, Ef, df, t0, Gd, Ls - sf, beta, vmaxs, pmaxs);
    return std::make_tuple(Ps, vs, vl);
}

/*
double ChFiber::solve_newton_raphson(double w, double sf, double Ef, double df, double t0, double Gd, double Ls, double
Ll, double beta) { auto f = [&](double x) { return error(x, w, sf, Ef, df, t0, Gd, Ls, Ll, beta); }; auto d_f =
[&](double x) { return (f(x + 1e-6) - f(x)) / 1e-6; };  // Numerical derivative

    double x0 = 0.5;  // Initial guess
    double tol = 1e-6;
    int max_iter = 1000;

    for (int i = 0; i < max_iter; ++i) {
        double fx = f(x0);
        double dfx = d_f(x0);

        if (std::isnan(dfx) || std::abs(dfx) < 1e-12) {
            std::cout << " Invalid derivative dfx. " << dfx << std::endl;
            throw std::runtime_error("Invalid derivative dfx.");
        }

        if (std::abs(fx) < tol) {
            break;
        }
        x0 = x0 - fx / dfx;
    }

    if (abs(f(x0)) > tol) {
        std::cout << " No root found for x within tolerance. " << std::endl;
        throw std::runtime_error("No root found for x within tolerance.");
        exit(EXIT_FAILURE);
    }

    double vs = x0 * (w - 2 * sf);
    if (std::isnan(vs)) {
        std::cout << " vs is NaN " << std::endl;
        throw std::runtime_error("vs is NaN");
    }

    double Ps = Pv(vs, Ef, df, t0, Gd, Ls - sf, beta);
    return Ps;
}
*/


}  // end of namespace ldpm
}  // end of namespace chrono
