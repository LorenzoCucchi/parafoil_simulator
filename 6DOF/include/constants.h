#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <Eigen/Dense>
#include "forces.h"
#include "tools.h"

namespace constants
{
    extern const double pi;
    extern const double g;
    extern const double rho;
    extern const double Mpay;
    extern const double Mpar;
    extern const double spay;
    extern const double spar;
    extern const double xb;
    extern const double yb;
    extern const double zb;
    extern const double b;
    extern const double c;
    extern const double t;
    extern const double a;
    extern const double Sw;
    extern const double Sp;
    extern const double Kc;
    extern const double Cc; 
    extern const Eigen::Matrix<double,3,3> eye;
    extern const Eigen::Matrix<double, 3, 3> I;
    extern const Eigen::Matrix<double, 3, 3> I_i;
    extern const Eigen::Matrix<double, 3, 3> Mf;
    extern const Eigen::Matrix<double, 3, 3> If;
    extern const Eigen::Vector3d Xgp;
    extern const Eigen::Vector3d Xgb;
    extern const Eigen::Matrix<double, 3, 3> Rgp;
    extern const Eigen::Matrix<double, 3, 3> Rgb;
    extern const Eigen::Matrix<double, 3, 3> A1;
    extern const Eigen::Matrix<double, 3, 3> A2;
    extern const Eigen::Matrix<double, 3, 3> A3;
    extern const Eigen::Matrix<double, 3, 3> A4;
    extern const Eigen::Matrix<double, 6, 6> A;

    // Aerodynamic coefficient
    extern const double CL0;
    extern const double CLa;
    extern const double CD0;
    extern const double CDa;
    extern const double Clp;
    extern const double Clphi;
    extern const double Cmq;
    extern const double Cm0;
    extern const double Cnr;
    extern const double Cma;
    extern const double CLda;
    extern const double CLds; 
    extern const double CDda; 
    extern const double CDds;
    extern const double Clda;
    extern const double Cnda;
}

#endif // CONSTANTS_H