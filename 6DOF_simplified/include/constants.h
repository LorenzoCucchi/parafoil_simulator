#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <C:\Program Files\eigen-3.4.0\Eigen\Dense>
#include "forces.h"
#include "tools.h"

namespace constants
{
    static constexpr double pi { 3.14159 };
    static constexpr double g {9.81};
    static constexpr double rho {1.2};
    static constexpr double Mpay {135};
    static constexpr double Mpar { 13};
    static constexpr double Mtot {Mpay+Mpar};
    static constexpr double IMtot {1/Mtot};
    static constexpr double spay { 0.5};
    static constexpr double spar { 7.5};
    static constexpr double xb {0.5};
    static constexpr double yb {0.5};
    static constexpr double zb {0.5};
    static constexpr double b { 7};
    static constexpr double c { 3};
    static constexpr double thick { 0.3};
    static constexpr double a { 0.8};
    static constexpr double Sw {21.0};
    static constexpr double Sp {0.5};
    static constexpr double Kc {1000000.0};
    static constexpr double Cc {0.5};
    static const Eigen::Matrix<double, 3,3> I {Ip()+Ib()};
    static const Eigen::Matrix<double, 3,3> I_i {I.inverse()};
    static const Eigen::Vector3d Xgp {-1.3, 0.0, -7.4};
    static const Eigen::Vector3d Xgb {0, 0, 0.66};
    static const Eigen::Matrix<double, 3, 3> Rgp {Omrot(Xgp)};
    static const Eigen::Matrix<double, 3, 3> Rgb {Omrot(Xgb)};

    // Aerodynamic coefficients
    static constexpr double CL0 { 0.4};
    static constexpr double CLa { 2.0};
    static constexpr double CD0 { 0.15};
    static constexpr double CDa { 1.0};
    static constexpr double Clp { -0.1};
    static constexpr double Clphi { -0.05};
    static constexpr double Cmq { -2.0};
    static constexpr double Cm0 { 0.018};
    static constexpr double Cnr { -0.07};
    static constexpr double Cma { -0.2};
    static constexpr double CLda { 0.0001};
    static constexpr double CLds { 0.21};
    static constexpr double CDda { 0.0001};
    static constexpr double CDds { 0.3};
    static constexpr double Clda { 0.0021};
    static constexpr double Cnda { 0.004};
}



#endif // CONSTANTS_H