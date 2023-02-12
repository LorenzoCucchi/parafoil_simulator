#include "constants.h"
#include <Eigen/Dense>
#include "forces.h"

using namespace Eigen;

namespace constants
{
    // actual global variables
    extern const double pi { 3.14159 };
    extern const double g {9.81};
    extern const double rho {1.2};
    extern const double Mpay {135};
    extern const double Mpar { 13};
    extern const double spay { 0.5};
    extern const double spar { 7.5};
    extern const double xb {0.5};
    extern const double yb {0.5};
    extern const double zb {0.5};
    extern const double b { 7};
    extern const double c { 3};
    extern const double t { 0.3};
    extern const double a { 0.8};
    extern const double Sw {21.0};
    extern const double Sp {0.5};
    extern const double Kc {1000000.0};
    extern const double Cc {0.5};
    extern const Eigen::Matrix<double, 3, 3> I {Ip()+Ib()};
    extern const Eigen::Matrix<double, 3,3> I_i {I.inverse()};
    extern const Eigen::Vector3d Xgp {-1.3, 0.0, -7.4};
    extern const Eigen::Vector3d Xgb {0, 0, 0.66};
    // Aerodynamic coefficients
    extern const double CL0 { 0.4};
    extern const double CLa { 2.0};
    extern const double CD0 { 0.15};
    extern const double CDa { 1.0};
    extern const double Clp { -0.1};
    extern const double Clphi { -0.05};
    extern const double Cmq { -2.0};
    extern const double Cm0 { 0.018};
    extern const double Cnr { -0.07};
    extern const double Cma { -0.2};
    extern const double CLda { 0.0001};
    extern const double CLds { 0.21};
    extern const double CDda { 0.0001};
    extern const double CDds { 0.3};
    extern const double Clda { 0.0021};
    extern const double Cnda { 0.004};
    
}