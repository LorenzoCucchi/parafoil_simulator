#include "forces.h"
#include <iostream>
#include "constants.h"
#include <C:\Program Files\eigen-3.4.0\Eigen\Dense>

using namespace std;
using namespace Eigen;
using namespace constants;

typedef Matrix<double, 3, 2> Matrix32d;

// Aerodynamic forces on the payload
Vector3d Fa_b( Vector3d& vec)
{
    Vector3d f;
    f = -0.5 * rho * Sp * vec.norm() * CD_b(vec) * vec;
    return f;
}

double CD_b( Vector3d& vec)
{
    double alpha = atan2(vec(2), vec(0));
    //if (alpha >= (12.0*pi/180.0))  alpha = (12.0*pi/180.0);
    double cd = CD0 + CDa * alpha*alpha;
    return cd;
}

// Weight forces on the payload
Vector3d W_b( Vector3d& vec)
{
    Vector3d W;
    W(0) = Mpay * g * (-sin(vec(1)));
    W(1) = Mpay * g * (sin(vec(0)) * cos(vec(1)));
    W(2) = Mpay * g * (cos(vec(0)) * cos(vec(1)));
    return W;
}

// Aerodynamic forces on the wing
Vector3d Fa_w( Vector3d& vec)
{
    double q = 0.5 * rho * Sw * vec.norm();
    Vector3d f;
    f(0) = q * (CL_w(vec) * vec(2) - CD_b(vec) * vec(0));
    f(1) = -q * CD_b(vec) * vec(1);
    f(2) = q * (-CL_w(vec) * vec(0) - CD_b(vec) * vec(2));
    return f;
}

double CL_w( Vector3d& vec)
{
    double alpha = atan2(vec(2), vec(0));
    //if (alpha >= (12.0*pi/180.0))  alpha = (12.0*pi/180.0);
    double cl = CL0 + CLa * alpha;
    return cl;
}

// Weight forces on the wing
Vector3d W_w( Vector3d& vec)
{
    Vector3d W;
    W(0) = Mpar * g * (-sin(vec(1)));
    W(1) = Mpar * g * (sin(vec(0)) * cos(vec(1)));
    W(2) = Mpar * g * (cos(vec(0)) * cos(vec(1)));
    return W;
}

// Aerodynamic moments
Vector3d Ma_w( Vector3d& vel, Vector3d& w, Vector3d rot){
    double q = 0.5*rho*Sw*vel.squaredNorm();
    double vnorm = vel.norm();
    Vector3d m;
    m(0) = q*(Clp*pow(b,2)*w(0)/(2*vnorm) + Clphi*b*rot(0));
    m(1) = q*(Cmq*pow(c,2)*w(1)/(2*vnorm) + Cm0*c + Cma*c*atan2(vel(2), vel(0)));
    m(2) = q*(Cnr*pow(b,2)*w(2)/(2*vnorm));
    return m;
}

// Apparent forces
Matrix3d MF()
{
    double kA = 0.848 * pi / 4;
    double kB = 0.339 * pi / 4;
    double kC = (b / c) / (1 + b / c) * pi/4;
    double A = kA * rho * thick*thick * b * (1+ 8.0/3.0 *a*a*a);
    double B = kB * rho * thick*thick + 2*a*a*(1-thick*thick)*c;
    double C = kC * rho * c*c * b * sqrt(1+2*a*a*(1-thick*thick));
    Matrix3d M;
    M << A, 0, 0,
          0, B, 0,
          0, 0, C;
        
    return M;
}

Matrix3d IF()
{
    double AR = b/c;
    double kA = 0.055*AR/(1+AR);
    double kB = 0.0308*AR/(1+AR);
    double kC = 0.0555;
    double A = kA * rho * pow(c, 2) * pow(b, 3);
    double B = kB * rho * pow(c, 4) * b * (1 + pi/6 * (1+AR)*AR*pow(a,2)*pow(thick,2));
    double C = kC * rho * pow(thick, 2) * pow(b, 3) * (1+8*pow(a,2));
    Matrix3d Ia;
    Ia << A, 0, 0,
          0, B, 0,
          0, 0, C;
        
    return Ia;
}

Matrix32d SFa(const Vector3d& vec, double delta){

    double sign;
    if (delta >= 0.0)  sign = 1;
    else sign = -1;

    double q = 0.5*rho*Sw*vec.norm();
    Matrix32d F;
    F(0,0) = q*(CLda*vec(2)-CDda*vec(0))*sign;
    F(0,1) = q*(CLds*vec(2)-CDds*vec(0));
    F(1,0) = q*(-CDda*vec(1)*sign);
    F(1,1) = q*(-CDds*vec(1));
    F(2,0) = q*(-CLda*vec(0)-CDda*vec(2))*sign;
    F(2,1) = q*(-CLds*vec(0)-CDds*vec(2));
    return F;
}

Matrix32d SMa(const Vector3d& vec){

    double q = 0.5*rho*Sw*vec.squaredNorm()*b/thick;
    Matrix32d M;
    M(0,0) = q*Clda;
    M(0,1) = 0.0;
    M(1,0) = 0.0;
    M(1,1) = 0.0;
    M(2,0) = q*Cnda;
    M(2,1) = 0.0;

    return M;
}

Matrix3d Ib(){
    Matrix3d Ia;
    Ia(0,0) = pow(yb,2)+pow(zb,2);
    Ia(0,1) = 0.0;
    Ia(0,2) = 0.0;
    Ia(1,0) = 0.0;
    Ia(1,1) = pow(zb,2)+pow(xb,2);
    Ia(1,2) = 0.0;
    Ia(2,0) = 0.0;
    Ia(2,1) = 0.0;
    Ia(2,2) = pow(yb,2)+pow(xb,2);
    Ia = Ia*Mpay/12;
    return Ia;
}

Matrix3d Ip(){
    Matrix3d I;
    I(0,0) = pow(b,2)+pow(thick,2);
    I(0,1) = 0.0;
    I(0,2) = 0.0;
    I(1,0) = 0.0;
    I(1,1) = pow(c,2)+pow(thick,2);
    I(1,2) = 0.0;
    I(2,0) = 0.0;
    I(2,1) = 0.0;
    I(2,2) = pow(b,2)+pow(c,2);
    I = I*Mpar/12;
    return I;
}


