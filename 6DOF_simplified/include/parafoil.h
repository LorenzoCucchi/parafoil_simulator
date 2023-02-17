#ifndef INC_6DOF_SIMPLIFIED_PARAFOIL_H
#define INC_6DOF_SIMPLIFIED_PARAFOIL_H

#include <boost\numeric\odeint.hpp>
#include <C:\Program Files\eigen-3.4.0\Eigen\Dense>
#include "constants.h"
#include "json.hpp"

using namespace std;
using namespace Eigen;
using namespace boost::numeric::odeint;
using json = nlohmann::json;

typedef Matrix<double, 15, 1> state_type;

class Parafoil {
public:
    Parafoil(const string &path);

    void print();

    Matrix<double, 3, 3> Inertia_calc() const;

    Matrix<double, 3, 3> Omrot(const Vector3d &vec);

    Vector3d Fa_b(const Vector3d &vec);

    Vector3d Fa_w(const Vector3d &vec);

    Vector3d W(const Vector3d &vec);

    Vector3d Ma_w(Vector3d &vel, Vector3d &w, Vector3d rot);

    double CD_b(const Vector3d &vec);

    double CL_w(const Vector3d &vec);

    Matrix<double, 3, 2> SFa(const Vector3d &vec, double delta);

    Matrix<double, 3, 2> SMa(const Vector3d &vec);

    Matrix<double, 3, 3> Wrot(const Vector3d &vec);

    Matrix<double, 3, 3> QuatToAtt(const Matrix<double, 4, 1> &q);

    Vector3d QuatToEuler(const Matrix<double, 4, 1> &q);

    Matrix<double, 4, 4> Omega(const Vector3d &vel);

    static Matrix<double, 4, 1> EulToQuat(const Vector3d &v);

    void my_system(const state_type &x, state_type &dxdt,double t);

    void simulate();

    void simulate_control();

private:
    double mpay;
    double mpar;
    double mtot;
    double imto;
    double xb;
    double zb;
    double b;
    double c;
    double th;
    double a;
    double sw;
    double sp;
    Eigen::Matrix<double, 3, 3> inv_inertia;
    Eigen::Vector3d xgp;
    Eigen::Vector3d xgb;
    Eigen::Matrix<double, 3, 3> rgp;
    Eigen::Matrix<double, 3, 3> rgb;
    double CL0;
    double CLa;
    double CD0;
    double CDa;
    double Clp;
    double Clphi;
    double Cmq;
    double Cm0;
    double Cnr;
    double Cma;
    double CLda;
    double CLds;
    double CDda;
    double CDds;
    double Clda;
    double Cnda;
    state_type initial_state;
    double dt;
    double fin_alt;
    Vector3d wind;
    Vector2d target;

};


#endif //INC_6DOF_SIMPLIFIED_PARAFOIL_H
