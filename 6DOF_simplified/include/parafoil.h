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

typedef Matrix<float, 13, 1> state_type;

class Parafoil {
public:
    Parafoil(const string &path);

    void print();

    Matrix<float, 3, 3> Inertia_calc();

    Matrix<float, 3, 3> Omrot(const Vector3f &vec);

    Vector3f Fa_b(const Vector3f &vec);

    Vector3f Fa_w(const Vector3f &vec);

    Vector3f W(const Vector3f &vec);

    Vector3f Ma_w(Vector3f &vel, Vector3f &w, Vector3f rot);

    float CD_b(const Vector3f &vec);

    float CL_w(const Vector3f &vec);

    Matrix<float, 3, 2> SFa(const Vector3f &vec, float delta);

    Matrix<float, 3, 2> SMa(const Vector3f &vec);

    Matrix<float, 3, 3> Wrot(const Vector3f &vec);

    Matrix<float, 3, 3> QuatToAtt(const Matrix<float, 4, 1> &q);

    Vector3f QuatToEuler(const Matrix<float, 4, 1> &q);

    Matrix<float, 4, 4> Omega(const Vector3f &vel);

    static Matrix<float, 4, 1> EulToQuat(const Vector3f &v);

    void my_system(const state_type &x, state_type &dxdt,double t);

    void simulate();

private:
    float mpay;
    float mpar;
    float mtot;
    float imto;
    float xb;
    float zb;
    float b;
    float c;
    float th;
    float a;
    float sw;
    float sp;
    Eigen::Matrix<float, 3, 3> inv_inertia;
    Eigen::Vector3f xgp;
    Eigen::Vector3f xgb;
    Eigen::Matrix<float, 3, 3> rgp;
    Eigen::Matrix<float, 3, 3> rgb;
    float CL0;
    float CLa;
    float CD0;
    float CDa;
    float Clp;
    float Clphi;
    float Cmq;
    float Cm0;
    float Cnr;
    float Cma;
    float CLda;
    float CLds;
    float CDda;
    float CDds;
    float Clda;
    float Cnda;
    state_type initial_state;
    float dt;
    float fin_alt;
    Vector3f wind;

};


#endif //INC_6DOF_SIMPLIFIED_PARAFOIL_H
