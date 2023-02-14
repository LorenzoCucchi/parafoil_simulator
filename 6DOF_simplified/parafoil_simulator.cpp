#include <iostream>
#include <boost\numeric\odeint.hpp>
#include <C:\Program Files\eigen-3.4.0\Eigen\Dense>
#include <chrono>
#include <fstream>
#include "tools.h"
#include "constants.h"
#include "forces.h"
#include "wind_estimation.h"


using namespace std;
using namespace std::chrono;
using namespace boost::numeric::odeint;
using namespace Eigen;
using namespace constants;

typedef Matrix<double, 13, 1> state_type;

void my_system(const state_type &x, state_type &dxdt, const double t) {

    Vector3d pos_c {x(0), x(1), x(2)};
    Matrix<double, 4,1> quat {x(3), x(4), x(5), x(6)};
    Vector3d vel_c {x(7), x(8), x(9)};
    Vector3d vel_rot {x(10), x(11), x(12)};

    Vector3d euler = QuatToEuler(quat);
    Vector3d vel_rot_eul = Wrot(euler) * vel_rot;


    Matrix3d T = QuatToAtt(quat);
    Matrix3d Sk_Om = Omrot(vel_rot);
    Vector3d Vw = vel_c - T * Wind();
    Vector3d vel_b = Vw + Sk_Om * Xgb;
    Vector3d vel_p = Vw + Sk_Om * Xgp;
    Vector2d sig {0.0, 0.0};

    if (t >= 10 && t < 1000) { sig << 40.0 * pi / 180, 0.0 * pi / 180; }
    if (sig(0) > sig(1)) { sig(0) = sig(0) - sig(1); }
    else {
        sig(1) = sig(0);
        sig(0) = sig(0) - sig(1);
    }

    Vector3d Fa = Fa_w(vel_p);
    Vector3d Fb = Fa_b(vel_b);
    Matrix<double, 3, 2> Sfa = SFa(vel_p, sig(0));
    dxdt.block<3, 1>(0, 0) = T.transpose() * vel_c;
    dxdt.block<4, 1>(3, 0) = 0.5 * Omega(vel_rot) * quat.normalized();
    dxdt.block<3, 1>(7, 0) = IMtot * (Fa + Fb + W_w(euler) + W_b(euler) -
                                                  Mtot * vel_rot.cross(vel_c) + Sfa * sig);
    dxdt.block<3, 1>(10, 0) = I_i * (Ma_w(vel_p, vel_rot_eul, euler) + Rgp * Fa + Rgb * Fb -
                                     Sk_Om * I * vel_rot + (SMa(vel_p) + Rgp * Sfa) * sig);
  
}

int main() {

    Vector3d euler;
    euler << 0.0, 0.0, 0.0;

    state_type x;

    x.block<3, 1>(0, 0) << 0, 0, -400;
    x.block<4, 1>(3, 0) = EulToQuat(euler);
    x.block<3, 1>(7, 0) << 13, 0, 2;
    x.block<3, 1>(10, 0) << 0, 0, 0;

    double t0 = 0.0, tf = 40.0;

    // Set the time step
    double dt = 0.01;
    double t = t0;

    vector<double> x_c, y_c, z_c, v_x, v_y, v_z, t_s, r_x, r_y, r_z;
    double z_land = 0.0;
    double z_ref = x(2);


    Matrix<double, 4, 1> quat;
    Vector3d eul,vel;

    auto start = high_resolution_clock::now();
    runge_kutta4<state_type> stepper;

    while (z_ref <= z_land) {
        t += dt;
        stepper.do_step(my_system, x, t, dt);
        x_c.push_back(x(0));
        y_c.push_back(x(1));
        z_c.push_back(-x(2));
        quat << x(3), x(4), x(5), x(6);
        eul = QuatToEuler(quat);
        r_x.push_back(eul(0) * 180 / pi);
        r_y.push_back(eul(1) * 180 / pi);
        r_z.push_back(eul(2) * 180 / pi);
        vel << x(7),x(8),x(9);
        vel = QuatToAtt(quat).transpose()*vel;
        v_x.push_back(vel(0));
        v_y.push_back(vel(1));
        v_z.push_back(vel(2));
        t_s.push_back(t);
        z_ref = x(2);
    }


    auto stop = high_resolution_clock::now();  
    auto duration = duration_cast<microseconds>(stop - start);
    cout << duration.count() << endl;

    WindEstimation(t_s, v_x, v_y);

    //fstream file;
    //file.open("simulation_results.txt", ios_base::out );
    //for (int i = 0; i<t_s.size(); i++)
    //{
    //    file<< t_s[i] << " " << x_c[i] << " " << y_c[i] << " " << z_c[i] << " " << r_x[i] << " " << r_y[i] << " " << r_z[i] << " " << v_x[i] << " " << v_y[i] << " " << v_z[i] << endl;
    //}
    //file.close();

//
//

    return 0;

}