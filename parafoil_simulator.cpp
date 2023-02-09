#include <iostream>
#include <boost/numeric/odeint.hpp>
#include <Eigen/Dense>
#include <Eigen/Core>
//#include <math.h>
#include "tools.h"
#include "constants.h"
#include "forces.h"
#include "joint.h"
#include "plotter.h"


using namespace std;
using namespace boost::numeric::odeint;
using namespace Eigen;
using namespace constants;

typedef Matrix<double, 21, 1> state_type;
typedef Matrix<double, 12, 12> state_matrix;


void my_system(const state_type &x, state_type &dxdt, const double t){
    Vector3d pos_c, rotb, rotp, vel_c, vel_rotb, vel_rotp, Fr, vel_b, vel_p;
    pos_c << x[0], x[1], x[2];
    rotb << x[3], x[4], x[5];
    rotp << x[6], x[7], x[8];
    vel_c << x[9], x[10], x[11];
    vel_rotb << x[12], x[13], x[14];
    vel_rotp << x[15], x[16], x[17];
    Fr << x[18], x[19], x[20];

    Vector3d Xcb, Xcp;
    Xcp <<  0,  0,  7.5;
    Xcb <<  0,  0,  -0.5;
    vel_b = Trot(rotb)*vel_c+Omrot(rotb)*Xcb;
    vel_p = Trot(rotp)*vel_c+Omrot(rotp)*Xcp;

    MatrixXd I = Matrix<double, 3, 3>::Identity();
    state_matrix A;
    A.setZero();
    A.block<3,3>(0,0) = (Mpar*I+MF())*Trot(rotp);
    A.block<3,3>(0,3) = -(Mpar*I+MF())*Omrot(Xcp);
    A.block<3,3>(0,9) = -Trot(rotp);
    A.block<3,3>(3,3) = Ip()+IF();
    A.block<3,3>(3,9) = Omrot(Xcp)*Trot(rotp);
    A.block<3,3>(6,0) = (Mpay*I)*Trot(rotb);
    A.block<3,3>(6,6) = -(Mpay*I)*Omrot(Xcb);
    A.block<3,3>(6,9) = Trot(rotb);
    A.block<3,3>(9,6) = Ib();
    A.block<3,3>(9,9) = -Omrot(Xcb)*Trot(rotb);
    Matrix<double, 12, 1> B;
    B.setZero();
    B.block<3,1>(0,0) =  Fa_w(vel_p)+W_w(rotp)-((Mpar*I+MF())*Omrot(vel_rotp)*Omrot(vel_rotp))*(Xcp) - Omrot(vel_rotp)*MF()*(Trot(rotp)*vel_c + Omrot(vel_rotp)*Xcp);
    B.block<3,1>(3,0) = Ma_w(vel_p,rotp,vel_rotp)-Trot(rotp)*(Trot(rotb).transpose())*Mc(rotb,rotp,vel_rotb,vel_rotp)-Omrot(vel_rotp)*(Ip()+IF())*vel_rotp - Omrot(vel_p)*MF()*vel_p;
    B.block<3,1>(6,0) = Fa_b(vel_b)+W_b(rotp)-Mpay*I*Omrot(vel_b)*Omrot(vel_b)*Xcb;
    B.block<3,1>(9,0) = Mc(rotb,rotp,vel_rotb,vel_rotp)-Omrot(vel_rotb)*Ip()*vel_rotb;

    Matrix<double, 12, 1> system_res;
    system_res = A.colPivHouseholderQr().solve(B);

    dxdt.block<3,1>(0,0) = vel_c;
    dxdt.block<3,1>(3,0) = vel_rotb;
    dxdt.block<3,1>(6,0) = vel_rotp;
    dxdt.block<12,1>(9,0) = system_res;
    
}



int main(){

state_type x;
x << 0.0, 0.0, -200.0, 
    0.0, 0.0, 0.0, 
    0.0, 0.0, 0.0, 
    10.0, 0.0, 3.0, 
    0.0, 0.0, 0.0,
    0.0, 0.0, 0.0,
    0.0, 0.0, 0.0;


double t0 = 0.0, tf = 40.0;

    // Set the time step
double dt = 0.001;
double t = t0;

vector<double> x_c,y_c,z_c,v_x,v_y,v_z,t_s;
double z_land = 0.0;
double z_ref = 200.0;

runge_kutta4<state_type> stepper;
    while (t <= 0.02) {
        t += dt;
        stepper.do_step(my_system, x, t, dt);
        x_c.push_back(x(0));
        y_c.push_back(x(1));
        z_c.push_back(x(2));
        v_x.push_back(x(10));
        v_y.push_back(x(11));
        v_z.push_back(x(12));
        t_s.push_back(t);
        z_ref = x(2);
    }

int p = plot_flight(x_c,y_c,z_c,v_x,v_y,v_z,t_s);


return 0;

}