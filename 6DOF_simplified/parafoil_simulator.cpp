#include <iostream>
#include <boost/numeric/odeint.hpp>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <chrono>
#include "tools.h"
#include "constants.h"
#include "forces.h"
#include "joint.h"
#include "plotter.h"


using namespace std;
using namespace std::chrono;
using namespace boost::numeric::odeint;
using namespace Eigen;
using namespace constants;

typedef Matrix<double, 13, 1> state_type;

void my_system(const state_type &x, state_type &dxdt, const double t){
    Vector3d pos_c, vel_c, vel_rot, vel_rot_eul, vel_b, vel_p, euler, Vw;
    Matrix<double, 4, 1> quat;
    pos_c << x(0), x(1), x(2);
    quat << x(3), x(4), x(5), x(6);
    vel_c << x(7), x(8), x(9);
    vel_rot << x(10), x(11), x(12);
    
    euler = QuatToEuler(quat);
    vel_rot_eul = Wrot(euler)*vel_rot;

    Matrix3d Sk_Om, T;
    T = QuatToAtt(quat);
    Sk_Om = Omrot(vel_rot);
    Vw = vel_c - T*Wind();
    vel_b = Vw + Sk_Om*Xgb;
    vel_p = Vw + Sk_Om*Xgp;
    Vector2d sig;
    sig.setZero();
    if (t>=10 && t<1000) {sig << 30.0*pi/180, 0.0*pi/180;}

    Matrix<double, 6, 1> B,Sd;
    B.block<3,1>(0,0) =  Fa_w(vel_p)+ Fa_b(vel_b) + W_w(euler)+W_b(euler) - (Mpay+Mpar)*vel_rot.cross(vel_c)  + SFa(vel_p,sig(0))*sig;
    B.block<3,1>(3,0) = Ma_w(vel_p,vel_rot_eul,euler) + Rgp*Fa_w(vel_p) + Rgb*Fa_b(vel_b) - Sk_Om*I*vel_rot+(SMa(vel_p,sig(0))+Rgp*SFa(vel_p,sig(0)))*sig;

    dxdt.block<3,1>(0,0) = T.transpose()*vel_c;
    dxdt.block<4,1>(3,0) = 0.5*Omega(vel_rot)*quat.normalized();
    dxdt.block<3,1>(7,0) = 1/(Mpar+Mpay) * B.block<3,1>(0,0);
    dxdt.block<3,1>(10,0) = I_i*B.block<3,1>(3,0);

}



int main(){

Vector3d euler;
euler << 0.0,0.0,0.0;


state_type x;
    
x.block<3,1>(0,0) << 0, 0, -400;
x.block<4,1>(3,0) = EulToQuat(euler);
x.block<3,1>(7,0) << 13, 0, 2;
x.block<3,1>(10,0) << 0, 0, 0;

double t0 = 0.0, tf = 40.0;

    // Set the time step
double dt = 0.01;
double t = t0;

vector<double> x_c,y_c,z_c,v_x,v_y,v_z,t_s,r_x,r_y,r_z;
double z_land = 0.0;
double z_ref = x(2);

cout << "Insert final time\n"<<endl;
cin >> tf;

Matrix<double, 4,1> quat;
Vector3d eul;

auto start = high_resolution_clock::now();
runge_kutta4<state_type> stepper;
    while (z_ref<=z_land) {
        t += dt;
        stepper.do_step(my_system, x, t, dt);
        x_c.push_back(x(0));
        y_c.push_back(x(1));
        z_c.push_back(-x(2));
        quat << x(3),x(4),x(5),x(6);
        eul = QuatToEuler(quat);
        r_x.push_back(eul(0)*180/pi);
        r_y.push_back(eul(1)*180/pi);
        r_z.push_back(eul(2)*180/pi);
        v_x.push_back(x(7));
        v_y.push_back(x(8));
        v_z.push_back(x(9));
        t_s.push_back(t);
        z_ref = x(2);
    }
auto stop = high_resolution_clock::now();
auto duration = duration_cast<milliseconds>(stop - start);
cout << duration.count() << endl;

int p = plot_flight(x_c,y_c,z_c,r_x,r_y,r_z,v_x,v_y,v_z,t_s);


return 0;

}