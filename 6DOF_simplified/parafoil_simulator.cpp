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

typedef Matrix<double, 12, 1> state_type;

void my_system(const state_type &x, state_type &dxdt, const double t){
    Vector3d pos_c, rot_c, vel_c, vel_rot, vel_b, vel_p, vel_rot_eul;
    pos_c << x(0), x(1), x(2);
    rot_c << x(3), x(4), x(5);
    vel_c << x(6), x(7), x(8);
    vel_rot << x(9), x(10), x(11);
    
    vel_rot = Wrot(rot_c)*vel_rot;

    Vector2d sig;
    sig.setZero();
 
    Matrix3d Om, Rgp, Rgb;
    Om = Omrot(vel_rot);
    Rgp = Omrot(Xgp);
    Rgb = Omrot(Xgb);
    
   
    vel_b = vel_c+Om*Xgb;
    vel_p = vel_c+Om*Xgp;
    
    if (t>=10 && t<1000) {sig << 30.0*pi/180, 0.0*pi/180;}


    Matrix<double, 6, 1> B,Sd;
    B.block<3,1>(0,0) =  Fa_w(vel_p)+W_w(rot_c)+W_b(rot_c)-(Mpay+Mpar)*Om*vel_c + Fa_b(vel_b) + SFa(vel_p,sig(0))*sig;
    B.block<3,1>(3,0) = Ma_w(vel_p,vel_rot,rot_c) + Rgp*Fa_w(vel_p) + Rgb*Fa_b(vel_b) - Om*(Ip()+Ib())*vel_rot+(SMa(vel_p,sig(0))+Rgp*SFa(vel_p,sig(0)))*sig;

      
    vel_c = Trot(rot_c).transpose()*vel_c;


    dxdt.block<3,1>(0,0) = vel_c;
    dxdt.block<3,1>(3,0) = Wrot(rot_c)*vel_rot;
    dxdt.block<3,1>(6,0) = 1/(Mpar+Mpay) * B.block<3,1>(0,0);
    dxdt.block<3,1>(9,0) = I_i*B.block<3,1>(3,0);

}



int main(){

state_type x;
x << 0.0, 0.0, -400.0, //position
    0.0, 0.0, 0.0*pi/180,    //rotation 
    13.0, 0.0, 2.0,  //velocity 
    0.0, 0.0, 0.0;    //vel rotation 
    


double t0 = 0.0, tf = 40.0;

    // Set the time step
double dt = 0.01;
double t = t0;

vector<double> x_c,y_c,z_c,v_x,v_y,v_z,t_s,r_x,r_y,r_z;
double z_land = 0.0;
double z_ref = x(2);

cout << "Insert final time\n"<<endl;
cin >> tf;

auto start = high_resolution_clock::now();
runge_kutta4<state_type> stepper;
    while (z_ref<=z_land) {
        t += dt;
        stepper.do_step(my_system, x, t, dt);
        x_c.push_back(x(0));
        y_c.push_back(x(1));
        z_c.push_back(-x(2));
        r_x.push_back(x(3)*180/pi);
        r_y.push_back(x(4)*180/pi);
        r_z.push_back(x(5)*180/pi);
        v_x.push_back(x(6));
        v_y.push_back(x(7));
        v_z.push_back(x(8));
        t_s.push_back(t);
        z_ref = x(2);
    }
auto stop = high_resolution_clock::now();
auto duration = duration_cast<milliseconds>(stop - start);
cout << duration.count() << endl;

int p = plot_flight(x_c,y_c,z_c,r_x,r_y,r_z,v_x,v_y,v_z,t_s);


return 0;

}