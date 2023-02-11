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

typedef Matrix<double, 18, 1> state_type;
typedef Matrix<double, 9, 9> state_matrix;


void my_system(const state_type &x, state_type &dxdt, const double t){
    Vector3d pos_c, rotb, rotp, vel_c, vel_rotb, vel_rotp, Fr, vel_b, vel_p;
    pos_c << x(0), x(1), x(2);
    rotb << x(3), x(4), x(5);
    rotp << x(6), x(7), x(8);
    vel_c << x(9), x(10), x(11);
    vel_rotb << x(12), x(13), x(14);
    vel_rotp << x(15), x(16), x(17);
    
    vel_rotb = Wrot(rotb)*vel_rotb;
    vel_rotp = Wrot(rotp)*vel_rotp;

    Vector3d Xcb, Xcp;
    Xcp <<  -1.3,  0,  -7.5;
    Xcb <<  0.0,  0,  0.5;
    vel_b = Trot(rotb)*vel_c+Omrot(vel_rotb)*Xcb;
    vel_p = Trot(rotp)*vel_c+Omrot(vel_rotp)*Xcp;
    
    
    
    MatrixXd I = Matrix<double, 3, 3>::Identity();
    state_matrix A;
    A.setZero();
    A.block<3,3>(0,0) = (Mpar*I+MF())*Trot(rotp);
    A.block<3,3>(0,3) = -(Mpar*I+MF())*Omrot(Xcp);
    A.block<3,3>(3,3) = Ip()+IF();
    A.block<3,3>(6,0) = Mpay*Trot(rotb);
    A.block<3,3>(6,6) = -Mpay*Omrot(Xcb);
    
    Matrix<double, 9, 1> B;
    B.setZero();
    B.block<3,1>(0,0) =  Fa_w(vel_p)+W_w(rotp)-((Mpar*I+MF())*Omrot(vel_rotp)*Omrot(vel_rotp))*(Xcp) - Omrot(vel_rotp)*MF()*vel_p;    
    B.block<3,1>(3,0) = Ma_w(vel_p,vel_rotp,rotp)-Trot(rotp)*(Trot(rotb).transpose())*Mc(rotb,rotp,vel_rotb,vel_rotp)-Omrot(vel_rotp)*(Ip()+IF())*vel_rotp - Omrot(vel_p)*MF()*vel_p;
    B.block<3,1>(6,0) = Fa_b(vel_b)+W_b(rotb)-Mpay*Omrot(vel_rotb)*Omrot(vel_rotb)*Xcb;


    Matrix<double, 9, 1> system_res;
    system_res = A.completeOrthogonalDecomposition().solve(B);

    dxdt.block<3,1>(0,0) = vel_c;
    dxdt.block<3,1>(3,0) = Wrot(rotb)*vel_rotb;
    dxdt.block<3,1>(6,0) = Wrot(rotp)*vel_rotp;
    dxdt.block<9,1>(9,0) = system_res;
    
    
}



int main(){

state_type x;
x << 0.0, 0.0, -200.0, //position
    0.0, 0.0, 0.0,    //rotation b
    0.0, 0.0, 0.0,    //rotation p
    13.0, 0.0, 2.0,  //velocity c
    0.0, 0.0, 0.0,    //vel rotation b
    0.0, 0.0, 0.0;    //vel rotation p
    


double t0 = 0.0, tf = 40.0;

    // Set the time step
double dt = 0.001;
double t = t0;

vector<double> x_c,y_c,z_c,v_x,v_y,v_z,t_s;
double z_land = 0.0;
double z_ref = 200.0;

cout << "Insert final time\n"<<endl;
cin >> tf;

runge_kutta4<state_type> stepper;
    while (t <= tf) {
        t += dt;
        stepper.do_step(my_system, x, t, dt);
        x_c.push_back(x(0));
        y_c.push_back(x(1));
        z_c.push_back(x(2));
        v_x.push_back(x(9));
        v_y.push_back(x(10));
        v_z.push_back(x(11));
        t_s.push_back(t);
        z_ref = x(2);
    }

int p = plot_flight(x_c,y_c,z_c,v_x,v_y,v_z,t_s);


return 0;

}