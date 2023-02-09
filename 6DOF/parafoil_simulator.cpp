#include <iostream>
#include <boost/numeric/odeint.hpp>
#include <Eigen/Dense>
#include <Eigen/Core>
#include "tools.h"
#include "constants.h"
#include "forces.h"
#include "joint.h"
#include "plotter.h"


using namespace std;
using namespace boost::numeric::odeint;
using namespace Eigen;
using namespace constants;

typedef Matrix<double, 12, 1> state_type;
typedef Matrix<double, 6, 6> state_matrix;


void my_system(const state_type &x, state_type &dxdt, const double t){
    Vector3d pos_c, rot_c, vel_c, vel_rot, vel_b, vel_p, vel_rot_eul;
    pos_c << x(0), x(1), x(2);
    rot_c << x(3), x(4), x(5);
    vel_c << x(6), x(7), x(8);
    vel_rot << x(9), x(10), x(11);
    
    //vel_rot_eul = Wrot(rot_c)*vel_rot;

    Vector2d sig;
    sig << 0.0, 0.0;

    Vector3d Xgb, Xgp;
    Xgp <<  -1.3,  0,  -7.4;
    Xgb <<  0.0,  0,  0.66;
    vel_b = vel_c+Omrot(vel_rot)*Xgb;
    vel_p = vel_c+Omrot(vel_rot)*Xgp;
    
    MatrixXd I = Matrix<double, 3, 3>::Identity();
    state_matrix A;
    A.setZero();
    A.block<3,3>(0,0) = ((Mpar+Mpay)*I+MF());
    A.block<3,3>(3,0) = Omrot(Xgp)*MF();
    A.block<3,3>(3,3) = Ip()+Ib()+IF();

    Matrix<double, 6, 1> B,Sd;
    B.setZero();
    B.block<3,1>(0,0) =  Fa_w(vel_p)+W_w(rot_c)+W_b(rot_c) + Fa_b(vel_b) - Omrot(vel_rot)*MF()*vel_p - ((Mpar+Mpay)*I+MF())*Omrot(vel_rot)*vel_c;
    B.block<3,1>(3,0) = Ma_w(vel_p,vel_rot,rot_c)-Omrot(vel_p)*MF()*vel_p + Omrot(Xgp)*Fa_w(vel_p) - Omrot(Xgp)*Omrot(vel_rot)*MF()*vel_p + Omrot(Xgb)*Fa_b(vel_b) - Omrot(vel_rot)*(Ip()+Ib()+IF())*vel_rot;

    Matrix<double, 6, 2> S;
    S.setZero();
    S.block<3,2>(0,0) =  SFa(vel_p,sig(0));
    S.block<3,2>(3,0) = SMa(vel_p,sig(0))+Omrot(Xgp)*SFa(vel_p,sig(0));
    //cout<<"S: "<<endl<<S*sig<<endl;
    if (t>=10 && t<30) {sig << 40.0*pi/180, 0.0*pi/180;}

    Sd = S*sig;
   
    Matrix<double, 6, 1> system_res;
    system_res = A.completeOrthogonalDecomposition().solve(B) + A.completeOrthogonalDecomposition().solve(Sd);
    
    vel_c = Trot(rot_c)*vel_c;
    dxdt.block<3,1>(0,0) = vel_c;
    dxdt.block<3,1>(3,0) = vel_rot;
    dxdt.block<6,1>(6,0) = system_res;
    
}



int main(){

state_type x;
x << 0.0, 0.0, -200.0, //position
    0.0, 0.0, 0.0*pi/180,    //rotation 
    13.0, 0.0, 2.0,  //velocity 
    0.0, 0.0, 0.0;    //vel rotation 
    


double t0 = 0.0, tf = 40.0;

    // Set the time step
double dt = 0.001;
double t = t0;

vector<double> x_c,y_c,z_c,v_x,v_y,v_z,t_s,r_x,r_y,r_z;
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
        r_x.push_back(x(3)*180/pi);
        r_y.push_back(x(4)*180/pi);
        r_z.push_back(x(5)*180/pi);
        v_x.push_back(x(6));
        v_y.push_back(x(7));
        v_z.push_back(x(8));
        t_s.push_back(t);
        z_ref = x(2);
    }
int p = plot_flight(x_c,y_c,z_c,r_x,r_y,r_z,v_x,v_y,v_z,t_s);


return 0;

}