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

typedef Matrix<double, 20, 1> state_type;
typedef Matrix<double, 9, 9> state_matrix;


void my_system(const state_type &x, state_type &dxdt, const double t){
    Vector3d pos_c, vel_c, vel_rotb, vel_rotp, Fr, vel_b, vel_p, euler_b,euler_p, vel_rotb_eul, vel_rotp_eul;
    Matrix<double, 4, 1> quat_b, quat_p;
    pos_c << x(0), x(1), x(2);
    quat_b << x(3), x(4), x(5), x(6);
    quat_p << x(7), x(8), x(9), x(10);
    vel_c << x(10), x(11), x(12);
    vel_rotb << x(13), x(14), x(15);
    vel_rotp << x(16), x(17), x(18);
    

    euler_b = QuatToEuler(quat_b);
    euler_p = QuatToEuler(quat_p);

    vel_rotb_eul = Wrot(euler_b)*vel_rotb;
    vel_rotp_eul = Wrot(euler_p)*vel_rotp;

    Matrix3d Sk_Om_b, Sk_Om_p, T_b, T_p;
    T_b = QuatToAtt(quat_b);
    T_p = QuatToAtt(quat_p);
    Sk_Om_b = Omrot(vel_rotb);
    Sk_Om_p = Omrot(vel_rotp);

    vel_b = T_b*vel_c+Sk_Om_b*Xgb;
    vel_p = T_p*vel_c+Sk_Om_p*Xgp;  
    
    state_matrix A;
    A.setZero();
    A.block<3,3>(0,0) = (Mpar*eye+Mf)*T_p;
    A.block<3,3>(0,3) = -(Mpar*eye+Mf)*Rgp;
    A.block<3,3>(3,3) = I_p+If;
    A.block<3,3>(6,0) = Mpay*T_b;
    A.block<3,3>(6,6) = -Mpay*Rgb;
    
    Matrix<double, 9, 1> B;
    B.setZero();
    B.block<3,1>(0,0) =  Fa_w(vel_p)+W_w(euler_p)-((Mpar*I+Mf)*Sk_Om_p*Sk_Om_p)*(Xgp) - Sk_Om_p*Mf*vel_p;    
    B.block<3,1>(3,0) = Ma_w(vel_p,vel_rotp_eul,euler_p)-T_p*(T_b.transpose())*Mc(euler_b,euler_p,vel_rotb_eul,vel_rotp_eul)-Sk_Om_p*(I_p+If)*vel_rotp - Omrot(vel_p)*Mf*vel_p;
    B.block<3,1>(6,0) = Fa_b(vel_b)+W_b(euler_b)-Mpay*Sk_Om_b*Sk_Om_b*Xgb;

    Matrix<double, 9, 1> system_res;
    system_res = A.completeOrthogonalDecomposition().solve(B);

    dxdt.block<3,1>(0,0) = vel_c;
    dxdt.block<4,1>(3,0) = 0.5*Omega(vel_rotb)*quat_b.normalized();
    dxdt.block<4,1>(7,0) = 0.5*Omega(vel_rotp)*quat_p.normalized();
    dxdt.block<9,1>(9,0) = system_res;
    
    
}



int main(){


Vector3d eulerp, eulerb;
eulerp << 0.0,0.0,0.0;
eulerb << 0.0,0.0,0.0;

state_type x;
    
x.block<3,1>(0,0) << 0, 0, -400;
x.block<4,1>(3,0) = EulToQuat(eulerb);
x.block<4,1>(7,0) = EulToQuat(eulerp);
x.block<3,1>(11,0) << 13, 0, 2;
x.block<6,1>(14,0) << 0, 0, 0, 0, 0, 0;


double t0 = 0.0, tf = 40.0;

    // Set the time step
double dt = 0.0001;
double t = t0;

vector<double> x_c,y_c,z_c,v_x,v_y,v_z,t_s,r_x,r_y,r_z;
double z_land = 0.0;
double z_ref = x(2);

cout << "Insert final time\n"<<endl;
cin >> tf;

Matrix<double, 4,1> quatb,quatp;
Vector3d eulp,eulb;

auto start = high_resolution_clock::now();
runge_kutta4<state_type> stepper;
    while (t <= tf) {
        t += dt;
        stepper.do_step(my_system, x, t, dt);
        x_c.push_back(x(0));
        y_c.push_back(x(1));
        z_c.push_back(-x(2));
        quatb << x(3),x(4),x(5),x(6);
        quatp << x(7),x(8),x(9),x(10);
        eulb = QuatToEuler(quatb);
        eulp = QuatToEuler(quatp);
        r_x.push_back(eulb(0)*180/pi);
        r_y.push_back(eulb(1)*180/pi);
        r_z.push_back(eulb(2)*180/pi);
        v_x.push_back(x(11));
        v_y.push_back(x(12));
        v_z.push_back(x(13));
        t_s.push_back(t);
        z_ref = x(2);
    }

auto stop = high_resolution_clock::now();
auto duration = duration_cast<milliseconds>(stop - start);
cout << duration.count() << endl;

int p = plot_flight(x_c,y_c,z_c,r_x,r_y,r_z,v_x,v_y,v_z,t_s);

return 0;

}