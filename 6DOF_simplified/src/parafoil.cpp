#include <boost\numeric\odeint.hpp>
#include <C:\Program Files\eigen-3.4.0\Eigen\Dense>
#include "constants.h"
#include "Parafoil.h"
#include "json.hpp"
#include "wind_estimation.h"
#include <fstream>
#include <functional>
#include "pid.h"

using namespace std;
using namespace boost::numeric::odeint;
using namespace Eigen;
using namespace constants;
using json = nlohmann::json;

typedef Matrix<double, 15, 1> state_type;

Parafoil::Parafoil(const string& path){
    std::ifstream f(path);
    json data = json::parse(f);
    json pay = data["Payload"];
    json par = data["Parafoil"];
    json aer = data["Aerodynamic Data"];
    json sim = data["Simulation Setup"];

    //Payload
    Parafoil::mpay = pay["Mass"];
    Parafoil::xb   = pay["Diameter"];
    Parafoil::zb   = pay["Length"];
    Parafoil::sp   = pay["Area"];
    double distb   = pay["Distance"];
    //Parafoil
    Parafoil::mpar = par["Mass"];
    Parafoil::b    = par["Span"];
    Parafoil::c    = par["Chord"];
    Parafoil::th   = par["Thickness"];
    Parafoil::a    = par["Height"];
    double distp = par["Distance"];
    double cant = par["Cant"];
    //Calculation to setup data
    Parafoil::sw   = Parafoil::b * Parafoil::c;
    Parafoil::mtot = Parafoil::mpay + Parafoil::mpar;
    Parafoil::imto = 1/Parafoil::mtot;
    Parafoil::xgp << -distp*sin(cant*pi/180),0.0,-distp*cos(cant*pi/180);
    Parafoil::xgb << distb*sin(cant*pi/180),0.0,distb*cos(cant*pi/180);
    Parafoil::inv_inertia = Inertia_calc();
    Parafoil::rgp = Omrot(Parafoil::xgp);
    Parafoil::rgb = Omrot(Parafoil::xgb);
    //Aerodynamics Data
    Parafoil::CL0   = aer["CL0"];
    Parafoil::CLa   = aer["CLa"];
    Parafoil::CD0   = aer["CD0"];
    Parafoil::CDa   = aer["CDa"];
    Parafoil::Clp   = aer["Clp"];
    Parafoil::Clphi = aer["Clphi"];
    Parafoil::Cmq   = aer["Cmq"];
    Parafoil::Cm0   = aer["Cm0"];
    Parafoil::Cnr   = aer["Cnr"];
    Parafoil::Cma   = aer["Cma"];
    Parafoil::CLda  = aer["CLda"];
    Parafoil::CLds  = aer["CLds"];
    Parafoil::CDda  = aer["CDda"];
    Parafoil::CDds  = aer["CDds"];
    Parafoil::Clda  = aer["Clda"];
    Parafoil::Cnda  = aer["Cnda"];
    //Simulation setup
    Parafoil::boolControl = sim["Control"];
    Parafoil::dt = sim["Step Time"];
    Parafoil::fin_alt = sim["Final Altitude"];
    state_type x;
    Vector3d eul;
    eul << sim["Attitude"].at(0), sim["Attitude"].at(1),sim["Attitude"].at(2);
    x.block<3, 1>(0, 0) <<  sim["Position"].at(0), sim["Position"].at(1),sim["Position"].at(2) ;
    x.block<4, 1>(3, 0) = EulToQuat(eul);
    x.block<3, 1>(7, 0) << sim["Speed"].at(0), sim["Speed"].at(1),sim["Speed"].at(2);
    x.block<3, 1>(10, 0) << sim["Angular Velocities"].at(0), sim["Angular Velocities"].at(1),sim["Angular Velocities"].at(2);
    x.block<2, 1>(13, 0).setZero();
    Parafoil::initial_state =  x;
    Parafoil::wind << sim["Wind Speed"].at(0), sim["Wind Speed"].at(1),sim["Wind Speed"].at(2);
    Parafoil::target << sim["Target"].at(0), sim["Target"].at(1);
    Parafoil::Kp = sim["Kp"];
    Parafoil::Ki = sim["Ki"];
}


void Parafoil::print(){
    cout<<mpay<<endl<<mpar<<endl<<mtot<<endl;
    cout<<imto<<endl<<xb<<endl<<zb<<endl;
    cout<<b<<endl<<c<<endl<<th<<endl;
    cout<<a<<endl<<sw<<endl<<sp<<endl;
    cout<<initial_state<<endl;
}


Matrix<double, 3, 3> Parafoil::Inertia_calc() const{
    Matrix3d in;
    in.setZero();
    in(0,0) = (pow(b,2)+pow(th,2))*mpar +(pow(xb,2)+pow(zb,2))*mpay;
    in(1,1) = (pow(c,2)+pow(th,2))*mpar +(pow(zb,2)+pow(xb,2))*mpay;
    in(2,2) = (pow(b,2)+pow(c,2))*mpar  +(2*pow(xb,2))*mpay;

    return (in/12).inverse();
}


Matrix3d Parafoil::Omrot(const Vector3d &vec){
    Matrix3d Mat;
    Mat << 0.0, -vec(2), vec(1), vec(2), 0.0, -vec(0), -vec(1), vec(0), 0.0;
    return Mat;
}


Vector3d Parafoil::Fa_b(const Vector3d &vec) {
    Vector3d f;
    f = -0.5 * rho * sp * vec.norm() * CD_b(vec) * vec;
    return f;
}


Vector3d Parafoil::Fa_w(const Vector3d &vec) {
    double q = 0.5 * rho * sw * vec.norm();
    Vector3d f;
    f(0) = q * (CL_w(vec) * vec(2) - CD_b(vec) * vec(0));
    f(1) = -q * CD_b(vec) * vec(1);
    f(2) = q * (-CL_w(vec) * vec(0) - CD_b(vec) * vec(2));
    return f;
}

Vector3d Parafoil::W(const Vector3d &vec) {
    Vector3d W;
    W(0) = g * (-sin(vec(1)));
    W(1) = g * (sin(vec(0)) * cos(vec(1)));
    W(2) = g * (cos(vec(0)) * cos(vec(1)));
    return W;
}

Vector3d Parafoil::Ma_w(Vector3d &vel, Vector3d &w, Vector3d rot) {
    double q = 0.5*rho*sw*vel.squaredNorm();
    double vnorm = vel.norm();
    Vector3d m;
    m(0) = q*(Clp*pow(b,2)*w(0)/(2*vnorm) + Clphi*b*rot(0));
    m(1) = q*(Cmq*pow(c,2)*w(1)/(2*vnorm) + Cm0*c + Cma*c*atan2(vel(2), vel(0)));
    m(2) = q*(Cnr*pow(b,2)*w(2)/(2*vnorm));
    return m;
}

double Parafoil::CD_b(const Vector3d &vec)
{
    double alpha = atan2(vec(2), vec(0));
    double cd = CD0 + CDa * alpha*alpha;
    return cd;
}


double Parafoil::CL_w(const Vector3d &vec)
{
    double alpha = atan2(vec(2), vec(0));
    double cl = CL0 + CLa * alpha;
    return cl;
}


Matrix<double, 3, 2> Parafoil::SFa(const Vector3d &vec, double delta){

    double sign;
    if (delta >= 0.0)  sign = 1;
    else sign = -1;

    double q = 0.5*rho*sw*vec.norm();
    Matrix<double, 3, 2> F;
    F(0,0) = q*(CLda*vec(2)-CDda*vec(0))*sign;
    F(0,1) = q*(CLds*vec(2)-CDds*vec(0));
    F(1,0) = q*(-CDda*vec(1)*sign);
    F(1,1) = q*(-CDds*vec(1));
    F(2,0) = q*(-CLda*vec(0)-CDda*vec(2))*sign;
    F(2,1) = q*(-CLds*vec(0)-CDds*vec(2));
    return F;
}

Matrix<double, 3, 2> Parafoil::SMa(const Vector3d &vec){

    double q = 0.5*rho*sw*vec.squaredNorm()*b/th;
    Matrix<double, 3, 2> M;
    M(0,0) = q*Clda;
    M(0,1) = 0.0;
    M(1,0) = 0.0;
    M(1,1) = 0.0;
    M(2,0) = q*Cnda;
    M(2,1) = 0.0;

    return M;
}


Matrix<double, 3, 3> Parafoil::Wrot(const Vector3d &vec){
    Matrix3d Mat;
    Mat.row(0) << 1.0, sin(vec(0))*tan(vec(1)), cos(vec(0))*tan(vec(1));
    Mat.row(1) << 0.0, cos(vec(0)), -sin(vec(0));
    Mat.row(2) << 0.0, sin(vec(0))/cos(vec(1)), cos(vec(0))/cos(vec(1));
    return Mat;
}


Matrix<double, 3, 3> Parafoil::QuatToAtt(const Matrix<double, 4, 1> &q){
    Matrix<double,4 ,1 > quat = q.normalized();
    Vector3d qv;
    qv = quat.block<3, 1>(0, 0);
    Matrix<double, 3, 3> rox;
    Matrix<double,3,3> eye = Matrix<double, 3, 3>::Identity();
    rox << 0, -qv(2), qv(1),
            qv(2), 0, -qv(0),
            -qv(1), qv(0), 0;

    return (pow(quat(3),2)-qv.squaredNorm())*eye + 2*(qv*qv.transpose())-2*quat(3)*rox;
}


Vector3d Parafoil::QuatToEuler(const Matrix<double, 4, 1> &q){

    Matrix<double,4 ,1 > quat = q.normalized();
    Vector3d euler;
    euler(0) = atan2((2*quat(0)*quat(3)+2*quat(1)*quat(2)),(-pow(quat(0),2)-pow(quat(1),2)+pow(quat(2),2)+pow(quat(3),2)));
    euler(1) = -asin(2*quat(0)*quat(2)-2*quat(1)*quat(3));
    euler(2) = atan2((2*quat(0)*quat(1)+2*quat(2)*quat(3)),(pow(quat(0),2)-pow(quat(1),2)-pow(quat(2),2)+pow(quat(3),2)));

    return euler;
}

Matrix<double, 4, 4> Parafoil::Omega(const Vector3d &vel){
    Matrix<double, 4,4> Om;
    Om <<          0,  vel(2), -vel(1), vel(0),
            -vel(2),           0,  vel(0), vel(1),
            vel(1), -vel(0),           0, vel(2),
            -vel(0), -vel(1), -vel(2),          0;
    return Om;
}

Matrix<double, 4, 1> Parafoil::EulToQuat(const Vector3d &v){
    Matrix<double, 4, 1> quat;
    Vector3d vec = v/2;
    double c0 = cos(vec(0)), s0 = sin(vec(0));
    double c1 = cos(vec(1)), s1 = sin(vec(1));
    double c2 = cos(vec(2)), s2 = sin(vec(2));

    quat(3) = c0*c1*c2-s0*s1*s2;
    quat(0) = s0*c1*c2+c0*s1*s2;
    quat(1) = c0*s1*c2+s0*c1*s2;
    quat(2) = c0*c1*s2-s0*s1*c2;
    return quat;
}



void Parafoil::my_system(const state_type &x, state_type &dxdt, const double t) {

    Vector3d pos_c {x(0), x(1), x(2)};
    Matrix<double, 4,1> quat {x(3), x(4), x(5), x(6)};
    Vector3d vel_c {x(7), x(8), x(9)};
    Vector3d vel_rot {x(10), x(11), x(12)};
    Vector3d euler = QuatToEuler(quat);
    Vector3d vel_rot_eul = Wrot(euler) * vel_rot;


    Matrix3d T = QuatToAtt(quat);
    Matrix3d Sk_Om = Omrot(vel_rot);
    Vector3d Vw = vel_c - T * wind;
    Vector3d vel_b = Vw + Sk_Om * xgb;
    Vector3d vel_p = Vw + Sk_Om * xgp;
    Vector2d sig {x(13), x(14)};

    sig = sig*2;

    Vector3d Fa = Fa_w(vel_p);
    Vector3d Fb = Fa_b(vel_b);
    Matrix<double, 3, 2> Sfa = SFa(vel_p, sig(0));
    dxdt.block<3, 1>(0, 0) = T.transpose() * vel_c;
    dxdt.block<4, 1>(3, 0) = 0.5 * Omega(vel_rot) * quat.normalized();
    dxdt.block<3, 1>(7, 0) = imto * (Fa + Fb + W(euler)*(mtot) -
                                      mtot * vel_rot.cross(vel_c) + Sfa * sig);
    dxdt.block<3, 1>(10, 0) = inv_inertia * (Ma_w(vel_p, vel_rot_eul, euler) + rgp * Fa + rgb * Fb -
                                     Sk_Om * (inv_inertia) * vel_rot + (SMa(vel_p) + rgp * Sfa) * sig);
    dxdt(13) = 0.0;
    dxdt(14) = 0.0;

}

void Parafoil::simulate(){
    double t0 = 0.0;
    double t = t0;

    vector<double> x_c, y_c, z_c, v_x, v_y, v_z, t_s, r_x, r_y, r_z;
    double z_land = 0.0;
    double z_ref = initial_state(2);

    state_type x = initial_state;
    Matrix<double, 4, 1> quat;
    Vector3d eul,vel;

    auto sys = [&capture0 = *this](auto && PH1, auto && PH2, auto && PH3) {
        capture0.my_system(std::forward<decltype(PH1)>(PH1),std::forward<decltype(PH2)>(PH2), std::forward<decltype(PH3)>(PH3));
    };

    auto start = std::chrono::high_resolution_clock::now();
    runge_kutta4<state_type> stepper;
       while (z_ref <= fin_alt) {
        t += dt;
        stepper.do_step(sys , x, t, dt);
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


      auto stop = std::chrono::high_resolution_clock::now();
      auto duration = chrono::duration_cast<std::chrono::microseconds>(stop - start);
      cout << duration.count() << endl;

    WindEstimation(t_s, v_x, v_y);

    fstream file;
    file.open("simulation_results.txt", ios_base::out );
    for (int i = 0; i<t_s.size(); i++)
    {
        file<< t_s[i] << " " << x_c[i] << " " << y_c[i] << " " << z_c[i] << " " << r_x[i] << " " << r_y[i] << " " << r_z[i] << " " << v_x[i] << " " << v_y[i] << " " << v_z[i] << endl;
    }
     file.close();

}


void Parafoil::simulate_control(){
    double t0 = 0.0;
    double t = t0;

    vector<double> x_c, y_c, z_c, v_x, v_y, v_z, t_s, r_x, r_y, r_z, u_s, e_s;
    double z_land = 0.0;
    double z_ref = initial_state(2);
    double errore = 0.0;
    double u = 0.0;
    state_type x=initial_state;
    Matrix<double, 4, 1> quat;
    Vector3d eul,vel;

    auto sys = [&capture0 = *this](auto && PH1, auto && PH2, auto && PH3) {
        capture0.my_system(std::forward<decltype(PH1)>(PH1),std::forward<decltype(PH2)>(PH2), std::forward<decltype(PH3)>(PH3));
    };

    Pid pid;

    auto start = std::chrono::high_resolution_clock::now();
    runge_kutta4<state_type> stepper;
    while (z_ref <= fin_alt) {

        double tin = 0.0;
        while(tin <= pid.dt){
            tin += dt;
            stepper.do_step(sys , x, tin, dt);
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
            t_s.push_back(tin+t);
            u_s.push_back(x(13));
            e_s.push_back(errore);
        }
        t += pid.dt;
        errore = pid.GuidanceSys(target, x.block<2,1>(0,0), eul(2));
        u = pid.PI(errore,Kp,Ki);
        u = pid.Saturation(u);

        x(13) = u;
        x(14) = 0;

        z_ref = x(2);
    }

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<std::chrono::microseconds>(stop - start);
    cout << duration.count() << endl;

    fstream file;
    file.open("simulation_results.txt", ios_base::out );
    for (int i = 0; i<t_s.size(); i++)
    {
        file<< t_s[i] << " " << x_c[i] << " " << y_c[i] << " " << z_c[i] << " " << r_x[i] << " " << r_y[i] << " " << r_z[i] << " " << v_x[i] << " " << v_y[i] << " " << v_z[i] <<" "<< u_s[i]<< endl;
    }
    file.close();
}