#include <boost\numeric\odeint.hpp>
#include <C:\Program Files\eigen-3.4.0\Eigen\Dense>
#include "constants.h"
#include "Parafoil.h"
#include "json.hpp"
#include "wind_estimation.h"
#include <fstream>
#include <functional>

using namespace std;
using namespace boost::numeric::odeint;
using namespace Eigen;
using namespace constants;
using json = nlohmann::json;

typedef Matrix<float, 13, 1> state_type;

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
    float distb   = pay["Distance"];
    //Parafoil
    Parafoil::mpar = par["Mass"];
    Parafoil::b    = par["Span"];
    Parafoil::c    = par["Chord"];
    Parafoil::th   = par["Thickness"];
    Parafoil::a    = par["Height"];
    float distp = par["Distance"];
    float cant = par["Cant"];
    //Calculation to setup data
    Parafoil::sw   = Parafoil::b * Parafoil::c;
    Parafoil::mtot = Parafoil::mpay + Parafoil::mpar;
    Parafoil::imto = 1/Parafoil::mtot;
    Parafoil::xgp << -distp*sin(cant*pi/180),0.0,-distp*cos(cant*pi/180);
    Parafoil::xgb << distb*sin(cant*pi/180),0.0,distb*cos(cant*pi/180);
    Parafoil::inv_inertia = Inertia_calc();
    Parafoil::rgp = Omrot(Parafoil::xgp);
    Parafoil::rgb = Omrot(Parafoil::xgb);
    cout<<xgp<<endl<<xgb<<endl;
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
    Parafoil::dt = sim["Step Time"];
    Parafoil::fin_alt = sim["Final Altitude"];
    state_type x;
    Vector3f eul;
    eul << sim["Attitude"].at(0), sim["Attitude"].at(1),sim["Attitude"].at(2);
    x.block<3, 1>(0, 0) <<  sim["Position"].at(0), sim["Position"].at(1),sim["Position"].at(2) ;
    x.block<4, 1>(3, 0) = EulToQuat(eul);
    x.block<3, 1>(7, 0) << sim["Speed"].at(0), sim["Speed"].at(1),sim["Speed"].at(2);
    x.block<3, 1>(10, 0) << sim["Angular Velocities"].at(0), sim["Angular Velocities"].at(1),sim["Angular Velocities"].at(2);
    Parafoil::initial_state =  x;
    Parafoil::wind << sim["Wind Speed"].at(0), sim["Wind Speed"].at(1),sim["Wind Speed"].at(2);
}


void Parafoil::print(){
    cout<<mpay<<endl<<mpar<<endl<<mtot<<endl;
    cout<<imto<<endl<<xb<<endl<<zb<<endl;
    cout<<b<<endl<<c<<endl<<th<<endl;
    cout<<a<<endl<<sw<<endl<<sp<<endl;
    cout<<initial_state<<endl;
}


Matrix<float, 3, 3> Parafoil::Inertia_calc(){
    Matrix3f in;
    in.setZero();
    in(0,0) = (pow(b,2)+pow(th,2))*mpar +(pow(xb,2)+pow(zb,2))*mpay;
    in(1,1) = (pow(c,2)+pow(th,2))*mpar +(pow(zb,2)+pow(xb,2))*mpay;
    in(2,2) = (pow(b,2)+pow(c,2))*mpar  +(2*pow(xb,2))*mpay;

    return (in/12).inverse();
}


Matrix3f Parafoil::Omrot(const Vector3f &vec){
    Matrix3f Mat;
    Mat << 0.0, -vec(2), vec(1), vec(2), 0.0, -vec(0), -vec(1), vec(0), 0.0;
    return Mat;
}


Vector3f Parafoil::Fa_b(const Vector3f &vec) {
    Vector3f f;
    f = -0.5 * rho * sp * vec.norm() * CD_b(vec) * vec;
    return f;
}


Vector3f Parafoil::Fa_w(const Vector3f &vec) {
    float q = 0.5 * rho * sw * vec.norm();
    Vector3f f;
    f(0) = q * (CL_w(vec) * vec(2) - CD_b(vec) * vec(0));
    f(1) = -q * CD_b(vec) * vec(1);
    f(2) = q * (-CL_w(vec) * vec(0) - CD_b(vec) * vec(2));
    return f;
}

Vector3f Parafoil::W(const Vector3f &vec) {
    Vector3f W;
    W(0) = g * (-sin(vec(1)));
    W(1) = g * (sin(vec(0)) * cos(vec(1)));
    W(2) = g * (cos(vec(0)) * cos(vec(1)));
    return W;
}

Vector3f Parafoil::Ma_w(Vector3f &vel, Vector3f &w, Vector3f rot) {
    double q = 0.5*rho*sw*vel.squaredNorm();
    double vnorm = vel.norm();
    Vector3f m;
    m(0) = q*(Clp*pow(b,2)*w(0)/(2*vnorm) + Clphi*b*rot(0));
    m(1) = q*(Cmq*pow(c,2)*w(1)/(2*vnorm) + Cm0*c + Cma*c*atan2(vel(2), vel(0)));
    m(2) = q*(Cnr*pow(b,2)*w(2)/(2*vnorm));
    return m;
}

float Parafoil::CD_b(const Vector3f &vec)
{
    float alpha = atan2(vec(2), vec(0));
    float cd = CD0 + CDa * alpha*alpha;
    return cd;
}


float Parafoil::CL_w(const Vector3f &vec)
{
    float alpha = atan2(vec(2), vec(0));
    float cl = CL0 + CLa * alpha;
    return cl;
}


Matrix<float, 3, 2> Parafoil::SFa(const Vector3f &vec, float delta){

    float sign;
    if (delta >= 0.0)  sign = 1;
    else sign = -1;

    float q = 0.5*rho*sw*vec.norm();
    Matrix<float, 3, 2> F;
    F(0,0) = q*(CLda*vec(2)-CDda*vec(0))*sign;
    F(0,1) = q*(CLds*vec(2)-CDds*vec(0));
    F(1,0) = q*(-CDda*vec(1)*sign);
    F(1,1) = q*(-CDds*vec(1));
    F(2,0) = q*(-CLda*vec(0)-CDda*vec(2))*sign;
    F(2,1) = q*(-CLds*vec(0)-CDds*vec(2));
    return F;
}

Matrix<float, 3, 2> Parafoil::SMa(const Vector3f &vec){

    float q = 0.5*rho*sw*vec.squaredNorm()*b/th;
    Matrix<float, 3, 2> M;
    M(0,0) = q*Clda;
    M(0,1) = 0.0;
    M(1,0) = 0.0;
    M(1,1) = 0.0;
    M(2,0) = q*Cnda;
    M(2,1) = 0.0;

    return M;
}


Matrix<float, 3, 3> Parafoil::Wrot(const Vector3f &vec){
    Matrix3f Mat;
    Mat.row(0) << 1.0, sin(vec(0))*tan(vec(1)), cos(vec(0))*tan(vec(1));
    Mat.row(1) << 0.0, cos(vec(0)), -sin(vec(0));
    Mat.row(2) << 0.0, sin(vec(0))/cos(vec(1)), cos(vec(0))/cos(vec(1));
    return Mat;
}


Matrix<float, 3, 3> Parafoil::QuatToAtt(const Matrix<float, 4, 1> &q){
    Matrix<float,4 ,1 > quat = q.normalized();
    Vector3f qv;
    qv = quat.block<3, 1>(0, 0);
    Matrix<float, 3, 3> rox;
    Matrix<float,3,3> eye = Matrix<float, 3, 3>::Identity();
    rox << 0, -qv(2), qv(1),
            qv(2), 0, -qv(0),
            -qv(1), qv(0), 0;

    return (pow(quat(3),2)-qv.squaredNorm())*eye + 2*(qv*qv.transpose())-2*quat(3)*rox;
}


Vector3f Parafoil::QuatToEuler(const Matrix<float, 4, 1> &q){

    Matrix<float,4 ,1 > quat = q.normalized();
    Vector3f euler;
    euler(0) = atan2((2*quat(0)*quat(3)+2*quat(1)*quat(2)),(-pow(quat(0),2)-pow(quat(1),2)+pow(quat(2),2)+pow(quat(3),2)));
    euler(1) = -asin(2*quat(0)*quat(2)-2*quat(1)*quat(3));
    euler(2) = atan2((2*quat(0)*quat(1)+2*quat(2)*quat(3)),(pow(quat(0),2)-pow(quat(1),2)-pow(quat(2),2)+pow(quat(3),2)));

    return euler;
}

Matrix<float, 4, 4> Parafoil::Omega(const Vector3f &vel){
    Matrix<float, 4,4> Om;
    Om <<          0,  vel(2), -vel(1), vel(0),
            -vel(2),           0,  vel(0), vel(1),
            vel(1), -vel(0),           0, vel(2),
            -vel(0), -vel(1), -vel(2),          0;
    return Om;
}

Matrix<float, 4, 1> Parafoil::EulToQuat(const Vector3f &v){
    Matrix<float, 4, 1> quat;
    Vector3f vec = v/2;
    float c0 = cos(vec(0)), s0 = sin(vec(0));
    float c1 = cos(vec(1)), s1 = sin(vec(1));
    float c2 = cos(vec(2)), s2 = sin(vec(2));

    quat(3) = c0*c1*c2-s0*s1*s2;
    quat(0) = s0*c1*c2+c0*s1*s2;
    quat(1) = c0*s1*c2+s0*c1*s2;
    quat(2) = c0*c1*s2-s0*s1*c2;
    return quat;
}



void Parafoil::my_system(const state_type &x, state_type &dxdt, const double t) {

    Vector3d pos_c {x(0), x(1), x(2)};
    Matrix<float, 4,1> quat {x(3), x(4), x(5), x(6)};
    Vector3f vel_c {x(7), x(8), x(9)};
    Vector3f vel_rot {x(10), x(11), x(12)};

    Vector3f euler = QuatToEuler(quat);
    Vector3f vel_rot_eul = Wrot(euler) * vel_rot;


    Matrix3f T = QuatToAtt(quat);
    Matrix3f Sk_Om = Omrot(vel_rot);
    Vector3f Vw = vel_c - T * Wind();
    Vector3f vel_b = Vw + Sk_Om * xgb;
    Vector3f vel_p = Vw + Sk_Om * xgp;
    Vector2f sig {0.0, 0.0};

    if (t >= 10 && t < 1000) { sig << 20.0 * pi / 180, 0.0 * pi / 180; }
    if (sig(0) > sig(1)) { sig(0) = sig(0) - sig(1); }
    else {
        sig(1) = sig(0);
        sig(0) = sig(0) - sig(1);
    }

    Vector3f Fa = Fa_w(vel_p);
    Vector3f Fb = Fa_b(vel_b);
    Matrix<float, 3, 2> Sfa = SFa(vel_p, sig(0));
    dxdt.block<3, 1>(0, 0) = T.transpose() * vel_c;
    dxdt.block<4, 1>(3, 0) = 0.5 * Omega(vel_rot) * quat.normalized();
    dxdt.block<3, 1>(7, 0) = imto * (Fa + Fb + W(euler)*(mtot) -
                                      mtot * vel_rot.cross(vel_c) + Sfa * sig);
    dxdt.block<3, 1>(10, 0) = inv_inertia * (Ma_w(vel_p, vel_rot_eul, euler) + rgp * Fa + rgb * Fb -
                                     Sk_Om * (inv_inertia) * vel_rot + (SMa(vel_p) + rgp * Sfa) * sig);

}

void Parafoil::simulate(){
    double t0 = 0.0;
    double t = t0;

    vector<float> x_c, y_c, z_c, v_x, v_y, v_z, t_s, r_x, r_y, r_z;
    float z_land = 0.0;
    float z_ref = initial_state(2);

    state_type x=initial_state;
    Matrix<float, 4, 1> quat;
    Vector3f eul,vel;

    auto sys = std::bind(&Parafoil::my_system, std::ref(*this) , std::placeholders::_1 , std::placeholders::_2 , std::placeholders::_3 );

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

    //WindEstimation(t_s, v_x, v_y);

    fstream file;
    file.open("simulation_results.txt", ios_base::out );
    for (int i = 0; i<t_s.size(); i++)
    {
        file<< t_s[i] << " " << x_c[i] << " " << y_c[i] << " " << z_c[i] << " " << r_x[i] << " " << r_y[i] << " " << r_z[i] << " " << v_x[i] << " " << v_y[i] << " " << v_z[i] << endl;
    }
     file.close();

}