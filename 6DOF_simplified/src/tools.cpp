#include "tools.h"
#include <cmath>
#include<stdlib.h>
#include <C:\Program Files\eigen-3.4.0\Eigen\Dense>

using namespace std;
using namespace Eigen;

[[maybe_unused]] Matrix3d Trot(const Vector3d& vec){
    Matrix3d Mat;
    Mat.row(0) << cos(vec(1))*cos(vec(2)), cos(vec(1))*sin(vec(2)), -sin(vec(1));
    Mat(1,0) = sin(vec(0))*sin(vec(1))*cos(vec(2)) - cos(vec(0))*sin(vec(2));
    Mat(1,1) = sin(vec(0))*sin(vec(1))*sin(vec(2)) + cos(vec(0))*sin(vec(2));
    Mat(1,2) = sin(vec(0))*cos(vec(1));
    Mat(2,0) = cos(vec(0))*sin(vec(1))*cos(vec(2)) + sin(vec(0))*sin(vec(2));
    Mat(2,1) = cos(vec(0))*sin(vec(1))*sin(vec(2)) - sin(vec(0))*cos(vec(2));
    Mat(2,2) = cos(vec(0))*cos(vec(1));
    return Mat;
}

Matrix3d Omrot(const Vector3d& vec){
    Matrix3d Mat;
    Mat << 0.0, -vec(2), vec(1), vec(2), 0.0, -vec(0), -vec(1), vec(0), 0.0;
    return Mat;
}

Matrix3d Wrot(const Vector3d& vec){
    Matrix3d Mat;
    Mat.row(0) << 1.0, sin(vec(0))*tan(vec(1)), cos(vec(0))*tan(vec(1));
    Mat.row(1) << 0.0, cos(vec(0)), -sin(vec(0));
    Mat.row(2) << 0.0, sin(vec(0))/cos(vec(1)), cos(vec(0))/cos(vec(1));
    return Mat;
}


Matrix3d QuatToAtt(const Matrix<double, 4, 1>& q){
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


Vector3d QuatToEuler(const Matrix<double, 4, 1>& q){

    Matrix<double,4 ,1 > quat = q.normalized();
    Vector3d euler;
    euler(0) = atan2((2*quat(0)*quat(3)+2*quat(1)*quat(2)),(-pow(quat(0),2)-pow(quat(1),2)+pow(quat(2),2)+pow(quat(3),2)));
    euler(1) = -asin(2*quat(0)*quat(2)-2*quat(1)*quat(3));
    euler(2) = atan2((2*quat(0)*quat(1)+2*quat(2)*quat(3)),(pow(quat(0),2)-pow(quat(1),2)-pow(quat(2),2)+pow(quat(3),2)));

    return euler;
}

Matrix<double, 4, 4> Omega(const Vector3d& vel){
    Matrix<double, 4,4> Om;
    Om <<          0,  vel(2), -vel(1), vel(0),
            -vel(2),           0,  vel(0), vel(1),
             vel(1), -vel(0),           0, vel(2),
            -vel(0), -vel(1), -vel(2),          0;
    return Om;
}

Matrix<double, 4, 1> EulToQuat(const Vector3d& v){
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

Vector3d Wind(){
    Vector3d w;
    //srand((unsigned) time(nullptr));
    w << -3, -2, 0;
    //w(0) = w(0) + (rand()%5)/10.0;
    //w(1) = w(1) + (rand()%5)/10.0;

    return w;
}