#include "tools.h"
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

Matrix3d Trot(Vector3d vec){
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

Matrix3d Omrot(Vector3d vec){
    Matrix3d Mat;
    Mat << 0.0, -vec(2), vec(1), vec(2), 0.0, -vec(0), -vec(1), vec(0), 0.0;
    return Mat;
}

Matrix3d Wrot(Vector3d vec){
    Matrix3d Mat;
    Mat.row(0) << 1.0, sin(vec(0))*tan(vec(1)), cos(vec(0))*tan(vec(1));
    Mat.row(1) << 0.0, cos(vec(0)), -sin(vec(0));
    Mat.row(2) << 0.0, sin(vec(0))/cos(vec(1)), cos(vec(0))/cos(vec(1));
    return Mat;
}


Matrix3d QuatToAtt(Matrix<double, 4, 1> quat){
    quat = quat.normalized();
    Vector3d qv;
    qv = quat.block<3, 1>(0, 0);
    Matrix<double, 3, 3> rox;
    MatrixXd eye = Matrix<double, 3, 3>::Identity();
    rox << 0, -qv(2), qv(1), 
        qv(2), 0, -qv(0), 
        -qv(1), qv(0), 0;

    return (pow(quat(3),2)-qv.squaredNorm())*eye + 2*(qv*qv.transpose())-2*quat(3)*rox;
}


Vector3d QuatToEuler(Matrix<double, 4, 1> quat){
    quat = quat.normalized();
    Vector3d euler;
    euler(0) = atan2((2*quat(0)*quat(3)+2*quat(1)*quat(2)),(-pow(quat(0),2)-pow(quat(1),2)+pow(quat(2),2)+pow(quat(3),2)));
    euler(1) = -asin(2*quat(0)*quat(2)-2*quat(1)*quat(3));
    euler(2) = atan2((2*quat(0)*quat(1)+2*quat(2)*quat(3)),(pow(quat(0),2)-pow(quat(1),2)-pow(quat(2),2)+pow(quat(3),2)));

    return euler;
}

Matrix<double, 4, 4> Omega(Vector3d vel){
    Matrix<double, 4,4> Om;
    Om <<          0,  vel(2), -vel(1), vel(0),
            -vel(2),           0,  vel(0), vel(1),
             vel(1), -vel(0),           0, vel(2),
            -vel(0), -vel(1), -vel(2),          0;
    return Om;
}

Matrix<double, 4, 1> EulToQuat(Vector3d vec){
    Matrix<double, 4, 1> quat;
    vec = vec/2;
    quat(0) = cos(vec(0))*cos(vec(1))*cos(vec(2))+sin(vec(0))*sin(vec(1))*sin(vec(2));
    quat(1) = sin(vec(0))*cos(vec(1))*cos(vec(2))-cos(vec(0))*sin(vec(1))*sin(vec(2));
    quat(2) = cos(vec(0))*sin(vec(1))*cos(vec(2))+sin(vec(0))*cos(vec(1))*sin(vec(2));
    quat(3) = cos(vec(0))*cos(vec(1))*sin(vec(2))-sin(vec(0))*sin(vec(1))*cos(vec(2));

    return quat;
}