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

