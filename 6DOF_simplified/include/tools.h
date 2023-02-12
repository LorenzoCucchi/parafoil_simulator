#ifndef TOOLS_H
#define TOOLS_H

#include <cmath>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

Matrix3d Trot(const Vector3d v);
Matrix3d Omrot(const Vector3d v);
Matrix3d Wrot(const Vector3d v);
Matrix3d QuatToAtt(const Matrix<double, 4, 1> v);
Vector3d QuatToEuler(const Matrix<double, 4, 1> v);
Matrix<double, 4, 4> Omega(const Vector3d v);
Matrix<double, 4, 1> EulToQuat(const Vector3d v);
Vector3d Wind();
#endif 