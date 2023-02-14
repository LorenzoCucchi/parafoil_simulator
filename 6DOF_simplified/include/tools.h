#ifndef TOOLS_H
#define TOOLS_H

#include <cmath>
#include<stdlib.h>
#include <C:\Program Files\eigen-3.4.0\Eigen\Dense>

using namespace std;
using namespace Eigen;

[[maybe_unused]] Matrix3d Trot(const Vector3d& v);
Matrix3d Omrot(const Vector3d& v);
Matrix3d Wrot(const Vector3d& v);
Matrix3d QuatToAtt(const Matrix<double, 4, 1>& v);
Vector3d QuatToEuler(const Matrix<double, 4, 1>& v);
Matrix<double, 4, 4> Omega(const Vector3d& v);
Matrix<double, 4, 1> EulToQuat(const Vector3d& v);
Vector3d Wind();
#endif 