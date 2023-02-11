#ifndef TOOLS_H
#define TOOLS_H

#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

Matrix3d Trot(const Vector3d v);
Matrix3d Omrot(const Vector3d v);
Matrix3d Wrot(const Vector3d v);

#endif 