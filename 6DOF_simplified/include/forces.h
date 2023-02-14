#ifndef FORCES_H
#define FORCES_H

#include <C:\Program Files\eigen-3.4.0\Eigen\Dense>
#include <iostream>

using namespace std;
using namespace Eigen;
typedef Matrix<double, 3,2> Matrix32d;

Vector3d Fa_b(Vector3d& vec);

Vector3d W_b(Vector3d& vec);

double CD_b(Vector3d& vec);

Vector3d Fa_w(Vector3d& vec);

double CL_w(Vector3d& vec);

Vector3d W_w(Vector3d& vec);

Vector3d Ma_w(Vector3d& vel,Vector3d& w, Vector3d rot);

Matrix3d MF();

Matrix3d IF();

Matrix32d SFa(const Vector3d& vec, double delta);

Matrix32d SMa(const Vector3d& vec);

Matrix3d Ib();

Matrix3d Ip();
#endif 