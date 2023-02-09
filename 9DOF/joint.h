#ifndef JOINT_H
#define JOINT_H

#include <Eigen/Dense>
#include "constants.h"

using namespace std;
using namespace Eigen;

Vector3d Mc(Vector3d vecB, Vector3d vecP, Vector3d vecwb, Vector3d vecwp);

Vector2d Eul_mod(Vector3d vec, Vector3d vec2);

#endif 