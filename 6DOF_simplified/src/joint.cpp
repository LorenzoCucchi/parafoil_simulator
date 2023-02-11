#include "joint.h"
#include "constants.h"
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;
using namespace constants;

Vector3d Mc(Vector3d vecB, Vector3d vecP, Vector3d vecwb, Vector3d vecwp){
    Vector3d m;
    Vector2d p, b;
    p = Eul_mod(vecP, vecwp);
    b = Eul_mod(vecB, vecwb);
    m << 0, 0, Kc*(p(0)-b(0))+Cc*(p(1)-b(1));  
    return m;
}

Vector2d Eul_mod(Vector3d vec, Vector3d vec2){
    double yaw = atan((sin(vec(0))*sin(vec(1))*sin(vec(2)) - cos(vec(0))*cos(vec(2)) )/(cos(vec(1)*cos(vec(2)))));
    double t_theta = (cos(vec(0))*sin(vec(1))*cos(vec(2)) + sin(vec(0))*sin(vec(2))*cos(yaw))/(cos(vec(1))*cos(vec(2)));
    double theta = atan(t_theta);
    Vector3d vec3;
    vec3 << -cos(theta)*t_theta, sin(theta)*t_theta, 1;
    double yaw_dot = vec3.transpose()*vec2;
    Vector2d res;
    res << yaw, yaw_dot;

    return res;
}

