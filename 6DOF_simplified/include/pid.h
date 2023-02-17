
#ifndef INC_6DOF_SIMPLIFIED_PID_H
#define INC_6DOF_SIMPLIFIED_PID_H

#include <C:\Program Files\eigen-3.4.0\Eigen\Dense>

using namespace Eigen;

class Pid {
private:
    double err = 0.0;
    double umax = 0.05;
public:

    double GuidanceSys(Vector2d targ, Vector2d pos, double head);

    double Proportional(double er, double kp);

    double PI(double er, double kp, double ki);

    double PID(double er, double kp, double ki, double kd);

    double Saturation(double u);

    double dt = 0.1;
};


#endif //INC_6DOF_SIMPLIFIED_PID_H
