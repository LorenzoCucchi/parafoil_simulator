#include "pid.h"
#include <C:\Program Files\eigen-3.4.0\Eigen\Dense>
#include "constants.h"
using namespace Eigen;
using namespace constants;

double Pid::GuidanceSys(Vector2d targ, Vector2d pos, double heading)
{
    Vector2d dir_corr = targ - pos;
    double psi_ref = atan2(dir_corr(1),dir_corr(0));
    double diff = heading - psi_ref;
    if (diff>pi){
        diff += -2*pi;
    }
    else if (diff < (-pi))
    {
        diff += 2*pi;
    }
    return -diff;
}

double Pid::Proportional(double er, double kp)
{
    Pid::err = er;
    return kp*er;
}

double Pid::PI(double er, double kp, double ki)
{
    double u = kp*er + ki*(err + er*dt);
    Pid::err = er;
    return u;
}

double Pid::PID(double er, double kp, double ki, double kd)
{
    double u = kp*er + ki*(err + er*dt) + kd*((er - err)/dt);
    err = er;
    return u;

}

double Pid::Saturation(double u)
{
    if (u > umax) u = umax;
    else if (u < -umax) u = -umax;
    else u = u;

    return u;
}