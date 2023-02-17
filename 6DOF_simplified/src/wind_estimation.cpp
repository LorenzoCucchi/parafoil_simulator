#include <C:\Program Files\eigen-3.4.0\Eigen\Dense>
#include <iostream>
#include "constants.h"
#include "wind_estimation.h"

using namespace std;
using namespace Eigen;

Vector3d Wind(){
    Vector3d w;
    //srand((unsigned) time(nullptr));
    w << 0, 0, 0;
    //w(0) = w(0) + (rand()%5)/10.0;
    //w(1) = w(1) + (rand()%5)/10.0;

    return w;
}

void we(vector<double> time, vector<double> vx, vector<double> vy)
{
    Matrix<double, 2000, 2> V;
    for (int i=0; i<time.size(); i++) {

        if (time[i] >= 20.0 && time.size()>(i+2000)) {
            V.setZero();
            for (int j = 0; j < 1999;) {
                V.block<1, 1>(j, 0) << vx[j + i];
                V.block<1, 1>(j, 1) << vy[j + i];
                j++;
            }
                i += 2000;
                we_calc(V);
            }

        }
}



void we_calc(Matrix<double, 2000, 2> V)
{
    double V_x_m = V.col(0).mean();
    double V_y_m = V.col(1).mean();
    Matrix<double, 2000, 1> Vx, Vy;
    Vx = Vx.setOnes()*V_x_m;
    Vy = Vy.setOnes()*V_y_m;
    Matrix<double, 2000, 1> V2;
    V2 = V.col(0).array().pow(2)+V.col(1).array().pow(2);
    double V2_mean = V2.mean();

    Matrix<double, 2000, 2> A;
    Matrix<double, 2000, 1> b;
    A.col(0) = V.col(0) - Vx;
    A.col(1) = V.col(1) - Vy;
    b = 0.5*(V2.array()-V2_mean);
    Vector2d W = (A.transpose()*A).inverse()*(A.transpose()*b);
    //cout<< "Velocita' X:: "<<W(0)<<"  Velocita' Y:: "<<W(1)<< endl;

}

void WindEstimation(vector<double> time, vector<double> gpsN, vector<double> gpsE)
{
    Eigen::Vector2d phi, wind;
    wind.setZero();
    Eigen::Matrix<double, 1, 2> phiT;
    Eigen::Matrix<double, 2, 2> funv;
    funv.col(0) << 1.0, 0.0;
    funv.col(1) << 0.0, 1.0;
    Eigen::Vector2d temp;
    int nSample=0;
    double vx, vy, v2, y;

    for (int i=0; i<time.size(); i++) {

        if (time[i] >= 20.0) {
            nSample++;
            vx = (vx * nSample + gpsN[i]) / (nSample + 1);
            vy = (vy * nSample + gpsE[i]) / (nSample + 1);
            v2 = (v2 * nSample + (gpsN[i] * gpsN[i] + gpsE[i] * gpsE[i])) / (nSample + 1);
            phi(0) = gpsN[i] - vx;
            phi(1) = gpsE[i] - vy;
            y = 0.5f * ((gpsN[i] * gpsN[i] + gpsE[i] * gpsE[i]) - v2);

            phiT = phi.transpose();
            funv = (funv - (funv * phi * phiT * funv) / (1 + (phiT * funv * phi)));
            temp = (0.5 * (funv + funv.transpose()) * phi) * (y - phiT * wind);
            wind = wind + temp;
        }
    }
    cout<< "Velocita' X:: "<<wind(0)<<"  Velocita' Y:: "<<wind(1)<< endl;

}