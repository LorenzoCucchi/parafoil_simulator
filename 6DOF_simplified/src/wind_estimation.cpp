#include <C:\Program Files\eigen-3.4.0\Eigen\Dense>
#include "constants.h"
#include "wind_estimation.h"

void we(vector<double> time, vector<double> vx, vector<double> vy)
{
    Matrix<double, 2000, 2> V;
    for (int i=0; i<time.size(); i++) {

        if (time[i] >= 20.0) {
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
    cout<< "Velocita' X:: "<<W(0)<<"  Velocita' Y:: "<<W(1)<< endl;

}