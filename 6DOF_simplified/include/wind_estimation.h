#ifndef WIND_ESTIMATION_H
#define WIND_ESTIMATION_H

#include "constants.h"
#include <C:\Program Files\eigen-3.4.0\Eigen\Dense>

using namespace Eigen;
using namespace std;

Vector3f Wind();
void we(vector<double> time, vector<double> vx, vector<double> vy);
void we_calc(Matrix<double, 2000, 2> V);
void WindEstimation(vector<double> time, vector<double> vx, vector<double> vy);

#endif
