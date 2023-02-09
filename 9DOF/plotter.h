#ifndef PLOTTER_H
#define PLOTTER_H

#include <sciplot/sciplot.hpp>

using namespace std;
using namespace sciplot;


int plot_flight(vector<double> x, vector<double> y, vector<double> z, vector<double> vx, vector<double> vy, vector<double> vz, vector<double> t);

#endif
