#include "plotter.h"
#include <sciplot/sciplot.hpp>

using namespace std;
using namespace sciplot;

int plot_flight(vector<double> x, vector<double> y, vector<double> z, vector<double> vx, vector<double> vy, vector<double> vz, vector<double> t)
{
    Plot2D plot0;
    plot0.xlabel("x [m]");
    plot0.ylabel("z [m]");
    plot0.drawCurve(x, z).label("Trajectory");

    Plot2D plot1;
    plot1.xlabel("t [s]");
    plot1.ylabel("th [degree]");
    plot1.drawCurve(t, vx).label("vx");

    Plot2D plot2;
    plot2.xlabel("t [s]");
    plot2.ylabel("velocity [m/s]");
    plot2.drawCurve(t, vx).label("vx");
    plot2.drawCurve(t, vy).label("vy");
    plot2.drawCurve(t, vz).label("vz");

    Plot2D plot3;
    plot3.xlabel("t [s]");
    plot3.ylabel("q [degree]");
    plot3.drawCurve(t, vz).label("vz");

    // Use the previous plots as sub-figures in a larger 2x2 figure.
    Figure fig = {{plot0, plot1},
                  {plot2, plot3}};

    fig.title("Rocket Flight");
    fig.palette("set1");

    // Create canvas to hold figure
    Canvas canvas = {{fig}};
    // Set canvas output size
    canvas.size(1600, 900);

    // Show the plot in a pop-up window
    canvas.show();

    // Save the figure to a PDF file
    canvas.save("example-multiplot.pdf");

    return 0;

}