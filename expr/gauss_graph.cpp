#include "../simulator/reconstruct.hpp"
#include <iostream>
#include <cmath>
#include <random>

#include <OpenADAPT/Plot/Canvas.h>

using std::array;
using std::vector;

int main()
{
    GaussSim::GGT<double, 1> ga([](array<double, 1> x) -> array<double, 1>
                                { return GaussSim::INVERSE(x); });

    size_t N = 3628800;
    size_t i;
    GaussSim::Torus<double, 1> init;

    vector<double> x(N, 0);
    vector<double> y(N, 0);

    for (i = 0; i < N; ++i)
    {
        init = GaussSim::Torus<double, 1>(array<double, 1>{(double)i / N});
        x[i] = init[0];
        y[i] = ga(init)[0];
    }

    adapt::Canvas2D canvas("gauss_graph.png");

    canvas.SetTitle("ガウス変換T_1のグラフ");
    canvas.SetXRange(0, 1);
    canvas.SetYRange(0, 1);
    canvas.SetXLabel("x");
    canvas.SetYLabel("T_1(x)");
    canvas.SetSizeRatio(1);

    canvas.PlotPoints(x, y, adapt::plot::notitle, adapt::plot::s_lines);
}