#include "../simulator/reconstruct.hpp"
#include <iostream>
#include <cmath>
#include <random>

#include <OpenADAPT/Plot/Canvas.h>

using namespace GaussSim;

double p = 0.52, q = 0.48, r = 0.45, s = 0.55;

double phi(array<double, 2> x)
{
    if (x[0] * x[1] == 0)
        return 0;
    return 1 / (std::pow(x[0], p) * std::pow(x[1], q));
}
double psi(array<double, 2> x)
{
    if (x[0] * x[1] == 0)
        return 0;
    return 1 / (std::pow(x[0], r) * std::pow(x[1], s));
}

// ./md_gauss p q r s filename
int main(int argc, char *argv[])
{
    std::string filename = "md_gauss_cm";
    if (argc >= 5)
    {
        p = std::stod(argv[1]);
        q = std::stod(argv[2]);
        r = std::stod(argv[3]);
        s = std::stod(argv[4]);
    }
    if (argc >= 6)
    {
        filename = argv[5];
    }

    GGT<double, 2> GT2D(
        [&](array<double, 2> x) -> array<double, 2>
        {
            return array<double, 2>{phi(x), psi(x)};
        });

    const size_t accuracy = 40;
    const size_t numOfPartition = 100;
    const size_t numOfIteration = 100000;
    const size_t numOfExperiments = 10;

    int i, j, k;

    double expectedArea, calculatedArea;
    double btm, lft, top, rit;
    Torus<double, 2> bl, tr;

    std::random_device seed;
    std::mt19937_64 mt(seed());
    std::uniform_real_distribution<double> ud(0, 1);

    adapt::Matrix<double> calcedPDF(numOfPartition, numOfPartition);
    for (i = 0; i < numOfPartition; ++i)
    {
        for (j = 0; j < numOfPartition; ++j)
        {
            calcedPDF[i][j] = 0;
        }
    }

    // experiment

    for (k = 0; k < numOfExperiments; ++k)
    {
    REDO:
        auto orbit = GT2D.orbit(array<double, 2>{ud(mt), ud(mt)}, numOfIteration);
        if (orbit[numOfIteration - 1][0] == 0 || orbit[numOfIteration - 1][1] == 0)
            goto REDO;
#pragma omp parallel for private(btm, lft, top, rit, bl, tr)
        for (i = 0; i < numOfPartition; ++i)
        {
#pragma omp parallel for private(btm, lft, top, rit, bl, tr)
            for (j = 0; j < numOfPartition; ++j)
            {
                btm = double(j) / numOfPartition;
                lft = double(i) / numOfPartition;
                top = double(j + 1) / numOfPartition;
                rit = double(i + 1) / numOfPartition;
                bl = Torus<double, 2>(array<double, 2>{lft, btm}, true);
                tr = Torus<double, 2>(array<double, 2>{rit, top}, true);

                calcedPDF[i][j] += GT2D.frequencyOfOrbit(bl, tr, orbit) * numOfPartition * numOfPartition;
                // std::cout << calcedPDF[i][j] << std::endl;
            }
        }
    }

    calcedPDF /= numOfExperiments;

    // plot

    adapt::CanvasCM canvas(filename + ".png");
    std::pair<double, double> xrange = {0, 1}, yrange = {0, 1};

    // canvas.SetPaletteDefined({{0, "yellow"}, {0.4, "red"}, {0.8, "black"}, {1.2, "blue"}, {1.6, "cyan"}});
    canvas.SetSizeRatio(-1);
    canvas.SetXLabel("x");
    canvas.SetYLabel("y");
    canvas.SetXRange(0, 1);
    canvas.SetYRange(0, 1);
    canvas.SetCBRange(0, 2);
    canvas.SetTitle("p = " + std::to_string(p) + " q = " + std::to_string(q) + " r = " + std::to_string(r) + " s = " + std::to_string(s));
    canvas.PlotColormap(calcedPDF, xrange, yrange, adapt::plot::notitle);
}