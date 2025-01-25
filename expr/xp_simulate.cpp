#include "../simulator/gauss.hpp"
#include "../simulator/helper/filter.hpp"
#include <array>
#include <vector>
#include <random>
#include <iomanip>
#include <sstream>

#include <OpenADAPT/Plot/Canvas.h>

#include "option.hpp"

Options options;

using std::array;
using std::vector;

array<double, 1> xp(array<double, 1> x)
{
    if (x[0] == 0)
        return array<double, 1>{0};
    return array<double, 1>{std::pow(x[0], -options.p)};
}

std::string to_string_with_precision(double value)
{
    std::ostringstream out;
    out << std::fixed << std::setprecision(options.precision) << value;
    return out.str();
}

int main(int argc, char *argv[])
{
    // apply commandline arguments

    options = GetOptions(argc, argv);

    // setting

    GaussSim::GGT<double, 1> ggt(xp);
    GaussSim::Torus<double, 1> tl(true), br(true);

    int i, j;
    const int numOfIteration = options.iterationRate * options.N;

    double interval_width = (double)1 / options.N;

    vector<double> times(options.N, 0), densitySum(options.N, 0), densityMean(options.N, 0);

    std::random_device seed;
    std::mt19937_64 mt(seed());
    std::uniform_real_distribution<double> ud(0, 1);

    for (j = 0; j < options.numOfExperiments; ++j)
    {
    REDO:
        auto orbit = ggt.orbit(array<double, 1>{ud(mt)}, numOfIteration);

        if (orbit[orbit.size() - 1][0] == 0)
            goto REDO;

        // calculate density

#pragma omp parallel for private(tl, br)
        for (i = 0; i < options.N; ++i)
        {
            tl[0] = (double)i / options.N;
            br[0] = (double)(i + 1) / options.N;
            times[i] = (tl[0] + br[0]) / 2;
            densitySum[i] += ggt.frequencyOfOrbit(tl, br, orbit) / interval_width;
        }
    }

    for (i = 0; i < options.N; ++i)
    {
        densityMean[i] = densitySum[i] / options.numOfExperiments;
    }

    if (options.filtered)
        densityMean = GaussSim::helper::filter::collection::SmoothingFilter7.applyFilter(densityMean);

    // plot density

    adapt::Canvas2D canvas("P=" + to_string_with_precision(options.p) + ".png");

    canvas.SetTitle("p = " + to_string_with_precision(options.p));
    canvas.SetYRangeMin(0);
    canvas.SetXLabel("[0, 1]");
    canvas.SetYLabel("density");

    canvas.PlotPoints(times, densityMean, adapt::plot::s_lines, adapt::plot::notitle);
}