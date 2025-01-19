#include "simulator/reconstruct.hpp"

#include <iostream>
#include <string>
#include <cmath>

using namespace GaussSim;

// ergodic measure of the transformation x -> {alpha / x}
double gaussmeasPDF(double x, double alpha)
{
    return 1 / (std::log(1 + 1 / alpha) * (x + alpha));
}

int main()
{
    ReconstructGGT<double, 2> GT2D(
        [&](array<double, 2> x) -> array<double, 2>
        {
            return INVERSE(array<double, 2>{x[0] * x[1], x[1]});
        },
        [&](array<double, 2> x) -> array<double, 2>
        {
            return array<double, 2>{
                (x[0] == 0) ? 0 : x[1] / x[0],
                (x[1] == 0) ? 0 : 1 / x[1]};
        });

    const size_t accuracy = 40;
    const size_t numOfPartition = 10;
    const size_t numOfIteration = 1000;
    const size_t numOfExperiments = 1;

    int i, j;

    double expectedArea, calculatedArea;
    double btm, lft, top, rit;
    Torus<double, 2> bl, tr;

    // experiment

    std::cout << "Experiment : verify whether any ergodic s" << std::endl
              << "Assume     : there exists some ergodic measures" << std::endl
              << "Method     : calculate measure of rectangles by Birkhoff ergodic theorem" << std::endl
              << std::endl;

    for (i = 0; i < numOfPartition; ++i)
    {
        for (j = 0; j < numOfPartition; ++j)
        {
            btm = double(j) / numOfPartition;
            lft = double(i) / numOfPartition;
            top = double(j + 1) / numOfPartition;
            rit = double(i + 1) / numOfPartition;
            bl = Torus<double, 2>(array<double, 2>{lft, btm}, true);
            tr = Torus<double, 2>(array<double, 2>{rit, top}, true);
            std::cout << "rectangle " << "[" << double(i) / numOfPartition << ", " << double(i + 1) / numOfPartition << "]"
                      << "x" << "[" << double(j) / numOfPartition << ", " << double(j + 1) / numOfPartition << "]: "
                      << GT2D.frequencyOfRandomOrbits(bl, tr, numOfIteration, 100) << std::endl;
        }
    }
}