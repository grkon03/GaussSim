#pragma once

#include "torus.hpp"
#include "natural.hpp"
#include <vector>
#include <functional>
#include <random>

namespace GaussSim
{
    using std::vector;

    // generalized gauss transformation
    template <Real R, size_t n>
    class GGT
    {
    protected:
        std::function<array<R, n>(array<R, n>)> originalTransformation;
        std::function<Torus<R, n>(Torus<R, n>)> transformation;

    public:
        GGT() {}
        GGT(std::function<array<R, n>(array<R, n>)> original);
        GGT(const GGT<R, n> &);

        // set N = numOfIteration, x = initial, then orbit() = {x, T(x), T^2(x), ..., T^{N-1}(x)} (N elements)
        vector<Torus<R, n>> orbit(Torus<R, n> initial, size_t numOfIteration) const;
        vector<array<NaturalNumber, n>> continuedFraction(Torus<R, n> target, size_t depth) const;

        double frequencyOfOrbit(Torus<R, n> rectBL, Torus<R, n> rectTR, Torus<R, n> initial, size_t depth) const;
        double frequencyOfRandomOrbits(Torus<R, n> rectBL, Torus<R, n> rectTR, size_t depth, size_t numOfExperiments) const;

        Torus<R, n> operator()(Torus<R, n> torus) const;
        array<R, n> operator()(array<R, n> coor) const;
    };

    template <Real R, size_t n>
    GGT<R, n>::GGT(std::function<array<R, n>(array<R, n>)> original)
        : originalTransformation(original)
    {
        transformation = [&](Torus<R, n> t)
        {
            return Torus<R, n>(originalTransformation(t.coordinate));
        };
    }

    template <Real R, size_t n>
    GGT<R, n>::GGT(const GGT<R, n> &ggt)
        : originalTransformation(ggt.originalTransformation), transformation(ggt.transformation) {}

    template <Real R, size_t n>
    vector<Torus<R, n>> GGT<R, n>::orbit(Torus<R, n> initial, size_t numOfIteration) const
    {
        vector<Torus<R, n>> orb;
        Torus<R, n> next(initial, false);

        int i;
        for (i = 0; i < numOfIteration; ++i)
        {
            orb.push_back(next);
            next = transformation(next);
        }

        return orb;
    }

    template <Real R, size_t n>
    vector<array<NaturalNumber, n>> GGT<R, n>::continuedFraction(Torus<R, n> arr, size_t depth) const
    {
        vector<array<NaturalNumber, n>> cf;

        auto orb = this->orbit(arr, depth);

        for (auto itr = orb.begin(); itr != orb.end(); ++itr)
        {
            cf.push_back(FloorArray<R, n>(originalTransformation(itr->coordinate)));
        }

        return cf;
    }

    template <Real R, size_t n>
    double GGT<R, n>::frequencyOfOrbit(Torus<R, n> rectBL, Torus<R, n> rectTR, Torus<R, n> initial, size_t depth) const
    {
        auto orb = this->orbit(initial, depth);
        int i;
        bool in;
        size_t timesOrbitComeToRect = 0;

        for (auto itr = orb.begin(); itr != orb.end(); ++itr)
        {
            in = true;
            for (i = 0; i < n; ++i)
            {
                if (!(rectBL[i] <= (*itr)[i] && (*itr)[i] <= rectTR[i]))
                {
                    // if not in the set
                    in = false;
                    break;
                }
            }
            if (in)
                ++timesOrbitComeToRect;
        }

        return (double)timesOrbitComeToRect / depth;
    }

    template <Real R, size_t n>
    double GGT<R, n>::frequencyOfRandomOrbits(Torus<R, n> rectBL, Torus<R, n> rectTR, size_t depth, size_t numOfExperiments) const
    {
        size_t i;
        int j;
        double sumOfFrequency = 0;
        std::mt19937 mt32;
        std::uniform_real_distribution<double> rndm(0, 1);
        Torus<R, n> x;
        for (i = 0; i < numOfExperiments; ++i)
        {
            for (j = 0; j < n; ++j)
            {
                x[j] = rndm(mt32);
            }

            sumOfFrequency += this->frequencyOfOrbit(rectBL, rectTR, x, depth);
        }

        return sumOfFrequency / numOfExperiments;
    }

    template <Real R, size_t n>
    Torus<R, n> GGT<R, n>::operator()(Torus<R, n> torus) const
    {
        return transformation(torus);
    }

    template <Real R, size_t n>
    array<R, n> GGT<R, n>::operator()(array<R, n> torus) const
    {
        return originalTransformation(torus);
    }

    const GGT<double, 1> normalGaussTransformation =
        GGT<double, 1>([](array<double, 1> x)
                       { return INVERSE<double, 1>(x); });
}