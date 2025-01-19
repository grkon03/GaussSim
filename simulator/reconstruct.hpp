#pragma once

#include "gauss.hpp"

namespace GaussSim
{
    template <Real R, size_t n>
    class ReconstructGGT : public GGT<R, n>
    {
    protected:
        std::function<array<R, n>(array<R, n>)> inverseOfOriginalTransformation;

    public:
        ReconstructGGT() {}
        ReconstructGGT(const GGT<R, n> ggt, std::function<array<R, n>(array<R, n>)> inverse)
            : GGT<R, n>(ggt), inverseOfOriginalTransformation(inverse) {}
        ReconstructGGT(
            std::function<array<R, n>(array<R, n>)> originalTransformation,
            std::function<array<R, n>(array<R, n>)> inverse)
            : GGT<R, n>(originalTransformation), inverseOfOriginalTransformation(inverse) {}

        array<R, n> reconstruct(Torus<R, n> value, size_t depthOfExpansion) const;
        array<R, n> reconstruct(vector<array<NaturalNumber, n>> expansion) const;
    };

    template <Real R, size_t n>
    array<R, n> ReconstructGGT<R, n>::reconstruct(Torus<R, n> value, size_t depthOfExpansion) const
    {
        auto cf = this->continuedFraction(value, depthOfExpansion);
        array<R, n> reconstructedValue;

        return reconstruct(cf);
    }

    template <Real R, size_t n>
    array<R, n> ReconstructGGT<R, n>::reconstruct(vector<array<NaturalNumber, n>> expansion) const
    {
        array<R, n> reconstructedValue;
        reconstructedValue.fill(static_cast<R>(0));
        for (auto ritr = expansion.rbegin(); ritr != expansion.rend(); ++ritr)
        {
            reconstructedValue = inverseOfOriginalTransformation(reconstructedValue + *ritr);
        }

        return reconstructedValue;
    }

    const ReconstructGGT<double, 1> reconstructNormalGT(
        normalGaussTransformation,
        [](array<double, 1> arr) -> array<double, 1>
        {
            return INVERSE<double, 1>(arr);
        });
}