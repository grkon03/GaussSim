#include <vector>

using std::vector;

namespace GaussSim::helper::filter
{
    // 1-dimentional local inner product filter
    class LIPFilter1D
    {
        vector<double> kernel;
        size_t kernelcenter;
        // if the filter is the smoothing filter of [1, 1, *1, 1, 1], then kernelcenter = 2
        // (which is the point of *, the center of the kernel matrix)

    public:
        LIPFilter1D(vector<double> kernel, size_t kernelcenter)
            : kernel(kernel), kernelcenter(kernelcenter) {}

        vector<double> applyFilter(vector<double> original) const;
    };

    namespace collection
    {
        const LIPFilter1D SmoothingFilter3({1, 1, 1}, 1);
        const LIPFilter1D SmoothingFilter5({1, 1, 1, 1, 1}, 2);
        const LIPFilter1D SmoothingFilter7({1, 1, 1, 1, 1, 1, 1}, 3);
        const LIPFilter1D GaussianFilter3({0.60653, 1, 0.60653}, 1);
        const LIPFilter1D GaussianFilter5({0.135335, 0.60653, 1, 0.60653, 0.135335}, 2);
    }

    vector<double> LIPFilter1D::applyFilter(vector<double> original) const
    {
        size_t ormin = 0, ormax = original.size() - 1;
        size_t i;
        long long kmin = -kernelcenter, kmax = kernel.size() - kernelcenter - 1;
        long long j, index;
        double weightsum;
        double liprod;

        vector<double> filtered(original.size(), 0);

        for (i = 0; i < original.size(); ++i)
        {
            weightsum = 0;
            liprod = 0;
            for (j = kmin; j <= kmax; ++j)
            {
                index = i + j;
                if (index < ormin || ormax < index)
                    continue;
                weightsum += kernel[j];
                liprod += kernel[j] * original[index];
            }
            filtered[i] = liprod / weightsum;
        }

        return filtered;
    }

    // 2-dimentional local inner product filter
    template <size_t kn, size_t km>
    class LIPFilter2D
    {
        Matrix<double, kn, km> kernel;
        std::pair<size_t, size_t> kernelcenter;

    public:
        LIPFilter2D(Matrix<double, kn, km> kernel, std::pair<size_t, size_t> kernelcenter)
            : kernel(kernel), kernelcenter(kernelcenter) {}

        template <size_t dn, size_t dm>
        Matrix<double, dn, dm> applyFilter(Matrix<double, dn, dm> original) const;
    };

    // あとで作る
    // template <size_t kn, size_t km>
    // template <size_t dn, size_t dm>
    // Matrix<double, dn, dm> LIPFilter2D<kn, km>::applyFilter(Matrix<double, dn, dm> original) const
    // {
    // }
}