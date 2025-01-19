#pragma once

#include <array>
#include <vector>

namespace GaussSim
{
    using std::array;
    using std::vector;
    namespace util
    {
        // extract i-th entry
        template <typename T, size_t n>
        vector<T> extract(vector<array<T, n>> vec, size_t i)
        {
            vector<T> res;
            for (auto itr = vec.begin(); itr != vec.end(); ++itr)
            {
                res.push_back((*itr)[i]);
            }

            return res;
        }
    }
}