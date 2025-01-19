#pragma once

#include "real.hpp"
#include <algorithm>

namespace GaussSim
{
    // n x m matrix (n: num of rows, m: num of columns)
    template <RealSubgroup G, size_t n, size_t m>
    struct Matrix
    {
        array<array<G, m>, n> entries; // entries[i][j] = the entry at i-th row and j-th column

        // multiplication

        array<G, n> operator*(array<G, m>) const;

        template <size_t k>
        const Matrix<G, n, k> operator*(const Matrix<G, m, k> &) const;

        // operation

        const Matrix<G, m, n> transpose() const;
    };

    template <RealSubgroup G, size_t n, size_t m>
    array<G, n> Matrix<G, n, m>::operator*(array<G, m> vec) const
    {
        int i, j;
        array<G, n> ans;

        for (i = 0; i < n; ++i)
        {
            ans[i] = static_cast<G>(0);
            for (j = 0; j < m; ++j)
            {
                ans[i] += entries[i][j] * vec[j];
            }
        }

        return ans;
    }

    template <RealSubgroup G, size_t n, size_t m>
    template <size_t k>
    const Matrix<G, n, k> Matrix<G, n, m>::operator*(const Matrix<G, m, k> &mat) const
    {
        Matrix<G, n, k> ans;
        int a, b, c;
        for (a = 0; a < n; ++a)
        {
            for (b = 0; b < k; ++b)
            {
                ans[a][b] = static_cast<G>(0);
                for (c = 0; c < m; ++c)
                {
                    ans[a][b] += entries[a][c] * mat.entries[c][b];
                }
            }
        }

        return ans;
    }

    template <RealSubgroup G, size_t n, size_t m>
    const Matrix<G, m, n> Matrix<G, n, m>::transpose() const
    {
        Matrix<G, m, n> res;

        int i, j;
        for (i = 0; i < n; ++i)
        {
            for (j = 0; j < m; ++j)
            {
                res.entries[j][i] = entries[i][j];
            }
        }

        return res;
    }

    template <RealSubgroup G, size_t n>
    struct SquareMatrix : public Matrix<G, n, n>
    {
        operator Matrix<G, n, n>() const;
    };

    template <RealSubgroup G, size_t n>
    SquareMatrix<G, n>::operator Matrix<G, n, n>() const
    {
        Matrix<G, n, n> mat;
        mat.entries = this->entries;

        return mat;
    }
}