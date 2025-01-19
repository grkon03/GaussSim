#pragma once

#include "real.hpp"
#include "matrix.hpp"

#include <functional>

namespace GaussSim
{
    template <Real R, size_t n>
    struct Torus
    {
    public:
        array<R, n> coordinate;
        bool noAdjust = false;

    public:
        Torus() { std::fill(coordinate.begin(), coordinate.end(), static_cast<R>(0)); }
        Torus(array<R, n> coor) : coordinate(coor) { adjust(); }
        Torus(array<R, n> coor, bool noAdjust) : coordinate(coor), noAdjust(noAdjust) { adjust(); }
        Torus(const Torus<R, n> &torus) : coordinate(torus.coordinate), noAdjust(torus.noAdjust) {}
        Torus(const Torus<R, n> &torus, bool noAdjust) : coordinate(torus.coordinate), noAdjust(noAdjust) {}

        void adjust();
        array<R, n> inverse() const;

        // get a coordinate

        R &operator[](size_t i);

        // group addition

        const Torus<R, n> operator+(const Torus<R, n> &) const;
        const Torus<R, n> operator-(const Torus<R, n> &) const;

        const Torus<R, n> &operator+=(const Torus<R, n> &);
        const Torus<R, n> &operator-=(const Torus<R, n> &);

        // parallel shift

        const Torus<R, n> operator+(R) const;
        const Torus<R, n> operator-(R) const;
        const Torus<R, n> &operator+=(R);
        const Torus<R, n> &operator-=(R);

        // scalar multiplication

        const Torus<R, n> operator*(int) const;
        const Torus<R, n> &operator*=(int);

    public:
        // static functions

        // the lebesgue measure of the rectangle with diagonal points a, b
        static R measure(Torus<R, n> a, Torus<R, n> b);

        // the integral on the rectangle with diagonal points a, b
        // N is the number of partitions (per a direction). As N be greater, the integral is more accurate but the calc speed is slower.
        static R integral(std::function<R(Torus<R, n>)> integrant, Torus<R, n> a, Torus<R, n> b, size_t N);
    };

    template <Real R, size_t n>
    void Torus<R, n>::adjust()
    {
        if (noAdjust)
            return;

        int i;
        for (i = 0; i < n; ++i)
        {
            coordinate[i] = MOD1<R>(coordinate[i]);
        }
    }

    template <Real R, size_t n>
    array<R, n> Torus<R, n>::inverse() const
    {
        array<R, n> res;

        int i;
        for (i = 0; i < n; ++i)
        {
            if (coordinate[i] == 0)
                res[i] = 0;
            else
                res[i] = 1 / coordinate[i];
        }
        return res;
    }

    template <Real R, size_t n>
    R &Torus<R, n>::operator[](size_t i)
    {
        return coordinate[i];
    }

    template <Real R, size_t n>
    const Torus<R, n> &Torus<R, n>::operator+=(const Torus<R, n> &t)
    {
        int i;
        for (i = 0; i < n; ++i)
        {
            coordinate[i] += t.coordinate[i];
        }

        adjust();
    }

    template <Real R, size_t n>
    const Torus<R, n> &Torus<R, n>::operator-=(const Torus<R, n> &t)
    {
        int i;
        for (i = 0; i < n; ++i)
        {
            coordinate[i] -= t.coordinate[i];
        }

        adjust();
    }

    template <Real R, size_t n>
    const Torus<R, n> Torus<R, n>::operator+(const Torus<R, n> &t) const
    {
        Torus<R, n> res = *this;
        return (res += t);
    }

    template <Real R, size_t n>
    const Torus<R, n> Torus<R, n>::operator-(const Torus<R, n> &t) const
    {
        Torus<R, n> res = *this;
        return (res -= t);
    }

    template <Real R, size_t n>
    const Torus<R, n> &Torus<R, n>::operator+=(R d)
    {
        int i;
        for (i = 0; i < n; ++i)
        {
            coordinate[i] += d;
        }

        adjust();
    }

    template <Real R, size_t n>
    const Torus<R, n> &Torus<R, n>::operator-=(R d)
    {
        int i;
        for (i = 0; i < n; ++i)
        {
            coordinate[i] -= d;
        }

        adjust();
    }

    template <Real R, size_t n>
    const Torus<R, n> Torus<R, n>::operator+(R d) const
    {
        Torus<R, n> res = *this;
        return (res += d);
    }

    template <Real R, size_t n>
    const Torus<R, n> Torus<R, n>::operator-(R d) const
    {
        Torus<R, n> res = *this;
        return (res -= d);
    }

    template <Real R, size_t n>
    const Torus<R, n> &Torus<R, n>::operator*=(int m)
    {
        int i;
        for (i = 0; i < n; ++i)
        {
            coordinate[i] *= m;
        }

        adjust();
    }

    template <Real R, size_t n>
    const Torus<R, n> Torus<R, n>::operator*(int m) const
    {
        Torus<R, n> res = *this;
        return (res *= m);
    }

    // toral homomorphism

    template <Real R, size_t n, size_t m>
    const Torus<R, m> operator*(Matrix<int, m, n> mat, Torus<R, n> tor)
    {
        int i, j;
        Torus<R, m> res;
        for (i = 0; i < m; ++i)
        {
            for (j = 0; j < n; ++j)
            {
                res.coordinate[i] += mat.entries[i][j] * tor.coordinate[j];
            }
        }

        res.adjust();

        return res;
    }

    template <Real R, size_t n>
    R Torus<R, n>::measure(Torus<R, n> a, Torus<R, n> b)
    {
        R area = static_cast<R>(1);
        int i;

        a.adjust();
        b.adjust();

        for (i = 0; i < n; ++i)
        {
            if (a[i] <= b[i])
                area = area * (b[i] - a[i]);
            else
                area = area * (b[i] + 1 - a[i]);
        }

        return area;
    }

    template <Real R, size_t n>
    R Torus<R, n>::integral(std::function<R(Torus<R, n>)> integrant, Torus<R, n> a, Torus<R, n> b, size_t N)
    {
        R value = static_cast<R>(0);
        Torus<R, n> searchingRectBL, searchingRectTR; // bottom-left, top-right
        array<size_t, n> searchingIndex;
        int i;

        searchingIndex.fill(0);

        std::function<void(void)> __SpecifyPoints = [&]()
        {
            for (i = 0; i < n; ++i)
            {
                searchingRectBL[i] =
                    (static_cast<R>(b[i]) * searchingIndex[i] + a[i] * static_cast<R>((N - searchingIndex[i]))) / static_cast<R>(N);
                searchingRectTR[i] =
                    (static_cast<R>(b[i]) * (searchingIndex[i] + 1) + a[i] * static_cast<R>((N - searchingIndex[i] - 1))) / static_cast<R>(N);
            }
            searchingRectBL.adjust();
            searchingRectTR.adjust();
        };

        std::function<void(void)> __CalcAprxIntegralOnRect = [&]()
        {
            R midvalue, area;
            Torus<R, n> midpoint;
            int i;
            for (i = 0; i < n; ++i)
            {
                midpoint[i] = (searchingRectBL[i] + searchingRectTR[i]) / 2;
            }

            midvalue = integrant(midpoint);
            area = Torus<R, n>::measure(searchingRectBL, searchingRectTR);

            value = value + midvalue * area;
        };

        std::function<void(void)>
            __NextRect = [&]()
        {
            ++searchingIndex[0];
            for (i = 0; i < n - 1; ++i)
            {
                if (searchingIndex[i] == N)
                {
                    searchingIndex[i] = 0;
                    ++searchingIndex[i + 1];
                }
            }
        };

        while (searchingIndex[n - 1] < N)
        {
            // specify the diagonal points
            __SpecifyPoints();

            // calculate the approximated integral on the rectangle and add to value
            __CalcAprxIntegralOnRect();

            // move to the next rectangle
            __NextRect();
        }

        return value;
    }

    template <size_t n>
    using DTorus = Torus<double, n>;

    namespace util
    {
        using std::vector;

        template <Real R, size_t n>
        vector<R> extract(vector<Torus<R, n>> vec, size_t i)
        {
            vector<R> res;
            for (auto itr = vec.begin(); itr != vec.end(); ++itr)
            {
                res.push_back((*itr).coordinate[i]);
            }

            return res;
        }
    }
}
