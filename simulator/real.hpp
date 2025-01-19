#pragma once

#include "natural.hpp"

#include <array>
#include <concepts>
#include <cmath>

namespace GaussSim
{

    int mod1(int a)
    {
        return 0;
    }
    double mod1(double a)
    {
        return std::fmod(a, 1.0);
    }

    NaturalNumber floor(int a)
    {
        return a;
    }
    NaturalNumber floor(double a)
    {
        return std::floor(a);
    }

    using std::array;

    template <typename T>
    concept Mod1DefinedO = requires(T a) { {mod1(a)}->std::same_as<T>; };
    template <typename T>
    concept Mod1DefinedI = requires(T a) { {a.mod1()}->std::same_as<T>; };
    template <typename T>
    concept Mod1Defined = Mod1DefinedO<T> || Mod1DefinedI<T>;

    template <typename T>
    concept FloorDefinedO = requires(T a) { {floor(a)} -> std::same_as<NaturalNumber>; };
    template <typename T>
    concept FloorDefinedI = requires(T a) { {a.floor()} -> std::same_as<NaturalNumber>; };
    template <typename T>
    concept FloorDefined = FloorDefinedO<T> || FloorDefinedI<T>;

    template <Mod1Defined T>
    T MOD1(T a)
    {
        if constexpr (Mod1DefinedO<T>)
        {
            return mod1(a);
        }
        else if constexpr (Mod1DefinedI<T>)
        {
            return a.mod1();
        }

        return a;
    }

    template <FloorDefined T>
    NaturalNumber FLOOR(T a)
    {
        if constexpr (FloorDefinedO<T>)
        {
            return floor(a);
        }
        else if constexpr (FloorDefinedI<T>)
        {
            return a.floor();
        }

        return 0;
    }

    template <FloorDefined T, size_t n>
    array<NaturalNumber, n> FloorArray(array<T, n> arr)
    {
        array<NaturalNumber, n> res;
        int i;
        for (i = 0; i < n; ++i)
        {
            res[i] = FLOOR<T>(arr[i]);
        }

        return res;
    }

    template <typename T>
    concept Arithmetic = requires(T a, T b) {
        { a + b } -> std::same_as<T>;
        { a - b } -> std::same_as<T>;
        { a *b } -> std::same_as<T>;
        { a / b } -> std::same_as<T>;
    };

    template <typename R>
    concept RealSubgroup =
        std::convertible_to<R, double> &&
        std::convertible_to<int, R> &&
        Mod1Defined<R> &&
        FloorDefined<R> &&
        Arithmetic<R>;

    template <typename R>
    concept Real =
        RealSubgroup<R> &&
        std::convertible_to<double, R>;

    template <RealSubgroup G, size_t n>
    array<G, n> INVERSE(array<G, n> arr)
    {
        array<G, n> res;
        int i;
        for (i = 0; i < n; ++i)
        {
            if (arr[i] == static_cast<G>(0))
                res[i] = static_cast<G>(0);
            else
                res[i] = static_cast<G>(1) / arr[i];
        }

        return res;
    }

    template <RealSubgroup G, size_t n>
    array<G, n> operator+(array<G, n> arr1, array<G, n> arr2)
    {
        int i;
        array<G, n> ans;
        for (i = 0; i < n; ++i)
        {
            ans[i] = arr1[i] + arr2[i];
        }

        return ans;
    }

    template <RealSubgroup G, size_t n>
    array<G, n> operator+(array<G, n> arr1, array<NaturalNumber, n> arr2)
    {
        int i;
        array<G, n> ans;
        for (i = 0; i < n; ++i)
        {
            ans[i] = arr1[i] + static_cast<G>(arr2[i]);
        }

        return ans;
    }

    template <RealSubgroup G, size_t n>
    array<G, n> operator*(array<G, n> arr, G scalar)
    {
        int i;
        array<G, n> ans;
        for (i = 0; i < n; ++i)
        {
            ans[i] = arr[i] * scalar;
        }

        return ans;
    }
}