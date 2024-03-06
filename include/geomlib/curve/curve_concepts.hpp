#pragma once 
#ifndef GEOMLIB_CURVE_CONCEPTS_HPP
#define GEOMLIB_CURVE_CONCEPTS_HPP

#include <cstddef>
#include <type_traits>
#include <concepts>

namespace geomlib {

template<typename T>
concept Curve = requires(const T& t) {
    typename T::PointType;
    typename T::VectorType;

    // Query to get points on the curve at the given normalized position.
    { t(std::declval<double>()) } -> std::convertible_to<typename T::PointType>;
    { t.point(std::declval<double>()) } -> std::convertible_to<typename T::PointType>;

    // Query to get the derivative of the curve at the specified position.
    { t.deriv(std::declval<double>()) } -> std::convertible_to<typename T::VectorType>;

    // Query to get the normalized tangent vector at the specified position.
    { t.tangent(std::declval<double>()) } -> std::convertible_to<typename T::VectorType>;
};

template<typename T>
concept LengthMeasurableCurve = Curve<T> &&
requires(const T& t) {
    // Query to get the length of the curve.
    { t.length() } -> std::convertible_to<double>;
    { t.length(std::declval<double>()) } -> std::convertible_to<double>;
};


template<Curve C>
constexpr auto Dim() noexcept -> size_t {
    return C::DIM;
}


}

#endif