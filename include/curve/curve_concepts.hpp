#pragma once 
#ifndef GEOMLIB_CURVE_CONCEPTS_HPP
#define GEOMLIB_CURVE_CONCEPTS_HPP

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

    // Query to get the length of the curve.
    { t.length() } -> std::convertible_to<double>;

    // Query to get the normalized tangent vector at the specified position.
    { t.tangent(std::declval<double>()) } -> std::convertible_to<typename T::VectorType>;
};

}

#endif