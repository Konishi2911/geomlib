#pragma once
#ifndef GEOMLIB_SURFACE_CONCEPTS_HPP
#define GEOMLIB_SURFACE_CONCEPTS_HPP

#include "lalib/vec.hpp"
#include <concepts>

namespace geomlib {

template<typename T>
concept Surface = requires(const T& t, const lalib::VecD<2>& s) {
    typename T::PointType;
    typename T::VectorType;
    { t.point(s) } -> std::convertible_to<typename T::PointType>;
    { t.normal(s) } -> std::convertible_to<typename T::VectorType>;
};

}

#endif