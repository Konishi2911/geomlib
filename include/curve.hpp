#pragma once
#ifndef GEOMLIB_CURVE_HPP
#define GEOMLIB_CURVE_HPP

#include "curve/curve_concepts.hpp"
#include "curve/cubic_spline.hpp"
#include "curve/polyline.hpp"
#include "curve/segment.hpp"
#include "../third_party/mathlib/include/roots/secant.hpp"

namespace geomlib {

/// @brief 
/// @tparam F 
/// @tparam N 
/// @param spline 
/// @param n 
/// @param f 
/// @return 
template<std::invocable<double> F, size_t N>
auto sample(const geomlib::CubicSpline<N>& spline, size_t n, F&& f) noexcept -> geomlib::Polyline<N> {
    auto root_solver = mathlib::Secant(1e-8);
    auto tot_len = spline.length();

    auto vec = std::vector<lalib::VecD<N>>();
    vec.reserve(n);
    for (auto k = 0u; k < n; ++k) {
        auto s = static_cast<double>(k) / (n - 1);
        auto len = tot_len * f(s);
        auto t = root_solver.find_root([&](double x){ return len - spline.length(x); }, 0.0, 1.0).value();

        vec.emplace_back(spline.point(t));
    }
    auto polyline = Polyline<N>(vec);
    return polyline;
}


}

#endif