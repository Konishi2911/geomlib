#pragma once
#ifndef GEOMLIB_CURVE_HPP
#define GEOMLIB_CURVE_HPP

#include "curve/curve_concepts.hpp"
#include "curve/cubic_spline.hpp"
#include "curve/polyline.hpp"
#include "curve/segment.hpp"
#include "../third_party/mathlib/include/roots/secant.hpp"
#include "../third_party/mathlib/include/integral/simpson.hpp"

namespace geomlib {

template<Curve C>
auto approx_length(double ss, double es, const C& curve, size_t n) -> double {
    auto integrator = mathlib::integral::Simpson();
    auto len = integrator.integrate(ss, es, [&](double s){ return curve.deriv(s).norm2(); }, n);
    return len;
}

/// @brief Creates a polyline object by sampling the point on the given curve with specified spacing function.
/// @tparam F   the type of spacing function 
/// @tparam C   the curve type that conforms to `LengthMeasurableCurve`
/// @param curve    a curve sampling the point 
/// @param n        the number of sampling point including endpoints. 
/// @param f        the spacing function invocable with double-type parameter that range is [0, 1]. 
/// @return 
template<LengthMeasurableCurve C, std::invocable<double> F>
auto sample(const C& curve, size_t n, F&& f) noexcept -> geomlib::Polyline<Dim<C>()> {
    constexpr size_t N = Dim<C>();
    auto root_solver = mathlib::Secant(1e-8);
    auto tot_len = curve.length();

    auto vec = std::vector<lalib::VecD<N>>();
    vec.reserve(n);
    for (auto k = 0u; k < n; ++k) {
        auto s = static_cast<double>(k) / (n - 1);
        auto len = tot_len * f(s);
        auto t = root_solver.find_root([&](double x){ return len - curve.length(x); }, 0.0, 1.0).value();

        vec.emplace_back(curve.point(t));
    }
    auto polyline = Polyline<N>(std::move(vec));
    return polyline;
}


}

#endif