#pragma once
#ifndef GEOMLIB_CURVE_HPP
#define GEOMLIB_CURVE_HPP

#include "curve/curve_concepts.hpp"
#include "curve/cubic_spline.hpp"
#include "curve/polyline.hpp"
#include "curve/segment.hpp"
#include "../third_party/mathlib/include/roots/secant.hpp"
#include "../third_party/mathlib/include/integral/simpson.hpp"
#include "../third_party/mathlib/include/nlp/nelder_mead.hpp"

namespace geomlib {

/// @brief  Returns the approximate length of the given curve.
/// @tparam C 
/// @param ss start point [0.0, 1.0]
/// @param es end point [0.0, 1.0], `ss < es`
/// @param curve 
/// @param n 
/// @return 
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
auto sample(const C& curve, size_t n, F&& f, double tol) noexcept -> geomlib::Polyline<Dim<C>()> {
    constexpr size_t N = Dim<C>();
    auto root_solver = mathlib::Secant(tol);
    auto tot_len = curve.length();

    auto vec = std::vector<lalib::VecD<N>>();
    vec.reserve(n);
    for (auto k = 0u; k < n; ++k) {
        auto s = static_cast<double>(k) / (n - 1);
        auto len = tot_len * f(s);
        auto t = root_solver.find_root([&](double x){ return len - curve.length(x); }, 0.0, 1.0);
        assert(t);

        vec.emplace_back(curve.point(t.value()));
    }
    auto polyline = Polyline<N>(std::move(vec));
    return polyline;
}

/// @brief Creates a polyline object by sampling the point on the given curve with specified spacing function.
/// @tparam F   the type of spacing function 
/// @tparam C   the curve type that conforms to `LengthMeasurableCurve`
/// @param curve    a curve sampling the point 
/// @param n        the number of sampling point including endpoints. 
/// @param f        the spacing function invocable with double-type parameter that range is [0, 1]. 
/// @return 
template<Curve C, std::invocable<double> F>
auto sample(const C& curve, size_t n, F&& f, double tol) noexcept -> geomlib::Polyline<Dim<C>()> {
    constexpr size_t N = Dim<C>();

    constexpr size_t n_len = 100;
    auto root_solver = mathlib::Secant(tol);
    auto tot_len = approx_length(0.0, 1.0, curve, n_len);

    auto vec = std::vector<lalib::VecD<N>>();
    vec.reserve(n);
    for (auto k = 0u; k < n; ++k) {
        auto s = static_cast<double>(k) / (n - 1);
        auto len = tot_len * f(s);
        auto t = root_solver.find_root([&](double x){ return len - approx_length(0.0, x, curve, n_len); }, 0.0, 1.0);
        assert(t);

        vec.emplace_back(curve.point(t.value()));
    }
    auto polyline = Polyline<N>(std::move(vec));
    return polyline;
}


/// @brief  Searches the nearest point from the given query, and returns the local position for it on the given curve.
/// @tparam C       the curve type
/// @param curve    the curve searching the nearest point
/// @param query    query
/// @param nsegs    the number of segments using at the first step of nearest point search.
/// @return 
template<Curve C>
auto nearest_on(const C& curve, const typename C::PointType& query, size_t nsegs) -> double {
    auto solver = mathlib::nlp::NelderMead(1e-6);

    double s = 0.0;
    double min_distance = std::numeric_limits<double>::infinity();
    for (auto i = 0u; i < nsegs; ++i) {
        auto sb = static_cast<double>(i) / nsegs;
        auto se = static_cast<double>(i + 1) / nsegs;

        auto cost = [&](auto s){ 
            double cost = 0.0;
            if (s < sb || s > se) {
                cost = std::abs(s - sb) / (se - sb);
            }
            s = std::clamp(s, sb, se);
            cost += (query - curve.point(s)).norm2();
            return cost;
        };

        auto sol = solver.solve(sb, se, std::move(cost), 100);
        assert(sol);
        if (sol.final_cost() < min_distance) {
            s = sol.sol();
            min_distance = sol.final_cost();
        }
    }
    return s;
}


}

#endif