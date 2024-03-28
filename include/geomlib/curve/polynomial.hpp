#pragma once
#ifndef GEOMLIB_CURVE_POLYNOMIAL_HPP
#define GEOMLIB_CURVE_POLYNOMIAL_HPP

#include "lalib/vec.hpp"
#include <cstddef>

namespace geomlib {

template<size_t N>
struct PolynomialCurve {
    using PointType = lalib::VecD<N>;
    using VectorType = lalib::VecD<N>;

    PolynomialCurve(std::vector<lalib::VecD<N>>&& coeffs) noexcept;

    auto n_orders() const noexcept -> size_t;

    auto operator()(double s) const noexcept -> PointType;
    auto point(double s) const noexcept -> PointType;
    auto deriv(double s) const noexcept -> VectorType;
    auto tangent(double s) const noexcept -> VectorType;

private:
    std::vector<lalib::VecD<N>> _coeffs;
};


template<size_t N>
PolynomialCurve<N>::PolynomialCurve(std::vector<lalib::VecD<N>>&& coeffs) noexcept: 
    _coeffs(std::move(coeffs))
{ }


template<size_t N>
auto PolynomialCurve<N>::n_orders() const noexcept -> size_t {
    return this->_coeffs.size() - 1;
}

template<size_t N>
inline auto PolynomialCurve<N>::operator()(double s) const noexcept -> PointType {
    auto p = this->_coeffs.front();
    for (auto i = 1u; i < this->_coeffs.size(); ++i) {
        lalib::axpy(std::pow(s, i), this->_coeffs[i], p);
    }
    return p;
}

template<size_t N>
inline auto PolynomialCurve<N>::point(double s) const noexcept -> PointType {
    return (*this)(s);
}

template<size_t N>
inline auto PolynomialCurve<N>::deriv(double s) const noexcept -> VectorType {
    auto deriv = lalib::VecD<N>::filled(0.0);
    for (auto i = 0u; i < this->_coeffs.size(); ++i) {
        lalib::axpy(i * std::pow(s, i - 1), this->_coeffs[i], deriv);
    }
    return deriv;
}

template<size_t N>
inline auto PolynomialCurve<N>::tangent(double s) const noexcept -> VectorType {
    auto deriv = this->deriv(s);
    deriv = deriv / deriv.norm2();
    return deriv;
}

}

#endif