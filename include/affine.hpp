#pragma once
#ifndef GEOMLIB_AFFINE_HPP
#define GEOMLIB_AFFINE_HPP

#include <cstdint>
#include "../third_party/lalib/include/vec.hpp"
#include "../third_party/lalib/include/mat.hpp"
#include "curve/curve_concepts.hpp"

namespace geomlib {

template<size_t N> struct Affine;

template<Curve C, size_t N>
struct AffineTransformedCurve {
public:
    using PointType = lalib::VecD<N>;
    using VectorType = lalib::VecD<N>;

    AffineTransformedCurve(C&& curve, Affine<N>&& affine) noexcept;

    auto operator()(double s) const noexcept -> PointType;
    auto point(double s) const noexcept -> PointType;
    auto tangent(double s) const noexcept -> VectorType;

private:
    C _curve;
    Affine<N> _affine;
};


template<size_t N>
struct Affine {
public:
    Affine(const lalib::MatD<N, N>& c, const lalib::VecD<N>& b) noexcept;
    Affine(const lalib::MatD<N, N>& c, lalib::VecD<N>&& b) noexcept;
    Affine(lalib::MatD<N, N>&& c, const lalib::VecD<N>& b) noexcept;
    Affine(lalib::MatD<N, N>&& c, lalib::VecD<N>&& b) noexcept;

    auto transform(const lalib::VecD<N>& vec, lalib::VecD<N>& rslt) const noexcept -> lalib::VecD<N>&;
    auto transform(const lalib::DynVecD& vec, lalib::DynVecD& rslt) const noexcept -> lalib::DynVecD&;

private:
    lalib::MatD<N, N> _c;
    lalib::VecD<N> _b;
};

template<size_t N, size_t AXIS>
auto rotate(double angle) noexcept -> Affine<N> = delete;

template<size_t N>
auto translate(const lalib::VecD<N>& p) noexcept -> Affine<N>;

template<size_t N>
auto translate(lalib::VecD<N>&& p) noexcept -> Affine<N>;


template<Curve C, size_t N>
auto transform(C&& curve, Affine<N>&& affine) noexcept -> AffineTransformedCurve<C, N>;


// ###### Implementations ####### //

template <size_t N>
inline Affine<N>::Affine(const lalib::MatD<N, N> &c, const lalib::VecD<N> &b) noexcept:
    _c(c), _b(b)
{ }

template <size_t N>
inline Affine<N>::Affine(lalib::MatD<N, N> &&c, const lalib::VecD<N> &b) noexcept:
    _c(std::move(c)), _b(b)
{ }

template <size_t N>
inline Affine<N>::Affine(const lalib::MatD<N, N> &c, lalib::VecD<N> &&b) noexcept:
    _c(c), _b(std::move(b))
{ }

template <size_t N>
inline Affine<N>::Affine(lalib::MatD<N, N> &&c, lalib::VecD<N> &&b) noexcept:
    _c(std::move(c)), _b(std::move(b))
{ }

template <size_t N>
inline auto Affine<N>::transform(const lalib::VecD<N>& vec, lalib::VecD<N>& rslt) const noexcept -> lalib::VecD<N>& {
    lalib::mul(1.0, this->_c, vec, 0.0, rslt);
    lalib::add(rslt, this->_b, rslt);
    return rslt;
}

template <size_t N>
inline auto Affine<N>::transform(const lalib::DynVecD& vec, lalib::DynVecD& rslt) const noexcept -> lalib::DynVecD& {
    lalib::mul(1.0, this->_c, vec, 0.0, rslt);
    lalib::add(rslt, this->_b, rslt);
    return rslt;
}


template<>
inline auto rotate<2, 0>(double angle) noexcept -> Affine<2> {
    auto mat = lalib::MatD<2, 2>({
        std::cos(angle),    -std::sin(angle),
        std::sin(angle),    std::cos(angle) 
    });
    return Affine<2>(std::move(mat), lalib::VecD<2>::filled(0.0));
}

template<>
inline auto rotate<3, 0>(double angle) noexcept -> Affine<3> {
    auto mat = lalib::MatD<3, 3>({
        1.0,    0.0,                0.0,
        0.0,    std::cos(angle),    -std::sin(angle),
        0.0,    std::sin(angle),    std::cos(angle) 
    });
    return Affine<3>(std::move(mat), lalib::VecD<3>::filled(0.0));
}

template<>
inline auto rotate<3, 1>(double angle) noexcept -> Affine<3> {
    auto mat = lalib::MatD<3, 3>({
        std::cos(angle),    0.0,    std::sin(angle),
        0.0,                1.0,    0.0,
        -std::sin(angle),   0.0,    std::cos(angle) 
    });
    return Affine<3>(std::move(mat), lalib::VecD<3>::filled(0.0));
}

template<>
inline auto rotate<3, 2>(double angle) noexcept -> Affine<3> {
    auto mat = lalib::MatD<3, 3>({
        std::cos(angle),    -std::sin(angle),   0.0,
        std::sin(angle),    std::cos(angle),    0.0,
        0.0,                0.0,                1.0,
    });
    return Affine<3>(std::move(mat), lalib::VecD<3>::filled(0.0));
}

template<size_t N>
inline auto translate(const lalib::VecD<N>& p) noexcept -> Affine<2> {
    auto affine = Affine<N>(lalib::MatD<N, N>::diag(1.0), p);
    return affine;
}

template<size_t N>
inline auto translate(lalib::VecD<N>&& p) noexcept -> Affine<2> {
    auto affine = Affine<N>(lalib::MatD<N, N>::diag(1.0), std::move(p));
    return affine;
}

template <Curve C, size_t N>
auto transform(C &&curve, Affine<N> &&affine) noexcept -> AffineTransformedCurve<C, N>
{
    auto trans_curve = AffineTransformedCurve<C, N>(std::move(curve), std::move(affine));
    return trans_curve;
}

template <Curve C, size_t N>
inline AffineTransformedCurve<C, N>::AffineTransformedCurve(C &&curve, Affine<N> &&affine) noexcept:
    _curve(std::move(curve)), _affine(std::move(affine))
{}

template <Curve C, size_t N>
inline auto AffineTransformedCurve<C, N>::operator()(double s) const noexcept -> PointType
{
    return this->point(s);
}

template <Curve C, size_t N>
inline auto AffineTransformedCurve<C, N>::point(double s) const noexcept -> PointType
{
    auto p = this->_curve(s);
    this->_affine.transform(p, p);
    return p;
}

template <Curve C, size_t N>
inline auto AffineTransformedCurve<C, N>::tangent(double s) const noexcept -> VectorType
{
    auto tan = this->_curve.tangent(s);
    auto trans_tan = this->_affine.transform(tan);
    return trans_tan;
}
}

#endif