#pragma once
#ifndef GEOMLIB_AFFINE_HPP
#define GEOMLIB_AFFINE_HPP

#include <cstdint>
#include "../third_party/lalib/include/vec.hpp"
#include "../third_party/lalib/include/mat.hpp"
#include "affine/affine_core.hpp"
#include "curve.hpp"

namespace geomlib {

template<typename T, size_t N>
concept AffineTransformableCurve = 
    Curve<T> &&
    requires(T& t) {
        { t.transform(std::declval<const Affine<N>&>()) };
    } &&
    requires(const T& t) {
        { t.transformed(std::declval<const Affine<N>&>()) } -> std::convertible_to<T>;
    };


template<Curve C, size_t N>
struct AffineTransformedCurve {
public:
    static constexpr size_t DIM = N;
    using PointType = lalib::VecD<N>;
    using VectorType = lalib::VecD<N>;

    AffineTransformedCurve(C&& curve, Affine<N>&& affine) noexcept;
    AffineTransformedCurve(AffineTransformedCurve<C, N>&&) noexcept = default;

    auto operator()(double s) const noexcept -> PointType;
    auto point(double s) const noexcept -> PointType;
    auto tangent(double s) const noexcept -> VectorType;

    auto transform(const Affine<N>& affine) noexcept -> AffineTransformedCurve<C, N>&;
    auto transformed(const Affine<N>& affine) const noexcept -> AffineTransformedCurve<C, N>;

    auto native() const noexcept -> C requires AffineTransformableCurve<C, N>;

    auto operator=(AffineTransformedCurve<C, N>&& other) noexcept -> AffineTransformedCurve<C, N>&;

private:
    C _curve;
    Affine<N> _affine;
};



/* ############################################################### *\
 * Helper functions to create the affine transformation object     *
\* ############################################################### */

template<size_t N, size_t AXIS>
auto rotate(double angle) noexcept -> Affine<N> = delete;

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

/// @brief  Creates a 2D-affine transformation object representing the rigid rotation.
/// @param angle    rotating angle (CCW is positive)
/// @param pivot    pivot position of the roattion 
/// @return 
inline auto rotate2d(double angle, const lalib::VecD<2>& pivot) noexcept -> Affine<2> {
    auto mat = lalib::MatD<2, 2>({
        std::cos(angle),    -std::sin(angle),
        std::sin(angle),    std::cos(angle) 
    });
    auto p = pivot - mat * pivot;

    return Affine<2>(std::move(mat), std::move(p));
}

/// @brief  Creates an affine transformation object representing isotropic scaling.
/// @tparam N   dimension
/// @param mag          magnification factor (e.g. 2 doubles geometries) 
/// @param translate    translation after scaling operation (default is zero)  
template<size_t N>
inline auto scale(double mag, lalib::VecD<N>&& translate = lalib::VecD<N>::filled(0.0)) noexcept -> Affine<N> {
    auto mat = lalib::MatD<N, N>::diag(mag);
    auto affine = Affine<N>(std::move(mat), std::move(translate));
    return affine;
}

/// @brief Creates an affine transformation object representing a translation.
/// @tparam N   dimension
/// @param p    translation vector 
template<size_t N>
inline auto translate(const lalib::VecD<N>& p) noexcept -> Affine<N> {
    auto affine = Affine<N>(lalib::MatD<N, N>::diag(1.0), p);
    return affine;
}

/// @brief Creates an affine transformation object representing a translation.
/// @tparam N   dimension
/// @param p    translation vector 
template<size_t N>
inline auto translate(lalib::VecD<N>&& p) noexcept -> Affine<N> {
    auto affine = Affine<N>(lalib::MatD<N, N>::diag(1.0), std::move(p));
    return affine;
}

/// @brief  Applies the given curve the specified transformation with lazy evaluation.
/// @tparam C   an arbitrary curve type that conforms to the concept `Curve`
/// @tparam N   dimension 
/// @param curve    a curve to be applied the transformation 
/// @param affine   an affine transformation object 
/// @return 
template <Curve C, size_t N>
auto transform_lazy(C &&curve, Affine<N> &&affine) noexcept -> AffineTransformedCurve<C, N>
{
    auto trans_curve = AffineTransformedCurve<C, N>(std::move(curve), std::move(affine));
    return trans_curve;
}

/// @brief    An overload to specialize in case the curve type is `AffineTransformedCurve`.
template <Curve C, size_t N>
auto transform_lazy(AffineTransformedCurve<C, N> &&curve, Affine<N> &&affine) noexcept -> AffineTransformedCurve<C, N>
{
    curve.transform(affine);
    return std::move(curve);
}

/// @brief  Applies the given curve the specified transformation with immidiate evaluation. The given curve object is overwrited after this operation completed.
/// @tparam C   an arbitrary curve type that conforms to the concept `Curve`
/// @tparam N   dimension 
/// @param curve    a curve to be applied the transformation 
/// @param affine   an affine transformation object 
/// @return 
template <Curve C, size_t N>
requires AffineTransformableCurve<C, N>
inline auto transform(C& curve, const Affine<N> &affine) noexcept -> C& {
    curve.transform(affine);
    return curve;
}

/// @brief  Applies the given curve the specified transformation with immidiate evaluation. The resultant curve is newly created and returned.
/// @tparam C   an arbitrary curve type that conforms to the concept `Curve`
/// @tparam N   dimension 
/// @param curve    a curve to be applied the transformation 
/// @param affine   an affine transformation object 
/// @return 
template <Curve C, size_t N>
requires AffineTransformableCurve<C, N>
inline auto transformed(const C& curve, const Affine<N> &affine) noexcept -> C {
    auto c = curve.transformed(affine);
    return c;
}

/// @brief  An overload for the `SegmentView`.
template <size_t N>
inline auto transformed(const SegmentView<N>& seg_view, const Affine<N>& affine) noexcept -> Segment<N> {
    auto seg = seg_view.transformed(affine);
    return seg;
}


/* ######################### *\
 * AffineTransformedCurve    *
\* ######################### */

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

template <Curve C, size_t N>
inline auto AffineTransformedCurve<C, N>::transform(const Affine<N> &affine) noexcept -> AffineTransformedCurve<C, N> &
{
    this->_affine.composite(affine);
    return *this;
}

template <Curve C, size_t N>
inline auto AffineTransformedCurve<C, N>::transformed(const Affine<N> &affine) const noexcept -> AffineTransformedCurve<C, N>
{
    auto curve = AffineTransformedCurve(*this).transform(affine);
    return curve;
}

template <Curve C, size_t N>
inline auto AffineTransformedCurve<C, N>::operator=(AffineTransformedCurve<C, N> &&other) noexcept -> AffineTransformedCurve<C, N> &
{
    this->_curve = std::move(other._curve);
    this->_affine = std::move(other._affine);
    return *this;
}

template <Curve C, size_t N>
inline auto AffineTransformedCurve<C, N>::native() const noexcept -> C
requires AffineTransformableCurve<C, N>
{
    auto native = this->_curve.transformed(this->_affine);
}
}

#endif