#pragma once 
#ifndef GEOMLIB_SEGMENT_HPP
#define GEOMLIB_SEGMENT_HPP

#include <array>
#include "../../third_party/lalib/include/vec.hpp"
#include "../affine/affine_core.hpp"

namespace geomlib {

template<size_t N> struct Segment;

template<size_t N>
struct SegmentView {
public:
    using PointType = lalib::SizedVec<double, N>;
    using VectorType = lalib::SizedVec<double, N>;

    SegmentView(const PointType& s, const PointType& e) noexcept;


    auto operator()(double s) const noexcept -> PointType;

    /// @brief  Returns the point on the segment specified by the given local position `s`
    /// @note   `s` is clamped with 0.0 or 1.0 if it is out of range.
    /// @param s    local position `s` which range is [0.0, 1.0].
    /// @return     a point on the segment
    auto point(double s) const noexcept -> PointType;

    /// @brief  Returns the length of the segment.
    /// @return length of the segment
    auto length() const noexcept -> double;

    /// @brief  Returns the length from starting point to the intermediate position `s`.
    /// @param s 
    auto length(double s) const noexcept -> double;

    /// @brief  Returns the tangent vector
    /// @param s 
    /// @return 
    auto tangent(double s) const noexcept -> VectorType;

    /// @brief transforms itself with given affine transformation object.
    /// @param affine 
    auto transformed(const Affine<N>& affine) const noexcept -> Segment<N>;

private:
    std::reference_wrapper<const PointType> _ps;
    std::reference_wrapper<const PointType> _pe;
};

template<size_t N>
struct Segment {
public:
    using PointType = lalib::SizedVec<double, N>;
    using VectorType = lalib::SizedVec<double, N>;

    Segment(const PointType& s, const PointType& e) noexcept;
    Segment(PointType&& s, PointType&& e) noexcept;

    auto operator()(double s) const noexcept -> PointType;

    /// @brief  Returns the point on the segment specified by the given local position `s`
    /// @note   `s` is clamped with 0.0 or 1.0 if it is out of range.
    /// @param s    local position `s` which range is [0.0, 1.0].
    /// @return     a point on the segment
    auto point(double s) const noexcept -> PointType;

    /// @brief  Returns the length of the segment.
    /// @return length of the segment
    auto length() const noexcept -> double;

    /// @brief  Returns the length from starting point to the intermediate position `s`.
    /// @param s 
    auto length(double s) const noexcept -> double;

    /// @brief  Returns the tangent vector
    /// @param s 
    /// @return 
    auto tangent(double s) const noexcept -> VectorType;

    
    /// @brief transforms itself with given affine transformation object.
    /// @param affine 
    auto transform(const Affine<N>& affine) noexcept -> Segment<N>&;

    /// @brief transforms itself with given affine transformation object.
    /// @param affine 
    auto transformed(const Affine<N>& affine) const noexcept -> Segment<N>;

private:
    std::array<lalib::SizedVec<double, N>, 2> _p;
    SegmentView<N> _segment_view;
};


// ==== Implementation ==== //

template <size_t N>
inline SegmentView<N>::SegmentView(const PointType &s, const PointType &e) noexcept:
    _ps(s), _pe(e)
{ }

template <size_t N>
inline auto SegmentView<N>::operator()(double s) const noexcept -> PointType
{
    return point(s);
}

template <size_t N>
inline auto SegmentView<N>::point(double s) const noexcept -> PointType
{
    s = std::clamp(s, 0.0, 1.0);
    auto p = s * this->_pe.get() + (1.0 - s) * this->_ps.get();
    return p;
}
template <size_t N>
inline auto SegmentView<N>::length() const noexcept -> double
{
    auto length = (_pe.get() - _ps.get()).norm2();
    return length;
}
template <size_t N>
inline auto SegmentView<N>::length(double s) const noexcept -> double
{
    auto len = s * this->length();
    return len;
}
template <size_t N>
inline auto SegmentView<N>::tangent(double) const noexcept -> VectorType
{
    auto t = (_pe.get() - _ps.get()) / this->length();
    return t;
}

template <size_t N>
inline auto SegmentView<N>::transformed(const Affine<N> &affine) const noexcept -> Segment<N>
{
    auto s = affine.transformed(this->_ps.get());
    auto e = affine.transformed(this->_pe.get());
    auto seg = Segment(std::move(s), std::move(e));
    return seg;
}

template <size_t N>
inline Segment<N>::Segment(const PointType &s, const PointType &e) noexcept:
    _p({s, e}), _segment_view(_p[0], _p[1])
{ }

template <size_t N>
inline Segment<N>::Segment(PointType &&s, PointType &&e) noexcept:
    _p({std::move(s), std::move(e)}), _segment_view(_p[0], _p[1])
{ }

template <size_t N>
inline auto Segment<N>::operator()(double s) const noexcept -> PointType
{
    return point(s);
}

template <size_t N>
inline auto Segment<N>::point(double s) const noexcept -> PointType
{
    return this->_segment_view.point(s);
}
template <size_t N>
inline auto Segment<N>::length() const noexcept -> double
{
    return this->_segment_view.length();
}
template <size_t N>
inline auto Segment<N>::length(double s) const noexcept -> double
{
    return this->_segment_view.length(s);
}
template <size_t N>
inline auto Segment<N>::tangent(double s) const noexcept -> VectorType
{
    return this->_segment_view.tangent(s);
}

template <size_t N>
inline auto Segment<N>::transform(const Affine<N> &affine) noexcept -> Segment<N> &
{
    affine.transform(this->_p[0]);
    affine.transform(this->_p[1]);
    this->_segment_view = SegmentView<N>(this->_p[0], this->_p[1]);
    return *this;
}
template <size_t N>
inline auto Segment<N>::transformed(const Affine<N> &affine) const noexcept -> Segment<N>
{
    auto seg = Segment(*this).transform(affine);
    return seg;
}
}

#endif