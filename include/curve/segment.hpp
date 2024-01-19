#pragma once 
#ifndef GEOMLIB_SEGMENT_HPP
#define GEOMLIB_SEGMENT_HPP

#include <array>
#include "../../third_party/lalib/include/vec.hpp"

namespace geomlib {

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

    /// @brief  Returns the tangent vector
    /// @param s 
    /// @return 
    auto tangent(double s) const noexcept -> VectorType;

private:
    const PointType& _ps;
    const PointType& _pe;
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

    /// @brief  Returns the tangent vector
    /// @param s 
    /// @return 
    auto tangent(double s) const noexcept -> VectorType;

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
    auto p = s * this->_pe + (1.0 - s) * this->_ps;
    return p;
}
template <size_t N>
inline auto SegmentView<N>::length() const noexcept -> double
{
    auto length = (_pe - _ps).norm2();
    return length;
}
template <size_t N>
inline auto SegmentView<N>::tangent(double) const noexcept -> VectorType
{
    auto t = (_pe - _ps) / this->length();
    return t;
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
inline auto Segment<N>::tangent(double s) const noexcept -> VectorType
{
    return this->_segment_view.tangent(s);
}

}

#endif