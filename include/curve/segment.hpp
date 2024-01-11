#pragma once 
#ifndef GEOMLIB_SEGMENT_HPP
#define GEOMLIB_SEGMENT_HPP

#include <array>
#include "../../third_party/lalib/include/vec.hpp"

namespace geomlib {

template<size_t N>
struct Segment {
public:
    using PointType = lalib::SizedVec<double, N>;
    using VectorType = lalib::SizedVec<double, N>;

    Segment(const PointType& s, const PointType& e) noexcept;


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
};


// ==== Implementation ==== //

template <size_t N>
inline Segment<N>::Segment(const PointType &s, const PointType &e) noexcept:
    _p({s, e})
{ }

template <size_t N>
inline auto Segment<N>::operator()(double s) const noexcept -> PointType
{
    return point(s);
}

template <size_t N>
inline auto Segment<N>::point(double s) const noexcept -> PointType
{
    s = std::clamp(s, 0.0, 1.0);
    auto p = s * this->_p[1] + (1.0 - s) * this->_p[0];
    return p;
}
template <size_t N>
inline auto Segment<N>::length() const noexcept -> double
{
    auto length = (_p[1] - _p[0]).norm2();
    return length;
}
template <size_t N>
inline auto Segment<N>::tangent(double) const noexcept -> VectorType
{
    auto t = (_p[1] - _p[0]) / this->length();
    return t;
}

}

#endif