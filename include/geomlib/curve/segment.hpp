#pragma once 
#ifndef GEOMLIB_SEGMENT_HPP
#define GEOMLIB_SEGMENT_HPP

#include <array>
#include "lalib/vec.hpp"
#include "geomlib/affine/affine_core.hpp"

namespace geomlib {

template<size_t N> struct Segment;

template<size_t N>
struct SegmentView {
public:
    static constexpr size_t DIM = N;
    using PointType = lalib::SizedVec<double, N>;
    using VectorType = lalib::SizedVec<double, N>;

    SegmentView(const PointType& s, const PointType& e) noexcept;


    auto operator()(double s) const noexcept -> PointType;

    /// @brief  Returns the point on the segment specified by the given local position `s`
    /// @note   `s` is clamped with 0.0 or 1.0 if it is out of range.
    /// @param s    local position `s` which range is [0.0, 1.0].
    /// @return     a point on the segment
    auto point(double s) const noexcept -> PointType;

    /// @brief  Returns the derivative of the segment.
    auto deriv(double s) const noexcept -> VectorType;

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

    /// @brief calculates the minimum distance to the query point
    /// @param query    query point
    auto distance(const VectorType& query) const noexcept -> double;

    /// @brief calculates a foot of a perpendicular
    auto foot(const VectorType& query) const noexcept -> VectorType;

    /// @brief calculates the position in the local coordinate from the given point
    auto local(const VectorType& query) const noexcept -> double;

    /// @brief checks if this segment includes the query point
    /// @note   This function DO NOT check whether the query point is on the line belonging to this segment.
    auto check_inclusion(const VectorType& query) const noexcept -> bool;


private:
    std::reference_wrapper<const PointType> _ps;
    std::reference_wrapper<const PointType> _pe;
};

template<size_t N>
struct Segment {
public:
    static constexpr size_t DIM = N;
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

    auto deriv(double s) const noexcept -> VectorType;

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
    

    /// @brief calculates the minimum distance to the query point
    /// @param query    query point
    auto distance(const VectorType& query) const noexcept -> double;

    /// @brief calculates a foot of a perpendicular
    auto foot(const VectorType& query) const noexcept -> VectorType;

    /// @brief calculates the position in the local coordinate from the given point
    auto local(const VectorType& query) const noexcept -> double;

    /// @brief checks if this segment includes the query point
    /// @note   This function DO NOT check whether the query point is on the line belonging to this segment.
    auto check_inclusion(const VectorType& query) const noexcept -> bool;

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
inline auto SegmentView<N>::deriv(double) const noexcept -> VectorType
{
    return (this->_pe.get() - this->_ps.get());
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
inline auto SegmentView<N>::distance(const VectorType &query) const noexcept -> double
{
    auto v = this->tangent(0.0);
    auto t = (query - this->_ps.get()).dot(v);
    if (t > 1) { 
        auto d = (query - this->_pe.get()).norm2();
        return d;
    } 
    else if (t < 0) {
        auto d = (query - this->_ps.get()).norm2();
        return d;
    }
    else {
        auto d = ((query - this->_ps.get()) - t * v).norm2();
        return d;
    }
}

template<size_t N>
inline auto SegmentView<N>::foot(const VectorType& query) const noexcept -> VectorType {
    auto v = this->tangent(0.0);
    auto pf = this->_ps.get() + ((query - this->_ps.get()).dot(v)) * v;
    return pf;
}

template<size_t N>
inline auto SegmentView<N>::local(const VectorType& query) const noexcept -> double {
    auto s = (query - this->_ps.get()).dot(this->tangent(0.0)) / this->length();
    return s;
}

template<size_t N>
inline auto SegmentView<N>::check_inclusion(const VectorType& query) const noexcept -> bool {
    auto sign = (query - this->_ps.get()).dot(query - this->_pe.get());
    return sign < 0;
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
inline auto Segment<N>::deriv(double s) const noexcept -> VectorType
{
    return this->_segment_view.deriv(s);
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
template <size_t N>
inline auto Segment<N>::distance(const VectorType &query) const noexcept -> double
{
    return this->_segment_view.distance(query);
}

template<size_t N>
inline auto Segment<N>::foot(const VectorType& query) const noexcept -> VectorType {
    return this->_segment_view.foot(query);
}

template<size_t N>
inline auto Segment<N>::local(const VectorType& query) const noexcept -> double {
    return this->_segment_view.local(query);
}

template<size_t N>
inline auto Segment<N>::check_inclusion(const VectorType& query) const noexcept -> bool {
    return this->_segment_view.check_inclusion(query);
}
}

#endif