#pragma once
#ifndef GEOMLIB_CURVE_POLYLINE_HPP
#define GEOMLIB_CURVE_POLYLINE_HPP

#include "../../third_party/lalib/include/vec.hpp"
#include "../affine/affine_core.hpp"
#include "segment.hpp"
#include <array>
#include <vector>
#include <cassert>

namespace geomlib {

template<size_t N>
struct Polyline {
public:
    using PointType = lalib::SizedVec<double, N>;
    using VectorType = lalib::SizedVec<double, N>;

    struct SegmentViewIterator {
        SegmentViewIterator() noexcept;
        SegmentViewIterator(const std::vector<PointType>& nodes) noexcept;
        SegmentViewIterator(const std::vector<PointType>& nodes, size_t cursor) noexcept;
        SegmentViewIterator(const SegmentViewIterator&) noexcept = default;
        SegmentViewIterator(SegmentViewIterator&&) noexcept = default;

        static auto sentinel(const std::vector<PointType>& nodes) noexcept -> SegmentViewIterator;

        // Conforming to the weakly incrementable
        using value_type = SegmentView<N>;
        using difference_type = ptrdiff_t;
        auto operator=(const SegmentViewIterator&) noexcept -> SegmentViewIterator&;
        auto operator=(SegmentViewIterator&&) noexcept -> SegmentViewIterator& = default;
        auto operator++() noexcept -> SegmentViewIterator&;
        auto operator++(int) noexcept -> SegmentViewIterator;

        // Conforming to the input_output_iterator
        auto operator*() const noexcept -> SegmentView<N>;

        auto operator==(const SegmentViewIterator& other) const noexcept -> bool;
        auto operator!=(const SegmentViewIterator& other) const noexcept -> bool;

    private:
        std::optional<std::reference_wrapper<const std::vector<PointType>>> _nodes;
        size_t _cursor;
    };

    Polyline(const std::vector<PointType>& nodes);
    Polyline(std::vector<PointType>&& nodes);

    /// @brief  Returns the number of segments in the polyline
    /// @return the number of segments
    auto n_segments() const noexcept -> size_t;

    auto begin() const noexcept -> SegmentViewIterator;
    auto end() const noexcept -> SegmentViewIterator;

    auto operator()(double s) const noexcept -> PointType;

    /// @brief  Returns the point on the segment specified by the given local position `s`
    /// @note   `s` is clamped with 0.0 or 1.0 if it is out of range.
    /// @param s    local position `s` which range is [0.0, 1.0].
    /// @return     a point on the segment
    auto point(double s) const noexcept -> PointType;

    /// @brief  Returns the length of the segment.
    /// @return length of the segment
    auto length() const noexcept -> double;

    auto dp(double s) const noexcept -> VectorType;

    /// @brief  Returns the tangent vector
    /// @param s 
    /// @return 
    auto tangent(double s) const noexcept -> VectorType;

    auto transform(const Affine<N>& affine) noexcept -> Polyline<N>&;

    auto transformed(const Affine<N>& affine) const noexcept -> Polyline<N>;

private:
    std::vector<double> _cumul_length;
    std::vector<PointType> _nodes;

    static auto __calc_cumul_length(const std::vector<PointType>& nodes) noexcept -> std::vector<double>;
    static auto __lerp(const PointType& p1, const PointType& p2, double t) noexcept -> PointType;
    auto __get_seg_id(double local_pos) const noexcept -> size_t;
    auto __calc_local_seg_pos(double local_pos, size_t seg_id) const noexcept -> double;
    auto __calc_seg_len(size_t seg_id) const noexcept -> double;
};

template <size_t N>
inline Polyline<N>::Polyline(const std::vector<PointType> &nodes):
    _cumul_length(__calc_cumul_length(nodes)), _nodes(nodes)
{ }

template <size_t N>
inline Polyline<N>::Polyline(std::vector<PointType> &&nodes):
    _cumul_length(__calc_cumul_length(nodes)), _nodes(std::move(nodes))
{ }

template <size_t N>
inline auto Polyline<N>::n_segments() const noexcept -> size_t
{
    auto n_seg = this->_nodes.size() - 1;
    return n_seg;
}

template <size_t N>
inline auto Polyline<N>::begin() const noexcept -> SegmentViewIterator
{
    return SegmentViewIterator(this->_nodes);
}

template <size_t N>
inline auto Polyline<N>::end() const noexcept -> SegmentViewIterator
{
    return SegmentViewIterator::sentinel(this->_nodes);
}

template <size_t N>
inline auto Polyline<N>::operator()(double s) const noexcept -> PointType
{
    return this->point(s);
}

template <size_t N>
inline auto Polyline<N>::point(double s) const noexcept -> PointType
{
    auto sid = this->__get_seg_id(s);
    auto seg_s = this->__calc_local_seg_pos(s, sid);
    auto p = this->__lerp(this->_nodes[sid], this->_nodes[sid + 1], seg_s);
    return p;
}

template <size_t N>
inline auto Polyline<N>::length() const noexcept -> double
{
    return this->_cumul_length.back();
}

template <size_t N>
inline auto Polyline<N>::dp(double s) const noexcept -> VectorType
{
    auto sid = this->__get_seg_id(s);
    auto seg_s = this->__calc_local_seg_pos(s, sid);
    auto dp = (this->_nodes[sid + 1] - this->_nodes[sid]) / this->__calc_seg_len(sid);
    return dp;
}

template <size_t N>
inline auto Polyline<N>::tangent(double s) const noexcept -> VectorType
{
    auto dp = this->dp(s);
    dp = dp / dp.norm2();
    return dp;
}

template <size_t N>
inline auto Polyline<N>::transform(const Affine<N> &affine) noexcept -> Polyline &
{
    auto n = this->_nodes.size();
    for (auto i = 0u; i < n; ++i) {
        affine.transform(this->_nodes[i]);
    }
    this->_cumul_length = __calc_cumul_length(this->_nodes);
    return *this;
}

template <size_t N>
inline auto Polyline<N>::transformed(const Affine<N> &affine) const noexcept -> Polyline
{
    auto polyline = Polyline(*this).transform(affine);
    return polyline;
}

template <size_t N>
inline auto Polyline<N>::__calc_cumul_length(const std::vector<PointType> &nodes) noexcept -> std::vector<double>
{
    double total_len = 0.0;
    auto cumul_len = std::vector<double>();
    cumul_len.emplace_back(0.0);
    for (auto i = 1u; i < nodes.size(); ++i) {
        auto len = (nodes[i] - nodes[i - 1]).norm2();
        total_len += len;
        cumul_len.emplace_back(total_len);
    }
    return cumul_len;
}

template <size_t N>
inline auto Polyline<N>::__lerp(const PointType &p1, const PointType &p2, double t) noexcept -> PointType
{
    return p1 + t * (p2 - p1);
}

template <size_t N>
inline auto Polyline<N>::__calc_local_seg_pos(double local_pos, size_t seg_id) const noexcept -> double
{
    assert(seg_id < this->n_segments());
    auto seg_len = this->__calc_seg_len(seg_id);
    auto s = (local_pos * this->length() - this->_cumul_length[seg_id]) / seg_len;
    return s;
}

template <size_t N>
inline auto Polyline<N>::__calc_seg_len(size_t seg_id) const noexcept -> double
{
    auto seg_len = this->_cumul_length[seg_id + 1] - this->_cumul_length[seg_id];
    return seg_len;
}

template <size_t N>
inline auto Polyline<N>::__get_seg_id(double local_pos) const noexcept -> size_t
{
    local_pos = std::clamp(local_pos, 0.0, 1.0);
    auto len = this->length();
    for (auto i = 1u; i < this->_cumul_length.size(); ++i) {
        if (local_pos * len < this->_cumul_length[i]) {
            return i - 1;
        }
    }
    return this->n_segments() - 1;  // the last segment
}

template <size_t N>
inline Polyline<N>::SegmentViewIterator::SegmentViewIterator() noexcept:
    _nodes(std::nullopt), _cursor(0)
{ }

template <size_t N>
inline Polyline<N>::SegmentViewIterator::SegmentViewIterator(const std::vector<PointType> &nodes) noexcept : 
    _nodes(nodes), _cursor(0)
{ }

template <size_t N>
inline Polyline<N>::SegmentViewIterator::SegmentViewIterator(const std::vector<PointType> &nodes, size_t cursor) noexcept : 
    _nodes(nodes), _cursor(cursor)
{ }

template <size_t N>
inline auto Polyline<N>::SegmentViewIterator::sentinel(const std::vector<PointType> &nodes) noexcept -> SegmentViewIterator
{
    return SegmentViewIterator(nodes, nodes.size() - 1);
}

template <size_t N>
inline auto Polyline<N>::SegmentViewIterator::operator=(const SegmentViewIterator &iter) noexcept -> SegmentViewIterator &
{
    this->_nodes = iter._nodes;
    this->_cursor = iter._cursor;
    return *this;
}

template <size_t N>
inline auto Polyline<N>::SegmentViewIterator::operator++() noexcept -> SegmentViewIterator &
{
    ++this->_cursor;
    return *this;
}

template <size_t N>
inline auto Polyline<N>::SegmentViewIterator::operator++(int) noexcept -> SegmentViewIterator
{
    auto tmp = *this;
    this->_cursor++;
    return tmp;
}

template <size_t N>
inline auto Polyline<N>::SegmentViewIterator::operator*() const noexcept -> SegmentView<N>
{
    auto view = SegmentView(this->_nodes->get()[this->_cursor], this->_nodes->get()[this->_cursor + 1]);
    return view;
}

template <size_t N>
inline auto Polyline<N>::SegmentViewIterator::operator==(const SegmentViewIterator &other) const noexcept -> bool
{
    auto b = this->_cursor == other._cursor;
    return b;
}

template <size_t N>
inline auto Polyline<N>::SegmentViewIterator::operator!=(const SegmentViewIterator &other) const noexcept -> bool
{
    return !((*this) == other);
}
}

#endif