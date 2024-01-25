#pragma once
#ifndef GEOMLIB_AFFINE_CORE_HPP
#define GEOMLIB_AFFINE_CORE_HPP

#include <cstdint>
#include "../../third_party/lalib/include/vec.hpp"
#include "../../third_party/lalib/include/mat.hpp"

namespace geomlib {

template<size_t N>
struct Affine {
public:
    Affine(const lalib::MatD<N, N>& c, const lalib::VecD<N>& b) noexcept;
    Affine(const lalib::MatD<N, N>& c, lalib::VecD<N>&& b) noexcept;
    Affine(lalib::MatD<N, N>&& c, const lalib::VecD<N>& b) noexcept;
    Affine(lalib::MatD<N, N>&& c, lalib::VecD<N>&& b) noexcept;

    /// @brief Returns the matrix associated to the transformation except for the translation.
    auto mat() const noexcept -> const lalib::MatD<N, N>&;

    /// @brief Returns the vector associated to the translation
    auto vec() const noexcept -> const lalib::VecD<N>&;

    auto transform(const lalib::VecD<N>& vec, lalib::VecD<N>& rslt) const noexcept -> lalib::VecD<N>&;
    auto transform(lalib::VecD<N>& vec) const noexcept -> lalib::VecD<N>&;
    auto transformed(const lalib::VecD<N>& vec) const noexcept -> lalib::VecD<N>;

    auto transform(const lalib::DynVecD& vec, lalib::DynVecD& rslt) const noexcept -> lalib::DynVecD&;
    auto transform(lalib::DynVecD& vec) const noexcept -> lalib::DynVecD&;
    auto transformed(const lalib::DynVecD& vec) const noexcept -> lalib::DynVecD;

    auto composite(const Affine<N>& a) noexcept -> Affine<N>&;

private:
    lalib::MatD<N, N> _c;
    lalib::VecD<N> _b;
};


// ### Implementations ### //

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
inline auto Affine<N>::mat() const noexcept -> const lalib::MatD<N, N> &
{
    return this->_c;
}

template <size_t N>
inline auto Affine<N>::vec() const noexcept -> const lalib::VecD<N> &
{
    return this->_b;
}

template <size_t N>
inline auto Affine<N>::transform(const lalib::VecD<N> &vec, lalib::VecD<N> &rslt) const noexcept -> lalib::VecD<N> &
{
    lalib::mul(1.0, this->_c, vec, 0.0, rslt);
    lalib::add(rslt, this->_b, rslt);
    return rslt;
}

template <size_t N>
inline auto Affine<N>::transform(lalib::VecD<N> &vec) const noexcept -> lalib::VecD<N> &
{
    this->transform(vec, vec);
    return vec;
}

template <size_t N>
inline auto Affine<N>::transformed(const lalib::VecD<N> &vec) const noexcept -> lalib::VecD<N>
{
    auto p = lalib::VecD<N>::uninit();
    this->transform(vec, p);
    return p;
}

template <size_t N>
inline auto Affine<N>::transform(const lalib::DynVecD& vec, lalib::DynVecD& rslt) const noexcept -> lalib::DynVecD& {
    lalib::mul(1.0, this->_c, vec, 0.0, rslt);
    lalib::add(rslt, this->_b, rslt);
    return rslt;
}

template <size_t N>
inline auto Affine<N>::transform(lalib::DynVecD &vec) const noexcept -> lalib::DynVecD &
{
    this->transform(vec, vec);
    return vec;
}

template <size_t N>
inline auto Affine<N>::transformed(const lalib::DynVecD &vec) const noexcept -> lalib::DynVecD
{
    auto p = lalib::DynVecD::uninit(N);
    this->transform(vec, p);
    return p;
}

template <size_t N>
inline auto Affine<N>::composite(const Affine<N> &a) noexcept -> Affine<N>&
{
    lalib::mul(1.0, a._c, this->_c, 0.0, this->_c);
    lalib::mul(1.0, a._c, this->_b, 0.0, this->_b);
    lalib::add(this->_b, a._b, this->_b);
    return *this;
}

}

#endif