#pragma once 
#include <cstdint>
#include <string>
#include <utility>
#include <vector>
#include <span>
#include <concepts>

namespace geomlib::io {

enum class ByteOrder { BigEndian, LittleEndian };

template<std::integral I>
auto as_binary(const I& val) noexcept -> std::span<const uint8_t, sizeof(I)> {
	return std::span<const uint8_t, sizeof(I)>(
		reinterpret_cast<const uint8_t*>(&val), sizeof(I)
	);
}

template<std::floating_point F>
auto as_binary(const F& val) noexcept -> std::span<const uint8_t, sizeof(F)> {
	return std::span<const uint8_t, sizeof(F)>(
		reinterpret_cast<const uint8_t*>(&val), sizeof(F)
	);
}

template<typename CharT>
auto as_binary(const std::basic_string<CharT>& str) noexcept -> std::span<const uint8_t> {
	return std::span<const uint8_t>(
		reinterpret_cast<const uint8_t*>(str.c_str()), str.size()
	);
}

template<typename T>
concept available_as_binary = requires(const T& a) { as_binary(a); };

template<typename T, ByteOrder B = ByteOrder::LittleEndian>
requires available_as_binary<T>
auto into_binary(const T& val) noexcept -> std::vector<uint8_t> {
	auto tmp = as_binary(val);
	auto vec = std::vector(tmp.begin(), tmp.end());
	if constexpr (B == ByteOrder::BigEndian) {
		std::reverse(vec.begin(), vec.end());
	}
	return vec;
}

template<typename T, ByteOrder B>
concept available_into_binary = requires(const T& a) { into_binary<T, B>(a); };

template<std::input_iterator It, ByteOrder B = ByteOrder::LittleEndian> 
requires available_into_binary<typename std::iterator_traits<It>::value_type, B>
auto into_binary(It begin, It end) noexcept -> std::vector<uint8_t> {
	auto buf = std::vector<uint8_t>();
	buf.reserve(std::distance(begin, end) * sizeof(typename std::iterator_traits<It>::value_type));

	for (; begin < end; ++begin) {
		auto tmp = into_binary<typename std::iterator_traits<It>::value_type, B>(*begin);
		buf.insert(buf.cend(), tmp.begin(), tmp.end());
	}
	return buf;
}

template<std::input_iterator It, ByteOrder B = ByteOrder::LittleEndian> 
requires std::input_iterator<decltype(std::begin(std::declval<typename std::iterator_traits<It>::value_type>()))>
auto into_binary(It begin, It end) noexcept -> std::vector<uint8_t> {
	auto buf = std::vector<uint8_t>();
	buf.reserve(std::distance(begin, end) * sizeof(typename std::iterator_traits<It>::value_type));

	for (; begin < end; ++begin) {
		using SubIt = decltype(std::begin(std::declval<typename std::iterator_traits<It>::value_type>()));
		auto tmp = into_binary<SubIt, B>((*begin).cbegin(), (*begin).cend());
		buf.insert(buf.cend(), tmp.cbegin(), tmp.cend());
	}
	return buf;
}

}
