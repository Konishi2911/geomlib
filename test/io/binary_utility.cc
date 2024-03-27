#include "geomlib/io/binary_utility.hpp"
#include <gtest/gtest.h>

TEST(BinaryUtilityTests, AsBinaryTest) {
	uint16_t uv16 = 65407;
	float fv32 = 1.0;
	std::string str = "abc";

	auto uv16_bin = geomlib::io::as_binary(uv16);	
	auto fv32_bin = geomlib::io::as_binary(fv32);
	auto str_bin = geomlib::io::as_binary(str);

	ASSERT_EQ(2, uv16_bin.size());
	EXPECT_EQ(0x7f, uv16_bin[0]);
	EXPECT_EQ(0xff, uv16_bin[1]);

	ASSERT_EQ(4, fv32_bin.size());
	EXPECT_EQ(0x00, fv32_bin[0]);
	EXPECT_EQ(0x00, fv32_bin[1]);
	EXPECT_EQ(0x80, fv32_bin[2]);
	EXPECT_EQ(0x3f, fv32_bin[3]);

	ASSERT_EQ(3, str_bin.size());
	EXPECT_EQ(0x61, str_bin[0]);
	EXPECT_EQ(0x62, str_bin[1]);
	EXPECT_EQ(0x63, str_bin[2]);
}

TEST(BinaryUtilityTests, IntoBinaryTest) {
	uint16_t uv16 = 65407;
	std::array<uint16_t, 3> arr = { 1024, 1025, 1026 };
	std::vector<uint16_t> vec = { 1024, 1025, 1026 };
	std::array<std::array<uint16_t, 2>, 2> h_arr = { 
		std::array<uint16_t, 2> {1024, 1025},
		std::array<uint16_t, 2> {1026, 1027}
	};
	std::vector<std::vector<uint16_t>> h_vec = { 
		std::vector<uint16_t> {1024, 1025},
		std::vector<uint16_t> {1026}
	};

	auto uv16_bin = geomlib::io::into_binary(uv16);
	auto uv16_bin_big = geomlib::io::into_binary<uint16_t, geomlib::io::ByteOrder::BigEndian>(uv16);
	auto arr_bin = geomlib::io::into_binary(arr.cbegin(), arr.cend());
	auto vec_bin = geomlib::io::into_binary(vec.cbegin(), vec.cend());
	auto h_arr_bin = geomlib::io::into_binary(h_arr.cbegin(), h_arr.cend());
	auto h_vec_bin = geomlib::io::into_binary(h_vec.cbegin(), h_vec.cend());

	ASSERT_EQ(2, uv16_bin.size());
	EXPECT_EQ(0x7f, uv16_bin[0]);
	EXPECT_EQ(0xff, uv16_bin[1]);

	ASSERT_EQ(2, uv16_bin_big.size());
	EXPECT_EQ(0xff, uv16_bin_big[0]);
	EXPECT_EQ(0x7f, uv16_bin_big[1]);

	ASSERT_EQ(6, arr_bin.size());
	EXPECT_EQ(0x00, arr_bin[0]);
	EXPECT_EQ(0x04, arr_bin[1]);
	EXPECT_EQ(0x01, arr_bin[2]);
	EXPECT_EQ(0x04, arr_bin[3]);
	EXPECT_EQ(0x02, arr_bin[4]);
	EXPECT_EQ(0x04, arr_bin[5]);

	ASSERT_EQ(6, vec_bin.size());
	EXPECT_EQ(0x00, vec_bin[0]);
	EXPECT_EQ(0x04, vec_bin[1]);
	EXPECT_EQ(0x01, vec_bin[2]);
	EXPECT_EQ(0x04, vec_bin[3]);
	EXPECT_EQ(0x02, vec_bin[4]);
	EXPECT_EQ(0x04, vec_bin[5]);

	ASSERT_EQ(8, h_arr_bin.size());
	EXPECT_EQ(0x00, h_arr_bin[0]);
	EXPECT_EQ(0x04, h_arr_bin[1]);
	EXPECT_EQ(0x01, h_arr_bin[2]);
	EXPECT_EQ(0x04, h_arr_bin[3]);
	EXPECT_EQ(0x02, h_arr_bin[4]);
	EXPECT_EQ(0x04, h_arr_bin[5]);
	EXPECT_EQ(0x03, h_arr_bin[6]);
	EXPECT_EQ(0x04, h_arr_bin[7]);

	ASSERT_EQ(6, h_vec_bin.size());
	EXPECT_EQ(0x00, h_vec_bin[0]);
	EXPECT_EQ(0x04, h_vec_bin[1]);
	EXPECT_EQ(0x01, h_vec_bin[2]);
	EXPECT_EQ(0x04, h_vec_bin[3]);
	EXPECT_EQ(0x02, h_vec_bin[4]);
	EXPECT_EQ(0x04, h_vec_bin[5]);
}
