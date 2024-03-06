add_executable(geomlib_segment_test curve/segment.cc)
target_link_libraries(geomlib_segment_test PRIVATE ${BLAS_LIBRARIES} ${OpenMP_CXX_LIBRARIES} GTest::GTest GTest::Main)
gtest_discover_tests(geomlib_segment_test)

add_executable(geomlib_cubic_spline_test curve/cubic_spline.cc)
target_link_libraries(geomlib_cubic_spline_test PRIVATE ${BLAS_LIBRARIES} ${OpenMP_CXX_LIBRARIES} GTest::GTest GTest::Main)
gtest_discover_tests(geomlib_cubic_spline_test)

add_executable(geomlib_polyline_test curve/polyline.cc)
target_link_libraries(geomlib_polyline_test PRIVATE ${BLAS_LIBRARIES} ${OpenMP_CXX_LIBRARIES} GTest::GTest GTest::Main)
gtest_discover_tests(geomlib_polyline_test)

add_executable(geomlib_affine_test ./affine.cc)
target_link_libraries(geomlib_affine_test PRIVATE ${BLAS_LIBRARIES} ${OpenMP_CXX_LIBRARIES} GTest::GTest GTest::Main)
gtest_discover_tests(geomlib_affine_test)

add_executable(geomlib_curve_common_test ./curve_common.cc)
target_link_libraries(geomlib_curve_common_test PRIVATE ${BLAS_LIBRARIES} ${OpenMP_CXX_LIBRARIES} GTest::GTest GTest::Main)
gtest_discover_tests(geomlib_curve_common_test)