add_executable(segment_test curve/segment.cc)
target_link_libraries(segment_test PRIVATE ${BLAS_LIBRARIES} ${OpenMP_CXX_LIBRARIES} GTest::GTest GTest::Main)
gtest_discover_tests(segment_test)

add_executable(cubic_spline_test curve/cubic_spline.cc)
target_link_libraries(cubic_spline_test PRIVATE ${BLAS_LIBRARIES} ${OpenMP_CXX_LIBRARIES} GTest::GTest GTest::Main)
gtest_discover_tests(cubic_spline_test)

add_executable(affine_test ./affine.cc)
target_link_libraries(affine_test PRIVATE ${BLAS_LIBRARIES} ${OpenMP_CXX_LIBRARIES} GTest::GTest GTest::Main)
gtest_discover_tests(affine_test)