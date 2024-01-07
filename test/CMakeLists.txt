add_executable(segment_test curve/segment.cc)
target_link_libraries(segment_test PRIVATE ${BLAS_LIBRARIES} ${OpenMP_CXX_LIBRARIES} GTest::GTest GTest::Main)
gtest_discover_tests(segment_test)
