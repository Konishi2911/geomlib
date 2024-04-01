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


## Surfaces

add_executable(geomlib_polygon_test surface/polygon.cc)
target_link_libraries(geomlib_polygon_test PRIVATE ${BLAS_LIBRARIES} ${OpenMP_CXX_LIBRARIES} GTest::GTest GTest::Main)
gtest_discover_tests(geomlib_polygon_test)

add_executable(geomlib_polygon_mesh_test surface/polygon_mesh.cc)
target_link_libraries(geomlib_polygon_mesh_test PRIVATE ${BLAS_LIBRARIES} ${OpenMP_CXX_LIBRARIES} GTest::GTest GTest::Main)
gtest_discover_tests(geomlib_polygon_mesh_test)


## IOs

add_executable(geomlib_binary_test io/binary_utility.cc)
target_link_libraries(geomlib_binary_test PRIVATE ${BLAS_LIBRARIES} ${OpenMP_CXX_LIBRARIES} GTest::GTest GTest::Main)
gtest_discover_tests(geomlib_binary_test)

add_executable(geomlib_stl_test io/stl.cc)
target_link_libraries(geomlib_stl_test PRIVATE ${BLAS_LIBRARIES} ${OpenMP_CXX_LIBRARIES} GTest::GTest GTest::Main)
gtest_discover_tests(geomlib_stl_test)


## Nearest Neighbour Searching

add_executable(geomlib_kdtree_test nns/kd_tree.cc)
target_link_libraries(geomlib_kdtree_test PRIVATE ${BLAS_LIBRARIES} ${OpenMP_CXX_LIBRARIES} GTest::GTest GTest::Main)
gtest_discover_tests(geomlib_kdtree_test)

