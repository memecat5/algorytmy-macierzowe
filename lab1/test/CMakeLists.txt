add_executable(matrix_mult_tests
    test_binet.cpp
)

target_link_libraries(matrix_mult_tests
    PRIVATE
        matrix_mult
        GTest::gtest_main
)

include(GoogleTest)
gtest_discover_tests(matrix_mult_tests)
