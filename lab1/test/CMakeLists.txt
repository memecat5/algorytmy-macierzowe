add_executable(binet_tests
    test_binet.cpp
)

target_link_libraries(binet_tests
    PRIVATE
        multiply_binet
        GTest::gtest_main
)

include(GoogleTest)
gtest_discover_tests(binet_tests)


# ===== STRASSEN TESTS =====
add_executable(strassen_tests test_strassen.cpp)

target_link_libraries(strassen_tests
    PRIVATE
        multiply_strassen
        GTest::gtest_main)

gtest_discover_tests(strassen_tests
    DISCOVERY_MODE PRE_TEST
    PROPERTIES LABELS "strassen"
)

add_test(NAME all_binet_tests COMMAND binet_tests)
add_test(NAME all_strassen_tests COMMAND strassen_tests)