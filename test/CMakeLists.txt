cmake_minimum_required(VERSION 3.15)
project(zonoopt_unit_tests)

function(create_test TEST_NAME SOURCES)
    add_executable(${TEST_NAME} ${SOURCES})
    target_link_libraries(${TEST_NAME} PRIVATE ZonoOpt)
    add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
endfunction()

# tests
create_test(test_get_leaves "cpp-tests/test_get_leaves.cpp")
create_test(test_safety_verification "cpp-tests/test_safety_verification.cpp")