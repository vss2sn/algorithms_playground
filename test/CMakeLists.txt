project(test)

find_package(Threads REQUIRED)
find_package(GTest REQUIRED)  # For GTest
enable_testing()

function(add_test_file TEST_SOURCE)
  get_filename_component(test_file_executable ${TEST_SOURCE} NAME_WE)
  add_executable(${test_file_executable})
  set_target_properties(${test_file_executable}
    PROPERTIES
      CXX_STANDARD 17
      CXX_STANDARD_REQUIRED YES
      CXX_EXTENSIONS NO
      CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -Wall -Werror"
  )
  target_sources(${test_file_executable} PRIVATE ${TEST_SOURCE})
  target_include_directories(${test_file_executable} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include ${GTEST_INCLUDE_DIRS} )
  target_link_libraries(${test_file_executable} PRIVATE algorithms_library gtest gtest_main Threads::Threads )
  add_test(${test_file_executable} ${test_file_executable})
endfunction()

add_test_file(${CMAKE_CURRENT_SOURCE_DIR}/src/sorting_tests.cpp)
