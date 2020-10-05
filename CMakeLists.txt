project(algorithms_playground)

# set( CMAKE_C_COMPILER /usr/bin/clang)
# set( CMAKE_CXX_COMPILER /usr/bin/clang++)

cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(algorithms_library)
target_include_directories(algorithms_library
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_sources(algorithms_library
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src/graph.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/graphAL.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/graphAM.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/binary_search_tree.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/binary_tree.cpp
    # sort is header only
)

set_target_properties(algorithms_library
  PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED YES
    CXX_EXTENSIONS NO
    CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -Werror"
)

add_executable(play ${CMAKE_CURRENT_SOURCE_DIR}/main/main.cpp)

set_target_properties(play
  PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED YES
    CXX_EXTENSIONS NO
    CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -Werror"
)

target_link_libraries(play PRIVATE algorithms_library)

option(RUN_TESTS "RUN_TESTS" ON)
if(RUN_TESTS)
  add_subdirectory(test)
  enable_testing()
endif(RUN_TESTS)
