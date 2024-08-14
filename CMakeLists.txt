cmake_minimum_required(VERSION 3.12)

# set the project name
project(2xornado VERSION 1.0.2)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(CMAKE_BUILD_TYPE)
else()
    set(CMAKE_BUILD_TYPE "Release")
endif()

# options

option(VERBOSITY "compile with verbosity outputs" ON)
if(VERBOSITY OR CMAKE_BUILD_TYPE STREQUAL "Debug") 
    add_definitions(-DVERBOSITY)
endif()

option(USE_LHGR "Use LHGR graph impl" ON)
if(USE_LHGR)
    message("using LHGR data structure for graph implementation!")
    add_definitions(-DUSE_LHGR)
endif()

option(USE_TRIE "Use trie for graph labels" ON)
if(USE_TRIE)
    message("using trie data structure for graph labels!")
    add_definitions(-DUSE_TRIE)
endif()

option(FULL_REDUCTION "always fully reduce vertex labels" ON) #DO NOT DEACTIVATE!
if(FULL_REDUCTION)
    message("use full-reduction for updating vertices")
    add_definitions(-DFULL_REDUCTION)
else()
    message( FATAL_ERROR "THERE ARE KNOWN BUGS when FULL_REDUCTION is OFF! -- DO NOT DEACTIVATE!")
endif()


#static linking to libstdc++ and libgcc in case of Release build; otherwise give warnings!
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall -Wextra -Wswitch -Wsign-compare -static-libgcc -static-libstdc++ -ltbb" )# -ftree-vectorizer-verbose=5 -Weffc++")


# add subdirectories
add_subdirectory(src)

# add the executable
add_executable(2xornado src/main.cpp src/argparse/argparse.hpp)

#find openmp package
find_package(OpenMP REQUIRED)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")


# link libs to 2xornado
target_link_libraries(2xornado PUBLIC ${EXTRA_LIBS} PRIVATE graph)
target_include_directories(
    2xornado PUBLIC
    "${PROJECT_BINARY_DIR}"
    ${EXTRA_INCLUDES}
    )

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-march=native" COMPILER_SUPPORTS_MARCH_NATIVE)
if(COMPILER_SUPPORTS_MARCH_NATIVE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -mtune=native")
endif()

##uncomment test and bench targets as test instances are not part of the repository
# set up unit tests with Catch2 - if found
find_package(Catch2 3 QUIET)
if(Catch2_FOUND)
    set_property(GLOBAL PROPERTY CTEST_TARGETS_ADDED ON)
    include(CTest)
    include(Catch)
    # These tests can use the Catch2-provided main
    add_executable(testing tests/test_xlit.cpp tests/test_xsys.cpp tests/test_graph.cpp tests/test_impl_graph.cpp tests/test_solve.cpp tests/test_trie.cpp)
    #add_executable(testing tests/test_trie.cpp)
    target_link_libraries(testing PRIVATE graph Catch2::Catch2WithMain graph)
    catch_discover_tests(testing)
    enable_testing()
endif()

## set up micro-benchmarks with google-bench
find_package(benchmark QUIET)
if(benchmark_FOUND)
    add_executable(bench tests/benchmark.cpp)
    target_link_libraries(bench PRIVATE benchmark::benchmark graph)
endif()

get_filename_component(BENCH_FILES "tests/2xnfs" ABSOLUTE CACHE)
add_definitions(-DBENCH_FILES="${BENCH_FILES}")




#add compile definitions to forward compile information to program
add_compile_definitions(__CMAKE_CXX_COMPILER_ID="${CMAKE_CXX_COMPILER_ID}")
add_compile_definitions(__CMAKE_CXX_COMPILER_VERSION="${CMAKE_CXX_COMPILER_VERSION}")
add_compile_definitions(__CMAKE_CXX_STANDARD="${CMAKE_CXX_STANDARD}")
string(TOUPPER ${CMAKE_BUILD_TYPE} CMAKE_BUILD_TYPE_UPPER)
add_compile_definitions(__CMAKE_CXX_FLAGS="${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE_UPPER}}")

add_compile_definitions(__PROJECT_NAME="${PROJECT_NAME}")
add_compile_definitions(__PROJECT_VERSION="${PROJECT_VERSION}")

