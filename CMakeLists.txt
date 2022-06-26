project(murk2)

set(CMAKE_CXX_STANDARD 20)

file(GLOB_RECURSE lib_source src/*.cpp)
file(GLOB_RECURSE lib_headers include/*.hpp)

find_package(Boost REQUIRED)
find_package(Threads REQUIRED)

add_subdirectory(c3lt)

add_library(murk2 ${lib_source})
target_include_directories(murk2 PUBLIC include)
target_link_libraries(murk2 PUBLIC c3lt Boost::headers Threads::Threads gmp)
#target_precompile_headers(murk2 PUBLIC ${lib_headers})

add_library(murk2::murk2 ALIAS murk2)

file(GLOB_RECURSE foo_source foo/*.cpp)
foreach(foo ${foo_source})
  file(RELATIVE_PATH foo_rel ${CMAKE_SOURCE_DIR}/foo ${foo})
  get_filename_component(foo_fname ${foo_rel} NAME_WE)
  get_filename_component(foo_dir ${foo_rel} DIRECTORY)
  string(REPLACE "/" "_" foo_dir ${foo_dir})

  set(foo_name ${foo_dir}_${foo_fname})

  add_executable(${foo_name} ${foo})

  target_link_libraries(${foo_name} murk2)
endforeach()

find_package(GTest REQUIRED)
include(GoogleTest)
enable_testing()
file(GLOB_RECURSE test_source test/*.cpp)
add_executable(murk2test ${test_source})
target_link_libraries(murk2test PUBLIC GTest::gtest_main murk2)
gtest_discover_tests(murk2test)
