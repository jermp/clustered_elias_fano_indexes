cmake_minimum_required(VERSION 2.8)
project(CPEF)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

if (UNIX)

  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-missing-braces")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ggdb")

  if (USE_SANITIZERS)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fno-omit-frame-pointer")
  endif ()

endif()

find_package(Boost 1.58.0 COMPONENTS iostreams unit_test_framework filesystem system thread REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS})

add_subdirectory(stxxl)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${STXXL_CXX_FLAGS}")
include_directories(${STXXL_INCLUDE_DIRS})

add_subdirectory(succinct EXCLUDE_FROM_ALL)
add_subdirectory(FastPFor EXCLUDE_FROM_ALL)

add_executable(create_freq_index create_freq_index.cpp)
target_link_libraries(create_freq_index
  ${Boost_LIBRARIES}
  FastPFor_lib
  )

add_executable(create_clustered_freq_index_sb create_clustered_freq_index_sb.cpp)
target_link_libraries(create_clustered_freq_index_sb
  ${Boost_LIBRARIES}
  FastPFor_lib
  )

add_executable(create_clustered_freq_index_fb create_clustered_freq_index_fb.cpp)
target_link_libraries(create_clustered_freq_index_fb
  ${Boost_LIBRARIES}
  FastPFor_lib
  )

add_executable(compute_clusters compute_clusters.cpp)
target_link_libraries(compute_clusters
  ${Boost_LIBRARIES}
  FastPFor_lib
  )

add_executable(queries queries.cpp)
target_link_libraries(queries
  ${Boost_LIBRARIES}
  FastPFor_lib
  )

add_executable(create_wand_data create_wand_data.cpp)
target_link_libraries(create_wand_data
  ${Boost_LIBRARIES}
  )
