if(TARGET igl::igl)
  return()
endif()

message(STATUS "contour-tesselation: Adding target igl::igl")


include(FetchContent)

## ================ IGL

FetchContent_Declare(
  libigl
  GIT_REPOSITORY https://github.com/libigl/libigl.git
  GIT_TAG c7324d3fa6fa94e55764b274ab15949ebbd2471f
)
FetchContent_GetProperties(libigl)
if(NOT libigl_POPULATED)
  FetchContent_Populate(libigl)
endif()

add_library(igl INTERFACE)
target_include_directories(igl SYSTEM INTERFACE ${libigl_SOURCE_DIR}/include)

# using gcc, requires pthread
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
target_link_libraries(igl INTERFACE Threads::Threads)

# dependencies
include(deps/eigen)
target_link_libraries(igl INTERFACE Eigen3::Eigen)
include(deps/triangle)
target_link_libraries(igl INTERFACE triangle::triangle)
include(deps/predicates)
target_link_libraries(igl INTERFACE predicates::predicates)
include(deps/glfw)
target_link_libraries(igl INTERFACE glfw::glfw)
include(deps/glad)
target_link_libraries(igl INTERFACE glad::glad)
#include(deps/embree)
#target_link_libraries(igl INTERFACE embree::embree)

# kill warnings
# if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang" OR
#   "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang"
# )
#   target_compile_options(igl PRIVATE "-Wno-unused-function")
# endif()

## ================ igl::igl

add_library(igl::igl ALIAS igl)