if(TARGET Eigen3::Eigen)
  return()
endif()

message(STATUS "contour-tesselation: Adding target Eigen3::Eigen")

include(FetchContent)
FetchContent_Declare(
  eigen
  GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
  GIT_TAG tags/3.3.7
  GIT_SHALLOW TRUE
)
FetchContent_GetProperties(eigen)
if(NOT eigen_POPULATED)
  FetchContent_Populate(eigen)
endif()
set(EIGEN_INCLUDE_DIRS ${eigen_SOURCE_DIR})


add_library(Eigen3::Eigen INTERFACE IMPORTED GLOBAL)

target_include_directories(Eigen3::Eigen SYSTEM INTERFACE ${EIGEN_INCLUDE_DIRS})
target_compile_definitions(Eigen3::Eigen INTERFACE
  -DEIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
)