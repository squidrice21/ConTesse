cmake_minimum_required(VERSION 3.15)
if(COMMAND cmake_policy)
    cmake_policy(SET CMP0054 NEW) # Only interpret if() arguments as variables or keywords when unquoted.
    cmake_policy(SET CMP0076 NEW) # target_sources() command converts relative paths to absolute.
    # Warnings about the version keyword
    # We need to use both commands to also silence it for subprojects
    cmake_policy(SET CMP0048 NEW)
    set(CMAKE_POLICY_DEFAULT_CMP0048 NEW)
    # Enable macos RPATH
    cmake_policy(SET CMP0042 NEW)
    set(CMAKE_POLICY_DEFAULT_CMP0042 NEW)
endif(COMMAND cmake_policy)

# Mac needs sdk 10.15 or newer.
# needs to go before project()
if(APPLE)
  set(CMAKE_OSX_DEPLOYMENT_TARGET "10.15" CACHE STRING "MacOSX sdk version")
endif()

project(ContourTessellation VERSION 0.1.1)

## OPTIONS -- use tbb static or dynamic
option(CONTESS_USE_TBB_STATIC    "Build with static TBB"   ON)

# Even though we are specifying 17 here, it seems that the linux builder 
# has an old g++ that does not have it.
set (CMAKE_CXX_STANDARD 17)
set (CMAKE_CXX_STANDARD_REQUIRED ON)
# Since we want to work on all dear OS's.
set (CMAKE_CXX_EXTENSIONS OFF)

# For FETCH CONTEXTS, put the dependencies in root folder to allow sharing them
# among debug/release builds
# If a value is specified when running cmake, it will override this.
set(FETCHCONTENT_BASE_DIR "${CMAKE_CURRENT_LIST_DIR}/_external" CACHE STRING "")

# Cmake will create a compile_commands.json file. It is useful for some editors
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

## Add the cmake/ folder to the module path so we can use the files in there.
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

## ============[ Compile with compiler warnings / for both us and deps ] =============

if (MSVC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP /D_USE_MATH_DEFINES")

  if (CMAKE_CXX_FLAGS MATCHES "/W[0-4]")
    string(REGEX REPLACE "/W[0-4]" "/W4" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
  else ()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
  endif ()
  # set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
  set(MSVC_RUNTIME "Dynamic")
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-unused-parameter")
  if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-gnu-anonymous-struct -Wno-c99-extensions -Wno-nested-anon-types")
  endif ()
endif ()

## ============[ Dependencies ] =============

include(deps/subdeval)
include(deps/tbb)
include(deps/cli11)
include(deps/spdlog)
include(deps/catch2)
include(deps/nlohmannjson)
include(deps/libigl)
include(deps/nanogui)
include(deps/opensubdiv)
include(deps/shor_van_wyck)

## ============[ Compile with compiler warnings / only for us ] =============

if(APPLE)
  set(CMAKE_MACOSX_RPATH ON)
  add_definitions(-DGL_SILENCE_DEPRECATION)
  add_definitions(-DCL_SILENCE_DEPRECATION)
elseif(WIN32)
  #add_definitions (/D "_CRT_SECURE_NO_WARNINGS")
  #what is this?
  #add_definitions (/D "__TBB_NO_IMPLICIT_LINKAGE")
endif()

## ================== [BaseContourTessellation] ============
# A base lib to be used by viewer and testss.
add_library(Base${PROJECT_NAME} STATIC)

file(GLOB CONTESS_SRCS 
  "src/PatchRotationAngle.cpp"
  "src/PatchRotationAngle.h"
  "src/aabb.h"
  "src/add_simplification_vertex.cpp"
  "src/add_simplification_vertex.h"
  "src/bvh.cpp"
  "src/bvh.h"
  "src/camera.cpp"
  "src/camera.h"
  "src/cdt_non_SO_component.cpp"
  "src/cdt_non_SO_component.h"
  "src/chain_contour.cpp"
  "src/chain_contour.h"
  "src/collapse_flipped_edges.cpp"
  "src/collapse_flipped_edges.h"
  "src/common.h"
  "src/compute_2d_rotation_indices.cpp"
  "src/compute_2d_rotation_indices.h"
  "src/cut_patch_to_disk.cpp"
  "src/cut_patch_to_disk.h"
  "src/decompose_wso_triangulation.cpp"
  "src/decompose_wso_triangulation.h"
  "src/default_shapes.h"
  "src/detect_contour_loops.cpp"
  "src/detect_contour_loops.h"
  "src/evaluate_radial_curvature.cpp"
  "src/evaluate_radial_curvature.h"
  "src/far_utils.cpp"
  "src/far_utils.h"
  "src/fix_flipped_faces.cpp"
  "src/fix_flipped_faces.h"
  "src/foldertools.cpp"
  "src/foldertools.h"
  "src/inflat_path.cpp"
  "src/inflat_path.h"
  "src/inflat_path_sparse.cpp"
  "src/inflat_path_sparse.h"
  "src/insert_interpolated_contours.cpp"
  "src/insert_interpolated_contours.h"
  "src/insert_planar_map_intersections.cpp"
  "src/insert_planar_map_intersections.h"
  "src/insert_subdivision_cusps.cpp"
  "src/insert_subdivision_cusps.h"
  "src/label_patch_boundary.cpp"
  "src/label_patch_boundary.h"
  "src/lift_min_surface.cpp"
  "src/lift_min_surface.h"
  "src/lift_patch_components.cpp"
  "src/lift_patch_components.h"
  "src/logger.cpp"
  "src/logger.h"
  "src/make_cut_feasible.cpp"
  "src/make_cut_feasible.h"
  "src/mesh.cpp"
  "src/mesh.h"
  "src/orient_contour.cpp"
  "src/orient_contour.h"
  "src/ray_cast_QI.cpp"
  "src/ray_cast_QI.h"
  "src/remove_cut_infeasible_disk.cpp"
  "src/remove_cut_infeasible_disk.h"
  "src/remove_wso_failure.cpp"
  "src/remove_wso_failure.h"
  "src/sew_components.cpp"
  "src/sew_components.h"
  "src/shape_utils.cpp"
  "src/shape_utils.h"
  "src/shrink_cut.cpp"
  "src/shrink_cut.h"
  "src/simplify_contours.cpp"
  "src/simplify_contours.h"
  "src/stitch_mesh.cpp"
  "src/stitch_mesh.h"
  "src/subdiv_common.cpp"
  "src/subdiv_common.h"
  "src/subdiv_lacewell.cpp"
  "src/subdiv_lacewell.h"
  "src/subdiv_osd.cpp"
  "src/subdiv_osd.h"
  "src/subdivide_contour_edges.cpp"
  "src/subdivide_contour_edges.h"
  "src/subdivide_contour_edges_even.cpp"
  "src/subdivide_contour_edges_even.h"
  "src/svg.cpp"
  "src/svg.h"
  "src/sweepLine.h"
  "src/tag_2d_intersections.cpp"
  "src/tag_2d_intersections.h"
  "src/tag_concave_edges.cpp"
  "src/tag_concave_edges.h"
  "src/tag_cusp_facing.cpp"
  "src/tag_cusp_facing.h"
  "src/tag_invalid_cut_edges.cpp"
  "src/tag_invalid_cut_edges.h"
  "src/tag_non_fish_tail.cpp"
  "src/tag_non_fish_tail.h"
  "src/tag_simplification_edges.cpp"
  "src/tag_simplification_edges.h"
  "src/tag_simplification_edges_cases.cpp"
  "src/tag_simplification_edges_cases.h"
  "src/tag_simplification_edges_small_loops.cpp"
  "src/tag_simplification_edges_small_loops.h"
  "src/tag_subdivision_contour_edges.cpp"
  "src/tag_subdivision_contour_edges.h"
  "src/testing.cpp"
  "src/testing.h"
  "src/testing_cases.cpp"
  "src/testing_cases.h"
  "src/triangulate_wso_patches.cpp"
  "src/triangulate_wso_patches.h"
  "src/viewgraph.h"
  "src/io/*.h"
  "src/io/*.cpp"
)
target_sources(Base${PROJECT_NAME} PRIVATE
    external/tri_tri_intersection.cpp
    external/2dintersection.cpp
    ${CONTESS_SRCS}
)

target_include_directories(Base${PROJECT_NAME} PUBLIC
    src/
    # Needed to include generated files
    ${CMAKE_CURRENT_BINARY_DIR}
    # Surface mesh should be its own target, but for now let's just dump it here
    # as an include directory
    ${CMAKE_CURRENT_SOURCE_DIR}/external
)

target_link_libraries(Base${PROJECT_NAME} PUBLIC
    TBB::tbb
    opensubdiv::opensubdiv
    nanogui::nanogui
    spdlog::spdlog
    shor_van_wyck::shor_van_wyck
    nlohmannjson::nlohmannjson
    subdeval::subdeval
)
target_compile_definitions(Base${PROJECT_NAME} PUBLIC 
  "CONTOURS_TESSELATION_DATA_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}/data\""
)

## ================== [TestsContourTessellation] ============
file(GLOB CONTESS_TEST_SRCS 
  "src/tests/dump_cases.cpp"
  "src/tests/main.cpp"
  "src/tests/test_pipeline.cpp"
)
add_executable(Tests${PROJECT_NAME}
  ${CONTESS_TEST_SRCS}
)
target_link_libraries(Tests${PROJECT_NAME} PUBLIC 
  Base${PROJECT_NAME}
  CLI11::CLI11 # command line parser
  Catch2::Catch2 # Unit test framework
)
