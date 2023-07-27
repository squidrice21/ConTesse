if(TARGET triangle::triangle)
  return()
endif()

message(STATUS "contour-tesselation: Adding target triangle::triangle")


include(FetchContent)

FetchContent_Declare(
  triangle
  GIT_REPOSITORY https://github.com/libigl/triangle.git
  GIT_TAG        d284c4a843efac043c310f5fa640b17cf7d96170
)
FetchContent_MakeAvailable(triangle)

FetchContent_GetProperties(triangle)
target_include_directories(triangle PUBLIC ${triangle_SOURCE_DIR})
add_library(triangle::triangle ALIAS triangle)
set_target_properties(triangle PROPERTIES FOLDER third_party)
