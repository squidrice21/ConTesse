if(TARGET nanovg::nanovg)
    return()
endif()

include(FetchContent)
FetchContent_Declare(
  nanovg
  GIT_REPOSITORY https://github.com/memononen/nanovg.git
  GIT_TAG        2bead03bea43b2418060aaa154f972829995e663  
)
message(STATUS "nanovg: creating target 'nanovg::nanovg'")
FetchContent_GetProperties(nanovg)
if(NOT nanovg_POPULATED)
  FetchContent_Populate(nanovg)
endif()

add_library(nanovg)
target_include_directories(nanovg PUBLIC ${nanovg_SOURCE_DIR}/src)
target_sources(nanovg PRIVATE 
  ${nanovg_SOURCE_DIR}/src/nanovg.c
)

add_library(nanovg::nanovg ALIAS nanovg)
set_target_properties(nanovg PROPERTIES FOLDER third_party)

