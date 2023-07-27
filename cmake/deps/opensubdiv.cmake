if(TARGET opensubdiv::opensubdiv)
    return()
endif()

message(STATUS "contess: creating target 'opensubdiv::opensubdiv'")

include(FetchContent)
FetchContent_Declare(
    opensubdiv
    GIT_REPOSITORY https://github.com/PixarAnimationStudios/OpenSubdiv.git
    GIT_TAG tags/v3_4_0
    GIT_SHALLOW TRUE
)

FetchContent_GetProperties(opensubdiv)
if(NOT opensubdiv_POPULATED)
    FetchContent_Populate(opensubdiv)
endif()

add_library(opensubdiv)
add_library(opensubdiv::opensubdiv ALIAS opensubdiv)

set_target_properties(opensubdiv PROPERTIES FOLDER third_party)

target_include_directories(opensubdiv PUBLIC ${opensubdiv_SOURCE_DIR})

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    target_compile_definitions(opensubdiv PUBLIC _USE_MATH_DEFINES)
endif()

file(GLOB SRC_FILES
    "${opensubdiv_SOURCE_DIR}/opensubdiv/far/*.h"
    "${opensubdiv_SOURCE_DIR}/opensubdiv/far/*.cpp"
    "${opensubdiv_SOURCE_DIR}/opensubdiv/sdc/*.h"
    "${opensubdiv_SOURCE_DIR}/opensubdiv/sdc/*.cpp"
    "${opensubdiv_SOURCE_DIR}/opensubdiv/vtr/*.h"
    "${opensubdiv_SOURCE_DIR}/opensubdiv/vtr/*.cpp"
    "${opensubdiv_SOURCE_DIR}/regression/shapes/*.h"
)
target_sources(opensubdiv PRIVATE ${SRC_FILES})

option(OPENSUBDIV_GREGORY_EVAL_TRUE_DERIVATIVES "Enable true derivative evaluation for Gregory basis patches" ON)
if(OPENSUBDIV_GREGORY_EVAL_TRUE_DERIVATIVES) 
  target_compile_definitions(opensubdiv PUBLIC "OPENSUBDIV_GREGORY_EVAL_TRUE_DERIVATIVES")
endif()

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang" OR
   "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    target_compile_options(opensubdiv PRIVATE
        "-Wno-unused-function"
    )
endif()