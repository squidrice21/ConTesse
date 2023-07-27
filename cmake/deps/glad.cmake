if(TARGET glad::glad)
    return()
endif()

include(FetchContent)
FetchContent_Declare(
  glad
  GIT_REPOSITORY https://github.com/libigl/libigl-glad.git
  GIT_TAG        09b4969c56779f7ddf8e6176ec1873184aec890f  
)
message(STATUS "glad: creating target 'glad::glad'")

FetchContent_MakeAvailable(glad)

add_library(glad::glad ALIAS glad)

set_target_properties(glad PROPERTIES FOLDER third_party)

