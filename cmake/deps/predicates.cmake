if(TARGET predicates::predicates)
  return()
endif()

message(STATUS "contour-tesselation: Adding target predicates::predicates")


include(FetchContent)

FetchContent_Declare(
  predicates
	GIT_REPOSITORY https://github.com/libigl/libigl-predicates.git
	GIT_TAG 5a1d2194ec114bff51d5a33230586cafb83adc86
)
FetchContent_MakeAvailable(predicates)
add_library(predicates::predicates ALIAS predicates)
set_target_properties(predicates PROPERTIES FOLDER third_party)
