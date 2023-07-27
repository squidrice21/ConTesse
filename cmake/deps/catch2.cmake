if(TARGET Catch2::Catch2)
    return()
endif()

include(FetchContent)
FetchContent_Declare(
    catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG v2.13.6
    GIT_SHALLOW TRUE
)
message(STATUS "contours-tesselation: creating target 'Catch2::Catch2'")
FetchContent_MakeAvailable(catch2)