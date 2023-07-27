cmake_minimum_required(VERSION 3.16)

if(COMMAND cmake_policy)
    cmake_policy(SET CMP0054 NEW) # Only interpret if() arguments as variables or keywords when unquoted.
    cmake_policy(SET CMP0076 NEW) # target_sources() command converts relative paths to absolute.
endif(COMMAND cmake_policy)

project(subdeval VERSION 0.0.1)


add_library(subdeval)
file(GLOB_RECURSE SOURCES_ "lib/subdeval/*.h" "lib/subdeval/*.C")
target_sources(subdeval PRIVATE ${SOURCES_})
target_include_directories(subdeval PUBLIC "lib")

add_executable(subdeval_example)
file(GLOB_RECURSE SOURCES_ "test/*.h" "test/*.C")
target_sources(subdeval_example PRIVATE ${SOURCES_})
target_link_libraries(subdeval_example PRIVATE subdeval)
