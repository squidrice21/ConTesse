add_library(subdeval)

set(SUBDEVAL_SOURCE_DIR "${PROJECT_SOURCE_DIR}/external/subdeval")

file(GLOB_RECURSE SOURCES_ "${SUBDEVAL_SOURCE_DIR}/lib/subdeval/*.h" "${SUBDEVAL_SOURCE_DIR}/lib/subdeval/*.C")
target_sources(subdeval PRIVATE ${SOURCES_})
target_include_directories(subdeval PUBLIC "${SUBDEVAL_SOURCE_DIR}/lib")

add_library(subdeval::subdeval ALIAS subdeval)
set_target_properties(subdeval PROPERTIES FOLDER third_party)
