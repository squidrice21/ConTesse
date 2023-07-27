## ================== [shor_van_wyck] ============
set(SHOR_SOURCE_DIR "${PROJECT_SOURCE_DIR}/external/shor_van_wyck")

add_library(shor_van_wyck STATIC)
target_sources(shor_van_wyck PRIVATE
  ${SHOR_SOURCE_DIR}/shor.cpp
  ${SHOR_SOURCE_DIR}/is_simple_polygon.cpp
  ${SHOR_SOURCE_DIR}/plot.cpp
)
target_link_libraries(shor_van_wyck PUBLIC
  igl::igl
)
target_include_directories(shor_van_wyck PUBLIC
  ${SHOR_SOURCE_DIR}
)
add_library(shor_van_wyck::shor_van_wyck ALIAS shor_van_wyck)
set_target_properties(shor_van_wyck PROPERTIES FOLDER third_party)
