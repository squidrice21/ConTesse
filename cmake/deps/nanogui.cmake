if(TARGET nanogui::nanogui)
    return()
endif()

if (NOT CMAKE_SIZEOF_VOID_P EQUAL 8)
  message(FATAL_ERROR "32 bit system is not supported.")
endif()

include(FetchContent)
FetchContent_Declare(
  nanogui
  GIT_REPOSITORY https://github.com/wjakob/nanogui.git
  GIT_TAG        e9ec8a1a9861cf578d9c6e85a6420080aa715c03
)
message(STATUS "contour-tesselation: creating target 'nanogui::nanogui'")
FetchContent_GetProperties(nanogui)
if(NOT nanogui_POPULATED)
  FetchContent_Populate(nanogui)
endif()
message(STATUS "contour-tesselation: Finished cloning nanogui")

add_library(nanogui)

# ========================[  Include dir  ] ===========================

target_include_directories(nanogui PUBLIC ${nanogui_SOURCE_DIR}/include)

# ========================[  Sources  ] ===========================

target_sources(nanogui PRIVATE 
  # Fonts etc.
  ${nanogui_SOURCE_DIR}/include/nanogui/glutil.h ${nanogui_SOURCE_DIR}/src/glutil.cpp
  ${nanogui_SOURCE_DIR}/include/nanogui/common.h ${nanogui_SOURCE_DIR}/src/common.cpp
  ${nanogui_SOURCE_DIR}/include/nanogui/widget.h ${nanogui_SOURCE_DIR}/src/widget.cpp
  ${nanogui_SOURCE_DIR}/include/nanogui/theme.h ${nanogui_SOURCE_DIR}/src/theme.cpp
  ${nanogui_SOURCE_DIR}/include/nanogui/layout.h ${nanogui_SOURCE_DIR}/src/layout.cpp
  ${nanogui_SOURCE_DIR}/include/nanogui/screen.h ${nanogui_SOURCE_DIR}/src/screen.cpp
  ${nanogui_SOURCE_DIR}/include/nanogui/label.h ${nanogui_SOURCE_DIR}/src/label.cpp
  ${nanogui_SOURCE_DIR}/include/nanogui/window.h ${nanogui_SOURCE_DIR}/src/window.cpp
  ${nanogui_SOURCE_DIR}/include/nanogui/popup.h ${nanogui_SOURCE_DIR}/src/popup.cpp
  ${nanogui_SOURCE_DIR}/include/nanogui/checkbox.h ${nanogui_SOURCE_DIR}/src/checkbox.cpp
  ${nanogui_SOURCE_DIR}/include/nanogui/button.h ${nanogui_SOURCE_DIR}/src/button.cpp
  ${nanogui_SOURCE_DIR}/include/nanogui/popupbutton.h ${nanogui_SOURCE_DIR}/src/popupbutton.cpp
  ${nanogui_SOURCE_DIR}/include/nanogui/combobox.h ${nanogui_SOURCE_DIR}/src/combobox.cpp
  ${nanogui_SOURCE_DIR}/include/nanogui/progressbar.h ${nanogui_SOURCE_DIR}/src/progressbar.cpp
  ${nanogui_SOURCE_DIR}/include/nanogui/slider.h ${nanogui_SOURCE_DIR}/src/slider.cpp
  ${nanogui_SOURCE_DIR}/include/nanogui/messagedialog.h ${nanogui_SOURCE_DIR}/src/messagedialog.cpp
  ${nanogui_SOURCE_DIR}/include/nanogui/textbox.h ${nanogui_SOURCE_DIR}/src/textbox.cpp
  ${nanogui_SOURCE_DIR}/include/nanogui/imagepanel.h ${nanogui_SOURCE_DIR}/src/imagepanel.cpp
  ${nanogui_SOURCE_DIR}/include/nanogui/imageview.h ${nanogui_SOURCE_DIR}/src/imageview.cpp
  ${nanogui_SOURCE_DIR}/include/nanogui/vscrollpanel.h ${nanogui_SOURCE_DIR}/src/vscrollpanel.cpp
  ${nanogui_SOURCE_DIR}/include/nanogui/colorwheel.h ${nanogui_SOURCE_DIR}/src/colorwheel.cpp
  ${nanogui_SOURCE_DIR}/include/nanogui/colorpicker.h ${nanogui_SOURCE_DIR}/src/colorpicker.cpp
  ${nanogui_SOURCE_DIR}/include/nanogui/graph.h ${nanogui_SOURCE_DIR}/src/graph.cpp
  ${nanogui_SOURCE_DIR}/include/nanogui/stackedwidget.h ${nanogui_SOURCE_DIR}/src/stackedwidget.cpp
  ${nanogui_SOURCE_DIR}/include/nanogui/tabheader.h ${nanogui_SOURCE_DIR}/src/tabheader.cpp
  ${nanogui_SOURCE_DIR}/include/nanogui/tabwidget.h ${nanogui_SOURCE_DIR}/src/tabwidget.cpp
  ${nanogui_SOURCE_DIR}/include/nanogui/glcanvas.h ${nanogui_SOURCE_DIR}/src/glcanvas.cpp
  ${nanogui_SOURCE_DIR}/include/nanogui/formhelper.h
  ${nanogui_SOURCE_DIR}/include/nanogui/toolbutton.h
  ${nanogui_SOURCE_DIR}/include/nanogui/opengl.h
  ${nanogui_SOURCE_DIR}/include/nanogui/nanogui.h
  ${nanogui_SOURCE_DIR}/include/nanogui/serializer/core.h
  ${nanogui_SOURCE_DIR}/include/nanogui/serializer/opengl.h
  ${nanogui_SOURCE_DIR}/include/nanogui/serializer/sparse.h
  ${nanogui_SOURCE_DIR}/src/serializer.cpp
)

# ========================[  Deps  ] ===========================


target_compile_definitions(nanogui PUBLIC -DNANOGUI_GLAD)
include(deps/glad)
target_link_libraries(nanogui PUBLIC glad::glad)
include(deps/glfw)
target_link_libraries(nanogui PUBLIC glfw::glfw)
include(deps/eigen)
target_link_libraries(nanogui PUBLIC Eigen3::Eigen)
include(deps/nanovg)
target_link_libraries(nanogui PUBLIC nanovg::nanovg)

if (WIN32)
  target_link_libraries(nanogui PUBLIC opengl32)
elseif (APPLE)
  find_library(cocoa_library Cocoa)
  find_library(opengl_library OpenGL)
  find_library(corevideo_library CoreVideo)
  find_library(iokit_library IOKit)
  target_link_libraries(nanogui PUBLIC ${cocoa_library} ${opengl_library} ${corevideo_library} ${iokit_library})
  target_sources(nanogui PRIVATE ${nanogui_SOURCE_DIR}/src/darwin.mm)
  target_compile_options(nanogui PRIVATE "-fobjc-arc")
elseif(CMAKE_SYSTEM MATCHES "Linux")
  target_link_libraries(nanogui PUBLIC GL Xxf86vm Xrandr Xinerama Xcursor Xi X11 pthread )
  target_link_libraries(nanogui PUBLIC rt)
endif()

# ========================[ Warning and flags ] ===========================

if (MSVC)
  # Disable annoying MSVC warnings (all targets)
  target_compile_options(nanogui PRIVATE "/D _CRT_SECURE_NO_WARNINGS")
endif()
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  # Quench annoying deprecation warnings when compiling GLFW on OSX
  target_compile_options(nanogui PRIVATE "-Wno-deprecated-declarations")
endif()


# ========================[ Resources ] ===========================

# Run simple cmake converter to put font files into the data segment

# Glob up resource files
file(GLOB resources "${nanogui_SOURCE_DIR}/resources/*.ttf")

# Concatenate resource files into a comma separated string
string (REGEX REPLACE "([^\\]|^);" "\\1," resources_string "${resources}")
string (REGEX REPLACE "[\\](.)" "\\1" resources_string "${resources_string}")

# Create command line for running bin2c cmake script
set(bin2c_cmdline
 "-DOUTPUT_C=nanogui_resources.cpp"
  "-DOUTPUT_H=nanogui_resources.h"
  "-DINPUT_FILES=${resources_string}"
  -P "${nanogui_SOURCE_DIR}/resources/bin2c.cmake"
)

# Run bin2c on resource files
add_custom_command(
  OUTPUT nanogui_resources.cpp nanogui_resources.h
  COMMAND ${CMAKE_COMMAND} ARGS ${bin2c_cmdline}
  DEPENDS ${resources}
  COMMENT "Running bin2c"
  PRE_BUILD VERBATIM
)

# Needed to generated files
target_include_directories(nanogui PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
target_sources(nanogui PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/nanogui_resources.cpp)


# ========================[ Double colon alias ] ===========================

add_library(nanogui::nanogui ALIAS nanogui)
set_target_properties(nanogui PROPERTIES FOLDER third_party)


