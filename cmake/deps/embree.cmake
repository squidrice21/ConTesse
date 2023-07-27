if(TARGET embree::embree)
    return()
endif()

include(FetchContent)
FetchContent_Declare(
    embree
    GIT_REPOSITORY https://github.com/embree/embree.git
    GIT_TAG        v3.9.0
    GIT_SHALLOW    TRUE
)

# Set Embree's default options
option(EMBREE_ISPC_SUPPORT "Build Embree with support for ISPC applications." OFF)
option(EMBREE_TUTORIALS    "Enable to build Embree tutorials"                 OFF)
option(EMBREE_STATIC_LIB   "Build Embree as a static library."                ON)
set(EMBREE_TESTING_INTENSITY 0         CACHE STRING "Intensity of testing (0 = no testing, 1 = verify and tutorials, 2 = light testing, 3 = intensive testing.")
set(EMBREE_TASKING_SYSTEM    "TBB"     CACHE STRING "Selects tasking system")
set(EMBREE_MAX_ISA           "DEFAULT" CACHE STRING "Selects highest ISA to support.")

# We want to compile Embree with TBB support, so we need to overwrite Embree's
# `find_package()` and provide variables. The following discussion provide some
# context on how to achieve this:
# - https://gitlab.kitware.com/cmake/cmake/issues/17735
# - https://crascit.com/2018/09/14/do-not-redefine-cmake-commands/
function(lagrange_import_embree)
    macro(ignore_package NAME)
        file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/${NAME}/${NAME}Config.cmake "")
        set(${NAME}_DIR ${CMAKE_CURRENT_BINARY_DIR}/${NAME} CACHE PATH "")
    endmacro()

    # Prefer Config mode before Module mode to prevent embree from loading its own FindTBB.cmake
    set(CMAKE_FIND_PACKAGE_PREFER_CONFIG TRUE)

    # Note that embree wants to be able to export() its target. If we use `TBB::tbb`
    # directly in `TBB_LIBRARIES` we will get an error. For now we hide tbb's target
    # from the exported dependencies of embree by going through an intermediate
    # IMPORTED library. If the client cares about where all this stuff is installed,
    # they will probably have their own way of installing embree/tbb.
    #
    # See relevant discussion here:
    # - https://gitlab.kitware.com/cmake/cmake/issues/17357
    # - https://gitlab.kitware.com/cmake/cmake/issues/15415
    ignore_package(TBB)
    include(deps/tbb)
    get_target_property(TBB_INCLUDE_DIRS TBB::tbb INTERFACE_INCLUDE_DIRECTORIES)
    add_library(lagrange_tbb INTERFACE IMPORTED)
    target_link_libraries(lagrange_tbb INTERFACE TBB::tbb)
    set(TBB_LIBRARIES lagrange_tbb)

    # Ready to include embree's atrocious CMake
    FetchContent_MakeAvailable(embree)

    # Now we need to do some juggling to propagate the include directory properties
    # along with the `embree` target
    add_library(embree::embree INTERFACE IMPORTED GLOBAL)
    target_include_directories(embree::embree SYSTEM INTERFACE ${embree_SOURCE_DIR}/include)
    target_link_libraries(embree::embree INTERFACE embree)
    add_library(embree::embree ALIAS embree)
endfunction()

# Call via a proper function in order to scope variables such as CMAKE_FIND_PACKAGE_PREFER_CONFIG and TBB_DIR
lagrange_import_embree()

# Some order for IDEs
set_target_properties(embree PROPERTIES FOLDER "third_party//embree")
set_target_properties(algorithms PROPERTIES FOLDER "third_party//embree")
set_target_properties(lexers PROPERTIES FOLDER "third_party//embree")
set_target_properties(math PROPERTIES FOLDER "third_party//embree")
set_target_properties(simd PROPERTIES FOLDER "third_party//embree")
set_target_properties(sys PROPERTIES FOLDER "third_party//embree")
set_target_properties(tasking PROPERTIES FOLDER "third_party//embree")
set_target_properties(uninstall PROPERTIES FOLDER "third_party//embree")