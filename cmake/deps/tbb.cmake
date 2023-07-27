####################################################################################################
# IMPORTANT
#
# This file defines a single ALIAS target `TBB::tbb`.
#
# Depending on the option CONTESS_USE_TBB_STATIC, this alias may point to either the dynamic version
# of TBB, or the static version. The official recommendation is to use the dynamic library:
#
# https://www.threadingbuildingblocks.org/faq/there-version-tbb-provides-statically-linked-libraries
# https://stackoverflow.com/questions/638278/how-to-statically-link-to-tbb
#
# For now we do not have a proper CMake workflow to deal with DLLs, so we default to tbb_static
####################################################################################################

if(TARGET TBB::tbb)
    return()
endif()

# Using wjakob's fork as it has a better cmake build system
# Change it back to intel's once they fix it
# https://github.com/intel/tbb/issues/6
include(FetchContent)
FetchContent_Declare(
    tbb
    GIT_REPOSITORY https://github.com/wjakob/tbb.git
    GIT_TAG 344fa84f34089681732a54f5def93a30a3056ab9
    GIT_SHALLOW FALSE
)

option(TBB_BUILD_SHARED          "Build TBB shared library" OFF)
option(TBB_BUILD_STATIC          "Build TBB static library" OFF)
option(TBB_BUILD_TBBMALLOC       "Build TBB malloc library" OFF)
option(TBB_BUILD_TBBMALLOC_PROXY "Build TBB malloc proxy library" OFF)
option(TBB_BUILD_TESTS           "Build TBB tests and enable testing infrastructure" OFF)
option(TBB_NO_DATE               "Do not save the configure date in the version string" ON)

# Make sure tbb or tbb_static is built, according to option
if(CONTESS_USE_TBB_STATIC)
    set(TBB_BUILD_STATIC ON CACHE BOOL "Build TBB static library" FORCE)
else()
    set(TBB_BUILD_SHARED ON CACHE BOOL "Build TBB shared library" FORCE)
endif()

FetchContent_MakeAvailable(tbb)

# Create an interface target following the upstream namespace `TBB::foo`.
# We use an imported target so that it doesn't appear in the export set of
# target (we assume that clients will have their own way of finding TBB).
add_library(TBB::tbb IMPORTED INTERFACE GLOBAL)
if(CONTESS_USE_TBB_STATIC)
    target_link_libraries(TBB::tbb INTERFACE tbb_static)
    get_target_property(TBB_INCLUDE_DIRS tbb_static INTERFACE_INCLUDE_DIRECTORIES)
else()
    target_link_libraries(TBB::tbb INTERFACE tbb)
    get_target_property(TBB_INCLUDE_DIRS tbb INTERFACE_INCLUDE_DIRECTORIES)
endif()
# Is this a bug in CMake? Transitive include dirs are not correctly propagated to
# the properties of the target TBB::tbb. Code seems to compile ok without it, but
# we need to be able to extract TBB's interface include directories for our embree.cmake
target_include_directories(TBB::tbb SYSTEM INTERFACE ${TBB_INCLUDE_DIRS})

foreach(name IN ITEMS tbb_def_files tbb_static tbb tbbmalloc tbbmalloc_static)
    if(TARGET ${name})
        set_target_properties(${name} PROPERTIES FOLDER third_party)
    endif()
endforeach()