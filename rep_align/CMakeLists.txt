add_library(rep_align rep_align.cc)
include_directories(${purity_SOURCE_DIR})
target_link_libraries(rep_align purity)

if(BUILD_TESTS)
    add_subdirectory(unit_tests)
endif()
