add_library(purity purity.cc)

if(BUILD_TESTS)
    add_subdirectory(unit_tests)
endif()
