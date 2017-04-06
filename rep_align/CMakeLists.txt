add_library(rep_align rep_align.cc)
target_link_libraries(rep_align purity common)

add_subdirectory(unit_tests)