add_library(input
        CatalogLoading.hh CatalogLoading.cpp
        ParameterLoading.hh ParameterLoading.cpp
        SampleStats.hh SampleStats.cpp)
target_link_libraries(input common stats locus_spec)
add_subdirectory(tests)