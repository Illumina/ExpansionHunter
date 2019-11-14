file(GLOB SOURCES "*.cpp")
add_library(output ${SOURCES})
target_link_libraries(output common stats genotyping locus_spec)
