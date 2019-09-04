add_library(locus STATIC
        Locus.hh Locus.cpp
        GraphLocus.hh GraphLocus.cpp)

target_link_libraries(locus reads alignment filtering common)