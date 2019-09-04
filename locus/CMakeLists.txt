add_library(locus STATIC
        Locus.hh Locus.cpp
        GraphLocus.hh GraphLocus.cpp StrFeature.cpp StrFeature.hh)

target_link_libraries(locus reads alignment filtering common)