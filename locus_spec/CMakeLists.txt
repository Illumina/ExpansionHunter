add_library(locus_spec
        CnvLocusSpec.hh CnvLocusSpec.cpp
        GraphLocusSpec.hh GraphLocusSpec.cpp
        LocusSpec.hh LocusSpec.cpp
        GraphBlueprint.hh GraphBlueprint.cpp
        RegionGraph.hh RegionGraph.cpp
        GraphLocusEncoding.hh GraphLocusEncoding.cpp)
target_link_libraries(locus_spec common graphtools)