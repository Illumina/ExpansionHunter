add_library(locus_spec
        CnvLocusSpec.hh CnvLocusSpec.cpp
        GraphLocusSpec.hh GraphLocusSpec.cpp
        LocusSpec.hh LocusSpec.cpp
        VariantSpec.hh VariantSpec.cpp
        )
target_link_libraries(locus_spec common graphtools)