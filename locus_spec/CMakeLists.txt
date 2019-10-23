add_library(locus_spec
        CnvLocusSpecification.hh CnvLocusSpecification.cpp
        GraphLocusSpecification.hh GraphLocusSpecification.cpp
        LocusSpecification.hh LocusSpecification.cpp
        VariantSpecification.hh VariantSpecification.cpp
        )
target_link_libraries(locus_spec common graphtools)