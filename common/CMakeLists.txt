add_library(common
        Common.hh Common.cpp
        CountTable.hh CountTable.cpp
        GenomicRegion.hh GenomicRegion.cpp
        HtsHelpers.hh HtsHelpers.cpp
        Parameters.hh Parameters.cpp
        Reference.hh Reference.cpp
        ReferenceContigInfo.hh ReferenceContigInfo.cpp
        WorkflowContext.hh WorkflowContext.cpp
        )
target_link_libraries(common region_spec ${htslib_static} ${zlib_static} ${Boost_LIBRARIES})
target_include_directories(common PUBLIC ${CMAKE_SOURCE_DIR}/thirdparty/spdlog/include)
add_subdirectory(tests)
