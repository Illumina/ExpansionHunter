add_library(filtering
        OrientationPredictor.hh OrientationPredictor.cpp
        BloomFilter.cpp BloomFilter.hh
        ${CMAKE_SOURCE_DIR}/thirdparty/murmur/MurmurHash3.h ${CMAKE_SOURCE_DIR}/thirdparty/murmur/MurmurHash3.cpp)
target_link_libraries(filtering common graphtools)

add_subdirectory(tests)
