cmake_minimum_required(VERSION 3.10.0)
PROJECT(FixOverlapVCF)
SET(CMAKE_SHARED_LINKER_FLAGS  -lz -lm -lbz2 -llzma -pthread -lcurl)

set(CMAKE_CXX_STANDARD 11)
include_directories(htslib)
link_libraries(pthread curl)
add_custom_command(
        OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/htslib/libhts.a
        COMMAND make
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/htslib/
)
add_custom_target(lhts DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/htslib/libhts.a )
add_executable(FixOverlapVCF ${CMAKE_CURRENT_SOURCE_DIR}/src/FixOverlapVCF.cpp)
target_link_libraries(FixOverlapVCF ${CMAKE_CURRENT_SOURCE_DIR}/htslib/libhts.a ${CMAKE_SHARED_LINKER_FLAGS})
add_dependencies(FixOverlapVCF lhts )