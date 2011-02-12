cmake_minimum_required(VERSION 2.8)

set(SOURCES
    example.c
    example.h
)

add_library(example ${SOURCES})
target_link_libraries(example m) 
