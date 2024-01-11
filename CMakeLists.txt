cmake_minimum_required(VERSION 3.5)

include_directories(samtools/htslib-1.18)
link_directories(samtools/htslib-1.18)

project(split_spliced_reads LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(split_spliced_reads main.cpp)

find_package(Threads REQUIRED)
target_link_libraries (split_spliced_reads PRIVATE libhts.a  z m bz2 lzma curl crypto deflate Threads::Threads)
