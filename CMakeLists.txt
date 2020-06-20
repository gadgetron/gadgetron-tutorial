cmake_minimum_required(VERSION 3.10)

project(GADGETRON-TUTORIAL LANGUAGES CXX C)

find_package(Gadgetron REQUIRED)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
add_library(gadgetron_tutorial SHARED accumulate_and_recon.cpp )
target_link_libraries(gadgetron_tutorial Gadgetron::gadgetron_core)