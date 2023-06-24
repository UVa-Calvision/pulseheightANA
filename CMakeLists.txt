cmake_minimum_required(VERSION 3.9)
project (pulsetest CXX)

# Look for ROOT. As long as an installation is available in the environment,
#this should succeed.
message(STATUS "Looking for ROOT")
find_package(ROOT REQUIRED COMPONENTS RIO)
message(STATUS "ROOT ${ROOT_VERSION} found at ${ROOT_BINDIR}")
include(${ROOT_USE_FILE})

add_executable(DRSPulseAna DRSPulseAna.cpp)
target_include_directories(DRSPulseAna PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})


# generate dictionaries and add them as a dependency of the executable
# (via the MODULE parameter):
ROOT_GENERATE_DICTIONARY(pulse_dict MODULE DRSPulseAna LINKDEF DRS_LinkDef.h)
target_link_libraries(DRSPulseAna ${ROOT_LIBRARIES})

### This version is NOT working yet
### I'm trying to generate a shared library to allow reading the
### the trees in python.
### https://root-forum.cern.ch/t/shared-object-with-cmake/18645/3
# or
# build a shared library first
#ROOT_GENERATE_DICTIONARY(DRS_Dict LINKDEF DRS_LinkDef.h)
#add_library(DRS SHARED DRS_Dict.cxx)
#target_link_libraries(pulsetest DRS ${ROOT_LIBRARIES})


