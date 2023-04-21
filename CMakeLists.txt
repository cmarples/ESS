cmake_minimum_required(VERSION 3.16)

# set the project name
project(
    EllipsoidSurfaceSampling 
	VERSION 0.0.1
    DESCRIPTION "Test ellipsoid surface sampling algorithms"
    LANGUAGES CXX)
	
# TRNG library code
add_subdirectory(lib/trng)

# compiled library code
add_subdirectory(src)

# create executable
add_executable(run_sampling surface_sampling.cpp)
target_link_libraries(run_sampling PUBLIC sampling trng4) 