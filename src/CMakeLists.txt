set(HEADER_LIST 
    "${PROJECT_SOURCE_DIR}/src/random_points_surface.hpp"
	"${PROJECT_SOURCE_DIR}/src/vectors.hpp"
	"${PROJECT_SOURCE_DIR}/src/rng.hpp"
	"${PROJECT_SOURCE_DIR}/src/ellipsoid_shape.hpp"
	"${PROJECT_SOURCE_DIR}/src/ellipsoid_patches.hpp"
	"${PROJECT_SOURCE_DIR}/src/contact_distributions.hpp"
	"${PROJECT_SOURCE_DIR}/src/numerical_integration.hpp"
	"${PROJECT_SOURCE_DIR}/src/surface_area.hpp"
	"${PROJECT_SOURCE_DIR}/src/run_tests.hpp")

# create library
add_library(sampling STATIC 
    random_points_surface.cpp
	ziggurat.cpp
	ellipsoid_shape.cpp
	ellipsoid_patches.cpp
	contact_distributions.cpp
	numerical_integration.cpp
	surface_area.cpp
	run_timing.cpp
	run_distributions.cpp
	run_areas.cpp
	run_sample.cpp
	${HEADER_LIST})

# include directory needed for building and for running
target_include_directories(sampling PUBLIC ${PROJECT_SOURCE_DIR}/src ${PROJECT_SOURCE_DIR}/lib)
target_compile_features(sampling PUBLIC cxx_std_11)