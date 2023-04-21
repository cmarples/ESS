// Declarations for the routines used to generate data.
// Author: Callum Marples

#ifndef TIMINGDEF
#define TIMINGDEF

#include "random_points_surface.hpp"
#include "ellipsoid_shape.hpp"
#include "ellipsoid_patches.hpp"
#include "rng.hpp"
#include "vectors.hpp"
#include <iomanip>
#include <fstream>
#include <iostream>
#include <cmath>
#include <ctime>
#include <vector>

// Timing
void RunTiming(int& N, std::string& rngString);

template <typename T>
double PerformTimeLoop(T& Sampler, int N);

// Distributions
void RunDistributions(int& N, std::string& rngString);

// Sample
void RunSamples(int& N, std::string& rngString);

// Patch Areas
void RunAreas();

template <typename T>
std::vector<int> PerformDistLoop(T& Sampler, int N, EllipsoidShape& E, PolarPatches& P);


#endif // TIMINGDEF
