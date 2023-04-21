// Run ellipsoid random sample calculation.
// Author: Callum Marples

#include "run_tests.hpp"
#include "random_points_surface.hpp"
#include "ellipsoid_shape.hpp"
#include "ellipsoid_patches.hpp"
#include "contact_distributions.hpp"
#include "rng.hpp"
#include "vectors.hpp"
#include <ctime>
#include <cmath>
#include <memory>
#include <vector>
#include <array>

using namespace std;

template <typename T>
vector<Vec3> PerformSampleLoop(T& Sampler, int N) {
    vector<Vec3> V;
    V.resize(N);
    for (int i=0; i!=N; ++i) {
        V[i] = Sampler.RandomPoint();
    }
    return V;
}

void WriteSample(vector<Vec3> V, string file, int N) {
    ofstream File(file);
    if (File.is_open()) {
        File << fixed << setprecision(15);
        for (int i=0; i!=N; ++i) {
            File << V[i][0] << "," << V[i][1] << "," << V[i][2] << "\n";
        }
        File.close();
    }


}


// Generate a small sample of random points on a sphere/ellipsoid and output each individual point.
// This data is used to illustrate the non-uniformity of the naive anisotropic scaling method.
void RunSamples(int& N, string& rngString) {
    cout << "\nRun loop of N = " << N << " random samples, recording all points.\n\n";

    double r1 = 1.0;
    double a3 = 3.0;
    double b3 = 3.0;
    double c3 = 1.5;

    // Initialise and 'warm up' random number generator
    shared_ptr<RNG> rng;
    unsigned long seed = 123456789;
    if (rngString == "laggedfib") {
        rng = shared_ptr<RNG> ( new RNG_lf(seed) );
    }
    else if (rngString == "yarn") {
        rng = shared_ptr<RNG> ( new RNG_yn(seed) );
    }
    else {
        rng = shared_ptr<RNG> ( new RNG_mt(seed) );
    }
    double x;
	for (int i=0; i!=N; ++i) {
		x = rng->RandNum();
	}

    // Sphere
    SphereMarsaglia R1(r1, rng);
    vector<Vec3> vSphere = PerformSampleLoop(R1, N);

    // Ellipsoid Naive
    EllipsoidNaive R2(a3, b3, c3, rng);
    vector<Vec3> vNaive = PerformSampleLoop(R2, N);

    // Ellipsoid Correct
    EllipsoidGradRej R3(a3, b3, c3, rng);
    vector<Vec3> vCorrect = PerformSampleLoop(R3, N);

    // Write to files
    string file1 = "data/sample_sphere_" + to_string(N) + ".csv";
    string file2 = "data/sample_naive_" + to_string(N) + ".csv";
    string file3 = "data/sample_correct_" + to_string(N) + ".csv";
    WriteSample(vSphere, file1, N);
    WriteSample(vNaive, file2, N);
    WriteSample(vCorrect, file3, N);

}
