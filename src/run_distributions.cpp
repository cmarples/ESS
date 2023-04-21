// Run random contact distribution calculation.
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
vector<int> PerformDistLoop(T& Sampler, int N, EllipsoidShape& E, PolarPatches& P) {
    ContactDistribution D(E, P);
    Vec3 v;
    for (int i=0; i!=N; ++i) {
        v = Sampler.RandomPoint();
        D.AddContact(v);
    }
    return D.GetContactDistribution();
}

template <typename T>
vector<int> PerformDistLoopCarts(T& Sampler, int N, EllipsoidShape& E, PolarPatches& P) {
    ContactDistribution D(E, P);
    Vec3 p;
    Polars v;
    for (int i=0; i!=N; ++i) {
        p = Sampler.RandomPoint();
        v = Sampler.Carts2Pols(p);
        D.AddContact(v.th, v.ph);
    }
    return D.GetContactDistribution();
}

template <typename T>
vector<int> PerformDistLoopPols(T& Sampler, int N, EllipsoidShape& E, PolarPatches& P) {
    ContactDistribution D(E, P);
    Polars v;
    for (int i=0; i!=N; ++i) {
        v = Sampler.RandomPointPols();
        D.AddContact(v.th, v.ph);
    }
    return D.GetContactDistribution();
}

// For each sampling algorithm, generate a 'contact distribution' on the chosen
// shapes. Output the resulting distributions to files. These distributions are
// then tested for uniformity.
void RunDistributions(int& N, string& rngString) {

    array<string, 6> sAlgos {"Cubic Rejection", "Marsaglia      ", "Trig           ", "Gaussian       ", "Cook           ", "Area Rejection "};
    array<string, 8> eAlgos {"Naive          ", "Gradient Reject", "Grad Rej (Trig)", "Grad Rej (Pols)", "Area Rejection ", "Area Rej (Pols)", "Area Rej (Mer) ", "Generic        "};
    array<vector<int>, 6> sSamples;
    array<vector<int>, 8> oSamples;
    array<vector<int>, 8> pSamples;
    array<vector<int>, 8> eSamples;

    cout << "\nRun loop of N = " << N << " random samples, recording the contact distribution.\n\n";

    double r1 = 1.0;
    double a2 = 3.0;
    double c2 = 1.5;
    double a3 = 3.0;
    double b3 = 2.0;
    double c3 = 1.0;

    EllipsoidShape E1(r1);            // Sphere    (1, 1, 1)
    EllipsoidShape E2(a2, c2);        // Oblate Spheroid  (3, 3, 1.5)
    EllipsoidShape E3(c2, a2);        // Prolate Spheroid  (1.5, 1.5, 3)
    EllipsoidShape E4(a3, b3, c3);    // Ellipsoid (3, 2, 1)

    PolarPatches P(181, 360);
    int noPatches = P.GetNoPatches();

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

    // Sphere algorithms (take r = 1) //
    cout << "Sphere Samplers\n\n";

    for (int i=0; i!=6; ++i) {
        // Select sampler and run
        if (i == 0) {
            SphereRejection R(r1, rng);
            sSamples[i] = PerformDistLoop(R, N, E1, P);
        }
        else if (i == 1) {
            SphereMarsaglia R(r1, rng);
            sSamples[i] = PerformDistLoop(R, N, E1, P);
        }
        else if (i == 2) {
            SphereTrig R(r1, rng);
            sSamples[i] = PerformDistLoop(R, N, E1, P);
        }
        else if (i == 3) {
            SphereGaussian R(r1, rng);
            sSamples[i] = PerformDistLoop(R, N, E1, P);
        }
        else if (i == 4) {
            SphereCook R(r1, rng);
            sSamples[i] = PerformDistLoop(R, N, E1, P);
        }
        else if (i == 5) {
            SphereAreaRejection R(r1, rng);
            sSamples[i] = PerformDistLoop(R, N, E1, P);
        }
        cout << sAlgos[i] << " complete\n";
    }

    cout << "\nOblate Spheroid Samplers\n\n";
    // Spheroid algorithms (take a, b, c = 3, 3, 1.5) //
    for (int i=0; i!=8; ++i) {
        // Select sampler and run
        if (i == 0) {
            EllipsoidNaive R(a2, a2, c2, rng);
            oSamples[i] = PerformDistLoop(R, N, E2, P);
        }
        else if (i == 1) {
            EllipsoidGradRej R(a2, a2, c2, rng);
            oSamples[i] = PerformDistLoop(R, N, E2, P);
        }
        else if (i == 2) {
            EllipsoidGradTrig R(a2, a2, c2, rng);
            oSamples[i] = PerformDistLoop(R, N, E2, P);
        }
        else if (i == 3) {
            EllipsoidGradRej R(a2, a2, c2, rng);
            oSamples[i] = PerformDistLoopCarts(R, N, E2, P);
        }
        else if (i == 4) {
            SpheroidAreaRej R(a2, c2, rng);
            oSamples[i] = PerformDistLoop(R, N, E2, P);
        }
        else if (i == 5) {
            SpheroidAreaRej R(a2, c2, rng);
            oSamples[i] = PerformDistLoopPols(R, N, E2, P);
        }
        else if (i == 6) {
            SpheroidAreaMer R(a2, c2, rng);
            oSamples[i] = PerformDistLoop(R, N, E2, P);
        }
        else if (i == 7) {
            EllipsoidGeneric R(a2, a2, c2, rng);
            oSamples[i] = PerformDistLoop(R, N, E2, P);
        }
        cout << eAlgos[i] << " complete\n";
    }

    cout << "\nProlate Spheroid Samplers\n\n";
    // Spheroid algorithms (take a, b, c = 1.5, 1.5, 3) //
    for (int i=0; i!=8; ++i) {
        // Select sampler and run
        if (i == 0) {
            EllipsoidNaive R(c2, c2, a2, rng);
            pSamples[i] = PerformDistLoop(R, N, E3, P);
        }
        else if (i == 1) {
            EllipsoidGradRej R(c2, c2, a2, rng);
            pSamples[i] = PerformDistLoop(R, N, E3, P);
        }
        else if (i == 2) {
            EllipsoidGradTrig R(c2, c2, a2, rng);
            pSamples[i] = PerformDistLoop(R, N, E3, P);
        }
        else if (i == 3) {
            EllipsoidGradRej R(c2, c2, a2, rng);
            pSamples[i] = PerformDistLoopCarts(R, N, E3, P);
        }
        else if (i == 4) {
            SpheroidAreaRej R(c2, a2, rng);
            pSamples[i] = PerformDistLoop(R, N, E3, P);
        }
        else if (i == 5) {
            SpheroidAreaRej R(c2, a2, rng);
            pSamples[i] = PerformDistLoopPols(R, N, E3, P);
        }
        else if (i == 6) {
            SpheroidAreaMer R(c2, a2, rng);
            pSamples[i] = PerformDistLoop(R, N, E3, P);
        }
        else if (i == 7) {
            EllipsoidGeneric R(c2, c2, a2, rng);
            pSamples[i] = PerformDistLoop(R, N, E3, P);
        }
        cout << eAlgos[i] << " complete\n";
    }


    cout << "\nEllipsoid Samplers\n\n";
    // Spheroid algorithms (take a, b, c = 3, 2, 1) //
    for (int i=0; i!=8; ++i) {
        // Select sampler and run
        if (i == 0) {
            EllipsoidNaive R(a3, b3, c3, rng);
            eSamples[i] = PerformDistLoop(R, N, E4, P);
        }
        else if (i == 1) {
            EllipsoidGradRej R(a3, b3, c3, rng);
            eSamples[i] = PerformDistLoop(R, N, E4, P);
        }
        else if (i == 2) {
            EllipsoidGradTrig R(a3, b3, c3, rng);
            eSamples[i] = PerformDistLoop(R, N, E4, P);
        }
        else if (i == 3) {
            EllipsoidGradRej R(a3, b3, c3, rng);
            eSamples[i] = PerformDistLoopCarts(R, N, E4, P);
        }
        else if (i == 4) {
            EllipsoidAreaRej R(a3, b3, c3, rng);
            eSamples[i] = PerformDistLoop(R, N, E4, P);
        }
        else if (i == 5) {
            EllipsoidAreaRej R(a3, b3, c3, rng);
            eSamples[i] = PerformDistLoopPols(R, N, E4, P);
        }
        else if (i == 6) {
            EllipsoidAreaMer R(a3, b3, c3, rng);
            eSamples[i] = PerformDistLoop(R, N, E4, P);
        }
        else if (i == 7) {
            EllipsoidGeneric R(a3, b3, c3, rng);
            eSamples[i] = PerformDistLoop(R, N, E4, P);
        }
        cout << eAlgos[i] << " complete\n";
    }

    // Write information to file
    string file = "data/distributions_" + rngString + ".csv";
    ofstream File(file);
    if (File.is_open()) {
        File << fixed;
        File << "Sampling distributions from each algorithm\n";
        File << "N = " << N << "\n";
        File << "Sphere,Oblate,Prolate,Triaxial\n";
        for (int i=0; i!=6; ++i) { File << sAlgos[i] << ","; }
        for (int i=0; i!=8; ++i) { File << eAlgos[i] << ","; }
        for (int i=0; i!=8; ++i) { File << eAlgos[i] << ","; }
        for (int i=0; i!=8; ++i) { File << eAlgos[i] << ","; }
        File << "\n";
        for (int j=0; j!=noPatches; ++j) {
            for (int i=0; i!=6; ++i) { File << sSamples[i][j] << ","; }
            for (int i=0; i!=8; ++i) { File << oSamples[i][j] << ","; }
            for (int i=0; i!=8; ++i) { File << pSamples[i][j] << ","; }
            for (int i=0; i!=8; ++i) { File << eSamples[i][j] << ","; }
            File << "\n";
        }
        File.close();
    }

}
