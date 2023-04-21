// Run ellipsoid sampling timing calculation.
// Author: Callum Marples

#include "run_tests.hpp"
#include "random_points_surface.hpp"
#include "rng.hpp"
#include "vectors.hpp"
#include <ctime>
#include <cmath>
#include <memory>

using namespace std;

template <typename T>
double PerformTimeLoop(T& Sampler, int N) {
    Vec3 p;
    clock_t t = clock();
    for (int i=0; i!=N; ++i) {
        p = Sampler.RandomPoint();
    }
    t = clock() - t;
    return (double) t / CLOCKS_PER_SEC;
}

template <typename T>
double PerformTimeLoopCart(T& Sampler, int N) {
    Vec3 p;
    Polars P;
    clock_t t = clock();
    for (int i=0; i!=N; ++i) {
        p = Sampler.RandomPoint();
        P = Sampler.Carts2Pols(p);
    }
    t = clock() - t;
    return (double) t / CLOCKS_PER_SEC;
}

template <typename T>
double PerformTimeLoopPols(T& Sampler, int N) {
    Polars P;
    clock_t t = clock();
    for (int i=0; i!=N; ++i) {
        P = Sampler.RandomPointPols();
    }
    t = clock() - t;
    return (double) t / CLOCKS_PER_SEC;
}

// This routine times each sampling algorithm for each shape and a given random number generator.
// Acceptance rates are also calculated.
// The results are output to files.
void RunTiming(int& N, string& rngString) {

    array<string, 6> sAlgos {"Cubic Rejection", "Marsaglia      ", "Trig           ", "Gaussian       ", "Cook           ", "Area Rejection "};
    array<double, 6> sTimes;
    array<int, 6> sAccpt;
    array<double, 6> sRates;

    array<string, 8> eAlgos {"Naive          ", "Gradient Reject", "Grad Rej (Trig)", "Grad Rej (Pols)", "Area Rejection ", "Area Rej (Pols)", "Area Rej (Mer) ", "Generic        "};
    array<double, 8> eTimes;
    array<int, 8> eAccpt;
    array<double, 8> eRates;

    array<string, 8> oAlgos {"Naive          ", "Gradient Reject", "Grad Rej (Trig)", "Grad Rej (Pols)", "Area Rejection ", "Area Rej (Pols)", "Area Rej (Mer) ", "Generic        "};
    array<double, 8> oTimes;
    array<int, 8> oAccpt;
    array<double, 8> oRates;

    array<string, 8> pAlgos {"Naive          ", "Gradient Reject", "Grad Rej (Trig)", "Grad Rej (Pols)", "Area Rejection ", "Area Rej (Pols)", "Area Rej (Mer) ", "Generic        "};
    array<double, 8> pTimes;
    array<int, 8> pAccpt;
    array<double, 8> pRates;

    cout << "\nRun loop of N = " << N << " random samples, using <ctime> to measure run-time.\n\n";
    cout << "Algorithm       : Run-time    Accepted     Acceptance rate\n\n";

    // Initialise and 'warm up' random number generator
    shared_ptr<RNG> rng;
    unsigned long seed = 987654321;
    if (rngString == "laggedfib") {
        rng = shared_ptr<RNG> ( new RNG_lf(seed) );
    }
    else if (rngString == "yarn") {
        rng = shared_ptr<RNG> ( new RNG_yn(seed) );
    }
    else {
        rng = shared_ptr<RNG> ( new RNG_mt(seed) );
    }
    // Warm up rng and measure run-time
    double x, t_rng;
    clock_t t = clock();
	for (int i=0; i!=N; ++i) {
		x = rng->RandNum();
	}
    t = clock() - t;
    t_rng = (double) t / CLOCKS_PER_SEC;
    cout << "RNG             : " << t_rng << " s\n\n";
    cout << "Sphere Samplers\n\n";

    // Sphere algorithms (take r = 1) //
    for (int i=0; i!=6; ++i) {
        // Select sampler and run
        if (i == 0) {
            SphereRejection R(1.0, rng);
            sTimes[i] = PerformTimeLoop(R, N);
            sAccpt[i] = R.GetNoAttempts();
            sRates[i] = ((double)N)/((double)sAccpt[i]);
        }
        else if (i == 1) {
            SphereMarsaglia R(1.0, rng);
            sTimes[i] = PerformTimeLoop(R, N);
            sAccpt[i] = R.GetNoAttempts();
            sRates[i] = ((double)N)/((double)sAccpt[i]);
        }
        else if (i == 2) {
            SphereTrig R(1.0, rng);
            sTimes[i] = PerformTimeLoop(R, N);
            sAccpt[i] = R.GetNoAttempts();
            sRates[i] = ((double)N)/((double)sAccpt[i]);
        }
        else if (i == 3) {
            SphereGaussian R(1.0, rng);
            sTimes[i] = PerformTimeLoop(R, N);
            sAccpt[i] = R.GetNoAttempts();
            sRates[i] = ((double)N)/((double)sAccpt[i]);
        }
        else if (i == 4) {
            SphereCook R(1.0, rng);
            sTimes[i] = PerformTimeLoop(R, N);
            sAccpt[i] = R.GetNoAttempts();
            sRates[i] = ((double)N)/((double)sAccpt[i]);
        }
        else if (i == 5) {
            SphereAreaRejection R(1.0, rng);
            sTimes[i] = PerformTimeLoop(R, N);
            sAccpt[i] = R.GetNoAttempts();
            sRates[i] = ((double)N)/((double)sAccpt[i]);
        }
        cout << sAlgos[i] << " : " << sTimes[i] << " s     " << sAccpt[i] << "    " << sRates[i] << "\n";
    }


    cout << "\nOblate Spheroid Samplers\n\n";
    // Spheroid algorithms (take a, b, c = 3, 3, 1.5) //
    for (int i=0; i!=8; ++i) {
        // Select sampler and run
        if (i == 0) {
            EllipsoidNaive R(3.0, 3.0, 1.5, rng);
            oTimes[i] = PerformTimeLoop(R, N);
            oAccpt[i] = R.GetNoAttempts();
            oRates[i] = ((double)N)/((double)oAccpt[i]);
        }
        else if (i == 1) {
            EllipsoidGradRej R(3.0, 3.0, 1.5, rng);
            oTimes[i] = PerformTimeLoop(R, N);
            oAccpt[i] = R.GetNoAttempts();
            oRates[i] = ((double)N)/((double)oAccpt[i]);
        }
        else if (i == 2) {
            EllipsoidGradTrig R(3.0, 3.0, 1.5, rng);
            oTimes[i] = PerformTimeLoop(R, N);
            oAccpt[i] = R.GetNoAttempts();
            oRates[i] = ((double)N)/((double)oAccpt[i]);
        }
        else if (i == 3) {
            EllipsoidGradRej R(3.0, 3.0, 1.5, rng);
            oTimes[i] = PerformTimeLoopCart(R, N);
            oAccpt[i] = R.GetNoAttempts();
            oRates[i] = ((double)N)/((double)oAccpt[i]);
        }
        else if (i == 4) {
            SpheroidAreaRej R(3.0, 1.5, rng);
            oTimes[i] = PerformTimeLoop(R, N);
            oAccpt[i] = R.GetNoAttempts();
            oRates[i] = ((double)N)/((double)oAccpt[i]);
        }
        else if (i == 5) {
            SpheroidAreaRej R(3.0, 1.5, rng);
            oTimes[i] = PerformTimeLoopPols(R, N);
            oAccpt[i] = R.GetNoAttempts();
            oRates[i] = ((double)N)/((double)oAccpt[i]);
        }
        else if (i == 6) {
            SpheroidAreaMer R(3.0, 1.5, rng);
            oTimes[i] = PerformTimeLoop(R, N);
            oAccpt[i] = R.GetNoAttempts();
            oRates[i] = ((double)N)/((double)oAccpt[i]);
        }
        else if (i == 7) {
            EllipsoidGeneric R(3.0, 3.0, 1.5, rng);
            oTimes[i] = PerformTimeLoop(R, N);
            oAccpt[i] = R.GetNoAttempts();
            oRates[i] = ((double)N)/((double)oAccpt[i]);
        }
        cout << oAlgos[i] << " : " << oTimes[i] << " s     " << oAccpt[i] << "    " << oRates[i] << "\n";
    }

    cout << "\nProlate Spheroid Samplers\n\n";
    // Spheroid algorithms (take a, b, c = 1.5, 1.5, 3) //
    for (int i=0; i!=8; ++i) {
        // Select sampler and run
        if (i == 0) {
            EllipsoidNaive R(1.5, 1.5, 3.0, rng);
            pTimes[i] = PerformTimeLoop(R, N);
            pAccpt[i] = R.GetNoAttempts();
            pRates[i] = ((double)N)/((double)pAccpt[i]);
        }
        else if (i == 1) {
            EllipsoidGradRej R(1.5, 1.5, 3.0, rng);
            pTimes[i] = PerformTimeLoop(R, N);
            pAccpt[i] = R.GetNoAttempts();
            pRates[i] = ((double)N)/((double)pAccpt[i]);
        }
        else if (i == 2) {
            EllipsoidGradTrig R(1.5, 1.5, 3.0, rng);
            pTimes[i] = PerformTimeLoop(R, N);
            pAccpt[i] = R.GetNoAttempts();
            pRates[i] = ((double)N)/((double)pAccpt[i]);
        }
        else if (i == 3) {
            EllipsoidGradRej R(1.5, 1.5, 3.0, rng);
            pTimes[i] = PerformTimeLoopCart(R, N);
            pAccpt[i] = R.GetNoAttempts();
            pRates[i] = ((double)N)/((double)pAccpt[i]);
        }
        else if (i == 4) {
            SpheroidAreaRej R(1.5, 3.0, rng);
            pTimes[i] = PerformTimeLoop(R, N);
            pAccpt[i] = R.GetNoAttempts();
            pRates[i] = ((double)N)/((double)pAccpt[i]);
        }
        else if (i == 5) {
            SpheroidAreaRej R(1.5, 3.0, rng);
            pTimes[i] = PerformTimeLoopPols(R, N);
            pAccpt[i] = R.GetNoAttempts();
            pRates[i] = ((double)N)/((double)pAccpt[i]);
        }
        else if (i == 6) {
            SpheroidAreaMer R(1.5, 3.0, rng);
            pTimes[i] = PerformTimeLoop(R, N);
            pAccpt[i] = R.GetNoAttempts();
            pRates[i] = ((double)N)/((double)pAccpt[i]);
        }
        else if (i == 7) {
            EllipsoidGeneric R(1.5, 1.5, 3.0, rng);
            pTimes[i] = PerformTimeLoop(R, N);
            pAccpt[i] = R.GetNoAttempts();
            pRates[i] = ((double)N)/((double)pAccpt[i]);
        }
        cout << pAlgos[i] << " : " << pTimes[i] << " s     " << pAccpt[i] << "    " << pRates[i] << "\n";
    }


    cout << "\nEllipsoid Samplers\n\n";
    // Ellipsoid algorithms (take a, b, c = 3, 2, 1) //
    for (int i=0; i!=8; ++i) {
        // Select sampler and run
        if (i == 0) {
            EllipsoidNaive R(3.0, 2.0, 1.0, rng);
            eTimes[i] = PerformTimeLoop(R, N);
            eAccpt[i] = R.GetNoAttempts();
            eRates[i] = ((double)N)/((double)eAccpt[i]);
        }
        else if (i == 1) {
            EllipsoidGradRej R(3.0, 2.0, 1.0, rng);
            eTimes[i] = PerformTimeLoop(R, N);
            eAccpt[i] = R.GetNoAttempts();
            eRates[i] = ((double)N)/((double)eAccpt[i]);
        }
        else if (i == 2) {
            EllipsoidGradTrig R(3.0, 2.0, 1.0, rng);
            eTimes[i] = PerformTimeLoop(R, N);
            eAccpt[i] = R.GetNoAttempts();
            eRates[i] = ((double)N)/((double)eAccpt[i]);
        }
        else if (i == 3) {
            EllipsoidGradRej R(3.0, 2.0, 1.0, rng);
            eTimes[i] = PerformTimeLoopCart(R, N);
            eAccpt[i] = R.GetNoAttempts();
            eRates[i] = ((double)N)/((double)eAccpt[i]);
        }
        else if (i == 4) {
            EllipsoidAreaRej R(3.0, 2.0, 1.0, rng);
            eTimes[i] = PerformTimeLoop(R, N);
            eAccpt[i] = R.GetNoAttempts();
            eRates[i] = ((double)N)/((double)eAccpt[i]);
        }
        else if (i == 5) {
            EllipsoidAreaRej R(3.0, 2.0, 1.0, rng);
            eTimes[i] = PerformTimeLoopPols(R, N);
            eAccpt[i] = R.GetNoAttempts();
            eRates[i] = ((double)N)/((double)eAccpt[i]);
        }
        else if (i == 6) {
            EllipsoidAreaMer R(3.0, 2.0, 1.0, rng);
            eTimes[i] = PerformTimeLoop(R, N);
            eAccpt[i] = R.GetNoAttempts();
            eRates[i] = ((double)N)/((double)eAccpt[i]);
        }
        else if (i == 7) {
            EllipsoidGeneric R(3.0, 2.0, 1.0, rng);
            eTimes[i] = PerformTimeLoop(R, N);
            eAccpt[i] = R.GetNoAttempts();
            eRates[i] = ((double)N)/((double)eAccpt[i]);
        }
        cout << eAlgos[i] << " : " << eTimes[i] << " s     " << eAccpt[i] << "    " << eRates[i] << "\n";
    }


    // Write to file
    string file = "data/timing_" + rngString + ".csv";
    ofstream File(file);
    if (File.is_open())
    {
        File << fixed;
        File << "Algorithm," << "Run-time (seconds)," << "Acceptance rate\n";
        File << "N = " << N << "\n";
        File << "RNG," << setprecision(3) << t_rng << "\n";
        File << "\nSphere\n\n";
        for (int i=0; i!=6; ++i) {
            File << sAlgos[i] << "," << setprecision(3) << sTimes[i] << "," << setprecision(15) << sRates[i] << "\n";
        }
        File << "\nOblate Spheroid\n\n";
        for (int i=0; i!=8; ++i) {
            File << oAlgos[i] << "," << setprecision(3) << oTimes[i] << "," << setprecision(15) << oRates[i] << "\n";
        }
        File << "\nProlate Spheroid\n\n";
        for (int i=0; i!=8; ++i) {
            File << pAlgos[i] << "," << setprecision(3) << pTimes[i] << "," << setprecision(15) << pRates[i] << "\n";
        }
        File << "\nEllipsoid\n\n";
        for (int i=0; i!=8; ++i) {
            File << eAlgos[i] << "," << setprecision(3) << eTimes[i] << "," << setprecision(15) << eRates[i] << "\n";
        }
        File.close();
    }

}
