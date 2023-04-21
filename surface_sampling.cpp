// Call the relevant sampling "Run" routine.
// Author: Callum Marples

#include "random_points_surface.hpp"
#include "rng.hpp"
#include "vectors.hpp"
#include "run_tests.hpp"
#include <iostream>

using namespace std;

int main() {

    string program;
    cout << "Specify program: ";
    cin >> program;

    if (program == "time") {
        // Number of samples
        int N;
        cout << "\nNumber of points to sample: ";
        cin >> N;
        // Random number generator
        string rngString;
        cout << "\nRandom number generator: ";
        cin >> rngString;
        RunTiming(N, rngString);
    }
    else if (program == "dist") {
        // Number of samples
        int N;
        cout << "\nNumber of points to sample: ";
        cin >> N;
        // Random number generator
        string rngString;
        cout << "\nRandom number generator: ";
        cin >> rngString;
        RunDistributions(N, rngString);
    }
    else if (program == "samp") {
        // Number of samples
        int N;
        cout << "\nNumber of points to sample: ";
        cin >> N;
        // Random number generator
        string rngString;
        cout << "\nRandom number generator: ";
        cin >> rngString;
        RunSamples(N, rngString);
    }
    else if (program == "area") {
        RunAreas();
    }


	return 0;
}
