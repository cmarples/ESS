/* Class methods for the "EllipsoidPatches" and "PolarPatches" classes.
   Also includes a routine for the binary search algorithm.
   Author: Callum Marples
*/

#include "ellipsoid_patches.hpp"
#include <vector>
#include <cmath>

using namespace std;
const double PI = 3.141592653589793238463;

// Construct by specifying the number of theta and phi points to be used.
EllipsoidPatches::EllipsoidPatches(const int& nThetaIn, const int& nPhiIn) {
    // Number of vertices in grid
    nTheta = nThetaIn;
    nPhi = nPhiIn;
    if (nPhi % 2 != 0) { // Require even nPhi so the (half) patches have the same symmetry as the ellipsoid
        nPhi += 1;
        cout << "Require even nPhi so the (half) patches have the same symmetry as the ellipsoid" << endl;
        cout << "Therefore, nPhi has been increased by 1 to give nPhi = " << nPhi << endl;
    }
    noPatches = (nTheta - 2)*nPhi + 2;
    // Spacing (in radians) of theta and phi
    deltaTheta = PI / (nTheta - 1);
    deltaPhi = 2.0*PI / (nPhi);
}

// If no input given, use some (small) default values.
EllipsoidPatches::EllipsoidPatches() {
    nTheta = 19;
    nPhi = 18;
    noPatches = (nTheta - 2)*nPhi + 2;
    deltaTheta = PI / (nTheta - 1);
    deltaPhi = 2.0*PI / (nPhi);
}


// Get patch index from theta and phi indices
int EllipsoidPatches::PatchIndex(const int& thetaIndex, const int& phiIndex) {
	if (thetaIndex > 0 && thetaIndex < nTheta-1) { return ( 1 + phiIndex + nPhi*(thetaIndex-1) ); }
    else if (thetaIndex == 0) { return 0; }
    else { return noPatches-1; } // (thetaIndex == nTheta-1)
}

void PolarPatches::InitialiseValues() {
    thVals.resize(nTheta, 0.0);
    phVals.resize(nPhi+1, 0.0);
    for (int i=0; i!=nTheta; ++i) { thVals[i] = i * deltaTheta; }
    for (int i=0; i!=nPhi+1; ++i) { phVals[i] = i * deltaPhi; }
}

int PolarPatches::FindPhiIndex(const double& ph) {
    int index = BinarySearch(phVals, ph);
    if (index == nPhi) { return 0; }
    else { return index; }
}

// Find index of closest value in vector v to number x, using the binary search algorithm
int BinarySearch(std::vector<double>& v, const double& x) {
    // Initialise
    int L = 0;
    int R = static_cast<int>(v.size()) - 1;
    int M;
    // Find index L such that : v[L] <= x < v[L+1]
    while (L != R) {
        M = ceil( (L+R)/2.0 );
        if (x < v[M]) { R = M-1; }
        else          { L = M; }
    }
    // Output the index whose value is closest to x
    double dL = x - v[L];
    double dR = v[L+1] - x;
    if (dL < dR) { return L; }
    else { return L+1; }
}
