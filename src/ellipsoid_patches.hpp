/* Ellipsoid patch class, "EllipsoidPatches"
   Author: Callum Marples
*/

#ifndef PATCH_DEF
#define PATCH_DEF

#include <iostream>
#include <vector>

// Find index of closest value in vector v to number x, using the binary search algorithm.
int BinarySearch(std::vector<double>& v, const double& x);

// This class stores the theta and phi coordinate information for a grid representing the
// surface of a generic ellipsoid (i.e. without shape information). This includes information
// such as the number of patches, number of theta and phi values and the increment in each coordinate.
class EllipsoidPatches {
    public:
        // Constructors
        EllipsoidPatches();
        EllipsoidPatches(const int& nThetaIn, const int& nPhiIn);
        // Getters
        int GetNTheta() { return nTheta; }
        int GetNPhi() { return nPhi; }
        int GetNoPatches() { return noPatches; }
        double GetDeltaTheta() { return deltaTheta; }
        double GetDeltaPhi() { return deltaPhi; }
        // Resize an input vector
        void InitialiseVector(std::vector<double>& v) { v.resize(noPatches, 0.0); }
        void InitialiseVector(std::vector<int>& v) { v.resize(noPatches, 0); }
        // Get index of patch with theta index i and phi index j
        int PatchIndex(const int& thetaIndex, const int& phiIndex);
    protected:
        // Data
        int nTheta;
        int nPhi;
        int noPatches;
        double deltaTheta;
        double deltaPhi;
};

// This variant of the patches class includes lists of the theta and phi values used to
// define an actual grid. Each patch can be represented using either one or two integers 
// (i.e. patch index or theta-phi index). Methods are included to convert between these index types.
class PolarPatches : public EllipsoidPatches {
    public:
        // Constructors
        PolarPatches() : EllipsoidPatches() { InitialiseValues(); }
        PolarPatches(EllipsoidPatches& E) : EllipsoidPatches(E) { InitialiseValues(); }
        PolarPatches(const int& nThetaIn, const int& nPhiIn) : EllipsoidPatches(nThetaIn, nPhiIn) { InitialiseValues(); }
        // Construct theta and phi value arrays
        void InitialiseValues();
        // Find indices given theta and phi values
        int FindThetaIndex(const double& th) { return BinarySearch(thVals, th); }
        int FindPhiIndex(const double& ph);
    protected:
        std::vector<double> thVals;
        std::vector<double> phVals;
};

#endif // PATCH_DEF
