/* Ellipsoid surface area routines.
   Author: Callum Marples
*/

#ifndef AREA_DEF
#define AREA_DEF

#include "ellipsoid_patches.hpp"
#include <vector>

// Calculate surface area using Legendre's formula
// with the elliptic integrals evaluated using the Carlson form.
double EllipsoidSurfaceArea(double a, double b, double c);

// Calculate approximate area using Knud Thomsen's formula
// (See http://www.numericana.com/answer/ellipsoid.htm#thomsen)
double EllipsoidApproxArea(double a, double b, double c, const double& p=1.6075);

// Calculate surface area numerically, using Romberg integration
double EllipsoidAreaNumerical(double a, double b, double c);

// Calculate patch areas on the surface of the ellipsoid
// The 8-fold symmetry of the triaxial ellipsoid is exploited to speed up the calculation
std::vector<double> EllipsoidPatchAreas(double a, double b, double c, EllipsoidPatches& pE);

// Calculate area of patch defined by corners (th0, th1, ph0, ph1)
double CalculatePatchArea(const double& a, const double& b, const double& c, const double& th0, const double& th1, double ph0, double ph1);

// Calculate surface area of a spheroid with repeated axis a and distinct axis c
double SpheroidSurfaceArea(double a, double c);

// Calculate area of spheroidal band defined by endpoints th0 and th1
double SpheroidBandArea(const double& a, const double& c, const double& th0, const double& th1);

#endif // AREA_DEF
