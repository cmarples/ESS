/* Implementation of ellipsoid surface area routines.
   Author: Callum Marples
*/

#include "surface_area.hpp"
#include "ellipsoid_patches.hpp"
#include "numerical_integration.hpp"
#include <cmath>
#include <functional>

using namespace std;
const double PI = 3.141592653589793238463;

//******************************************************************************
// Surface area of entire ellipsoid
//******************************************************************************

// Calculate surface area using Legendre's formula
// with the elliptic integrals evaluated using the Carlson form.
double EllipsoidSurfaceArea(double a, double b, double c) {
    double cos2Phi = c / a;        // cos(phi)
    double phi = acos(cos2Phi);
    cos2Phi = SQR(cos2Phi);        // cos^2(phi)
    double sin2Phi = 1.0 - cos2Phi;
    double k = sqrt( (a*a*(b*b - c*c)) / (b*b*(a*a - c*c)) );
    return 2.0*PI*( c*c + a*b*(sin2Phi*Legendre2nd(phi, k) + cos2Phi*Legendre1st(phi, k)) / sqrt(sin2Phi) );
}


// Calculate approximate area using Knud Thomsen's formula
// (See http://www.numericana.com/answer/ellipsoid.htm#thomsen)
double EllipsoidApproxArea(double a, double b, double c, const double& p) {
    return 4.0*PI*pow( ( pow(a, p)*pow(b, p) + pow(a, p)*pow(c, p) + pow(b, p)*pow(c, p) ) / 3.0, 1.0/p);
}


// Calculate surface area numerically
double EllipsoidAreaNumerical(double a, double b, double c) {
    double delta = 1.0 - c*c/(a*a);
    double epsln = 1.0 - c*c/(b*b);
    // Define integrand
    function<double(double)> I = [&delta, &epsln] (double s) {
        double k = sqrt(epsln*(1.0-s*s) / (1.0 - delta*s*s));
        return sqrt(1.0 - delta*s*s) * Legendre2nd(k);
    };
    return 8.0*a*b*RombergIntegration(I, 0.0, 1.0);
}


// Calculate surface area of a spheroid with repeated axis a and distinct axis c
double SpheroidSurfaceArea(double a, double c) {
    double q = 1.0 - a*a/(c*c);
    if (a < c) { // Prolate
        q = sqrt(q);
        return 2.0*PI*a*( a + c*asin(q)/q );
    }
    else {       // Oblate
        q = sqrt(-q);
        return 2.0*PI*a*( a + c*asinh(q)/q );
    }
}


// Calculate area of spheroidal band defined by endpoints th0 and th1
double SpheroidBandArea(const double& a, const double& c, const double& th0, const double& th1) {
    double q = 1.0 - a*a/(c*c);
    double z0 = c * cos(th0);
    double z1 = c * cos(th1);
    double rootq;
    double z0c = z0/c;
    double z1c = z1/c;
    if (a < c) { // Prolate
        rootq = sqrt(q);
        return PI*a*c*( (asin(rootq*z0c) - asin(rootq*z1c))/rootq + z0c*sqrt(1.0 - q*z0c*z0c) - z1c*sqrt(1.0 - q*z1c*z1c) );
    }
    else {       // Oblate
        rootq = sqrt(-q);
        return PI*a*c*( (asinh(rootq*z0c) - asinh(rootq*z1c))/rootq + z0c*sqrt(1.0 - q*z0c*z0c) - z1c*sqrt(1.0 - q*z1c*z1c) );
    }
}


//******************************************************************************
// Surface area of patches
//******************************************************************************

// Calculate patch areas on the surface of the ellipsoid
// The 8-fold symmetry of the triaxial ellipsoid is exploited to speed up the calculation
vector<double> EllipsoidPatchAreas(double a, double b, double c, EllipsoidPatches& pE) {
    // Define output vector
    vector<double> patchAreas;
    pE.InitialiseVector(patchAreas);
    // Get patch data
    int nTheta = pE.GetNTheta();
    int nPhi = pE.GetNPhi();
    int noPatches = pE.GetNoPatches();
    double deltaTheta = pE.GetDeltaTheta();
    double deltaPhi = pE.GetDeltaPhi();

    int thetaEnd = (nTheta - 1) / 2;
    int phiEnd = (nPhi - 2) / 4;
    bool isThetaOdd = true;
    if (nTheta % 2 == 0) { // nTheta even
        isThetaOdd = false;
    }
    int m = 0;
    bool isPhi4 = false;
    if (nPhi % 4 == 0) { // nPhi multiple of 4
        isPhi4 = true;
        m = (nPhi - 4) / 4; // Number of patches between a and b polar patches
        phiEnd = m + 1;
    }
    double th0, th1, ph0, ph1, lmda;
    double theta = 0.0;
    double phi;
    int k = 0;

    for (int i=0; i!=thetaEnd+1; ++i) {
        phi = 0.0;
        for (int j=0; j!=phiEnd+1; ++j) {
            // Calculate patch area
            if (k == 0) { // Polar Cap Areas
                patchAreas[0] = 4.0 * CalculatePatchArea(a, b, c, 0.0, 0.5*deltaTheta, 0.0, 0.5*PI);
                patchAreas[noPatches-1] = patchAreas[0];
                k += 1;
                break;
            }
            else {
                lmda = 1.0;
                // Theta limits
                if (i == thetaEnd && isThetaOdd) {
                    th0 = 0.5*(PI - deltaTheta);
                    th1 = 0.5*PI;
                    lmda *= 2.0;
                }
                else {
                    th0 = theta - 0.5*deltaTheta;
                    th1 = theta + 0.5*deltaTheta;
                }
                // Phi limits
                if (j == phiEnd && isPhi4) {
                    ph0 = 0.5*(PI - deltaPhi);
                    ph1 = 0.5*PI;
                    lmda *= 2.0;
                }
                else if (j == 0) {
                    ph0 = 0.0;
                    ph1 = 0.5*deltaPhi;
                    lmda *= 2.0;
                }
                else {
                    ph0 = phi - 0.5*deltaPhi;
                    ph1 = phi + 0.5*deltaPhi;
                }
                // Perform patch area calculation
                patchAreas[k] = lmda * CalculatePatchArea(a, b, c, th0, th1, ph0, ph1);
                // Reflect patch area to obtain area for equivalent patches
                if (isPhi4 == false) {
                    if (j == 0) {
                        patchAreas[pE.PatchIndex(i, nPhi/2)] = patchAreas[k];
                        if (!isThetaOdd || i != thetaEnd) {
                            patchAreas[pE.PatchIndex(nTheta-1-i, 0)] = patchAreas[k];
                            patchAreas[pE.PatchIndex(nTheta-1-i, nPhi/2)] = patchAreas[k];
                        }
                    }
                    else {
                        patchAreas[pE.PatchIndex(i, phiEnd+j)] = patchAreas[k];
                        patchAreas[pE.PatchIndex(i, nPhi/2+j)] = patchAreas[k];
                        patchAreas[pE.PatchIndex(i, nPhi/2+phiEnd+j)] = patchAreas[k];
                        if (!isThetaOdd || i != thetaEnd) {
                            patchAreas[pE.PatchIndex(nTheta-1-i, j)] = patchAreas[k];
                            patchAreas[pE.PatchIndex(nTheta-1-i, phiEnd+j)] = patchAreas[k];
                            patchAreas[pE.PatchIndex(nTheta-1-i, nPhi/2+j)] = patchAreas[k];
                            patchAreas[pE.PatchIndex(nTheta-1-i, nPhi/2+phiEnd+j)] = patchAreas[k];
                        }
                    }
                }
                else {
                    if (j == 0) {
                        patchAreas[pE.PatchIndex(i, nPhi/2)] = patchAreas[k];
                        if (!isThetaOdd || i != thetaEnd) {
                            patchAreas[pE.PatchIndex(nTheta-1-i, 0)] = patchAreas[k];
                            patchAreas[pE.PatchIndex(nTheta-1-i, nPhi/2)] = patchAreas[k];
                        }
                    }
                    else if (j == phiEnd) {
                        patchAreas[pE.PatchIndex(i, nPhi/2 + j)] = patchAreas[k];
                        if (!isThetaOdd || i != thetaEnd) {
                            patchAreas[pE.PatchIndex(nTheta-1-i, j)] = patchAreas[k];
                            patchAreas[pE.PatchIndex(nTheta-1-i, nPhi/2 + j)] = patchAreas[k];
                        }
                    }
                    else {
                        patchAreas[pE.PatchIndex(i, 2*phiEnd-j)] = patchAreas[k];
                        patchAreas[pE.PatchIndex(i, 2*phiEnd+j)] = patchAreas[k];
                        patchAreas[pE.PatchIndex(i, 4*phiEnd-j)] = patchAreas[k];
                        if (!isThetaOdd || i != thetaEnd) {
                            patchAreas[pE.PatchIndex(nTheta-1-i, j)] = patchAreas[k];
                            patchAreas[pE.PatchIndex(nTheta-1-i, 2*phiEnd-j)] = patchAreas[k];
                            patchAreas[pE.PatchIndex(nTheta-1-i, 2*phiEnd+j)] = patchAreas[k];
                            patchAreas[pE.PatchIndex(nTheta-1-i, 4*phiEnd-j)] = patchAreas[k];
                        }
                    }
                }
                // Advance counter k (pixel number)
                if (j == phiEnd) { k = pE.PatchIndex(i+1, 0); }
                else { k += 1; }
            }
            phi += deltaPhi;
        }
        theta += deltaTheta;
    }
    return patchAreas;
}

// Calculate area of patch defined by corners (th0, th1, ph0, ph1)
double CalculatePatchArea(const double& a, const double& b, const double& c, const double& th0, const double& th1, double ph0, double ph1) {
    if (a == c && b == c) {
        return (a*a*(ph1 - ph0)*(cos(th0) - cos(th1)));
    }
    else {
        double t0 = cos(th0);
        double t1 = cos(th1);
        double gFrac, kFrac, preFac;
        if (a > b) {
            gFrac = (b*b)/(c*c) - 1.0;
            kFrac = 1.0 - (b*b)/(a*a);
            preFac = a*c;
            ph0 = ph0 - 0.5*PI;
            ph1 = ph1 - 0.5*PI;
        }
        else {
            gFrac = (a*a)/(c*c) - 1.0;
            kFrac = 1.0 - (a*a)/(b*b);
            preFac = b*c;
        }
        // Define Integrand
        function<double(double)> I = [&gFrac, &kFrac, &ph0, &ph1] (double t) {
            double g = sqrt(1.0 + gFrac*t*t);
            double k = sqrt(kFrac*(1.0-t*t)) / g;
            return g*(Legendre2nd(ph1, k) - Legendre2nd(ph0, k));
        };
        return preFac*RombergIntegration(I, t1, t0);
    }
}
