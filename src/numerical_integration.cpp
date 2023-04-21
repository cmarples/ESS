/* Numerical Integration Implementation
   Author: Callum Marples
*/

#include "numerical_integration.hpp"
#include <functional>
#include <iostream>
#include <cmath>

using namespace std;
const double PI = 3.141592653589793238463;


// Integrate function Fun between limits a and b using Romberg integration
double RombergIntegration(std::function<double(double)> Fun, double a, double b) {
    double h = b - a;
    double eps = 1.0e-10;
    double TCurr[25]; // Stores values in Romberg table ...
    double TPrev[25]; // with jMax = 25 hardcoded
    bool conFlag = false; // Used to record whether convergence has been reached
    Trapz T(Fun, a, b);
    T.Refine();
    TPrev[0] = 0.5*h*T.sum; // Using endpoints only to give T(1,1)
    for (int j=1; j<=25; ++j) {
        T.Refine();
        TCurr[0] = T.sum; // Calculate T(j+1, 1) by refinement of the trapezium rule
        for (int k=2; k<=j; ++k) {
            // Calculate T(j+1, k) using Richardson extrapolation
            TCurr[k-1] = TCurr[k-2] + (TCurr[k-2] - TPrev[k-2]) / (pow(4.0, j-1.0) - 1.0);
        }
        // Check for convergence (two steps in a row must agree within the tolerance)
        if (fabs(TCurr[j-1] - TPrev[j-2]) < eps) {
            if (conFlag) { return TCurr[j-1]; }
            else { conFlag = true; }
        }
        else if (conFlag) { conFlag = false; }
        h *= 0.5;
        // Reassign TCurr to TPrev
        for (int i=0; i<j; ++i) {
            TPrev[i] = TCurr[i];
        }
    }
    cout << "Maximum number of trapezium rule refinements reached!" << endl;
    return TCurr[24];
}



//************************************************************************************************************
// Trapezium Rule Implementation
//************************************************************************************************************

// Trapz constructor
Trapz::Trapz(std::function<double(double)> FunIn, double aLim, double bLim) {
    Fun = FunIn;
    a = aLim;
    b = bLim;
    n = 0;
    sum = 0.0;
}

// Perform one refinement step of the extended trapezium rule.
// This increases n by one, takes the new points and recalculates the sum approximating the integral.
// Based on the Numerical Recipes routine 'trapzd'.
void Trapz::Refine() {
    n += 1;
    if (n == 1) {
        sum = 0.5*(b - a)*(Fun(a) + Fun(b));
    }
    else {
        int it, j;
        double x, tnm, sumAdd, del;
        for (it=1, j=1; j<n-1; j++) { it <<=1; } // Number of new points (is a power of 2)
        tnm = it;
        del = (b-a) / tnm;                       // New spacing
        x = a + 0.5*del;
        for (sumAdd=0.0, j=1; j<=it; ++j, x+=del) { sumAdd += Fun(x); }
        sum = 0.5*(sum + (b-a)*sumAdd/tnm);
    }
}

// Use a Trapz object to integrate function Fun between limits a and b.
double TrapeziumRule(std::function<double(double)> Fun, double a, double b, double eps, int jMax) {
    Trapz T(Fun, a, b);
    T.Refine();
    double oldSum;
    for (int j=0; j<jMax; ++j) {
        oldSum = T.sum;
        T.Refine();
        if (j > 5) {
            if (fabs(T.sum - oldSum) < eps*fabs(oldSum) || (T.sum == 0.0 && oldSum == 0.0)) {
                return T.sum;
            }
        }
    }
    cout << "Maximum number of trapezium rule refinements reached!" << endl;
    return T.sum;
}


//************************************************************************************************************
// Elliptic Integrals
//************************************************************************************************************

// Function to numerically evaluate Carlson's integral of the first kind.
double CarlsonRF(double x, double y, double z) {
    double r = 0.0025; // Used to define tolerance
    double third = 1.0 / 3.0;
    double rtx, rty, rtz, lmda, delx, dely, delz, A;
    do {
        rtx = sqrt(x);
        rty = sqrt(y);
        rtz = sqrt(z);
        lmda = rtx*(rty + rtz) + rty*rtz;
        x = 0.25 * (x + lmda);
        y = 0.25 * (y + lmda);
        z = 0.25 * (z + lmda);
        A = third * (x+y+z);
        delx = (A - x) / A;
        dely = (A - y) / A;
        delz = (A - z) / A;
    } while (Max(Max(fabs(delx), fabs(dely)), fabs(delz)) > r);
    double E2 = delx*dely - delz*delz;
    double E3 = delx*dely*delz;
    return (1.0 - E2/10.0 + E3/14.0 + E2*E2/24.0 - 3.0*E2*E3/44.0) / sqrt(A);
}


// Function to numerically evaluate the degenerate Carlson integral of the first kind.
double CarlsonRC(double x, double y) {
    double r = 0.0012;
    double third = 1.0 / 3.0;
    double lmda, s, A;
    do {
        lmda = 2.0*sqrt(x)*sqrt(y) + y;
        x = 0.25 * (x + lmda);
        y = 0.25 * (y + lmda);
        A = third * (x+y+y);
        s = (y - A) / A;
    } while (fabs(s) > r);
    return (1.0 + s*s*(0.3 + s*(1.0/7.0 + s*(0.375 + s*9.0/22.0)))) / sqrt(A);
}


// Function to numerically evaluate Carlson's integral of the second kind.
const double D1 = 3.0/14.0;
const double D2 = 1.0/6.0;
const double D3 = 9.0/22.0;
const double D4 = 3.0/26.0;
const double D5 = 0.25*D3;
const double D6 = 1.5*D4;
double CarlsonRD(double x, double y, double z) {
    double rtx, rty, rtz, lmda, delx, dely, delz, A;
    double r = 0.0015;
    double sum = 0.0;
    double fac = 1.0;
    do {
        rtx = sqrt(x);
        rty = sqrt(y);
        rtz = sqrt(z);
        lmda = rtx*(rty + rtz) + rty*rtz;
        sum += fac/(rtz*(z+lmda));
        fac = 0.25*fac;
        x = 0.25 * (x + lmda);
        y = 0.25 * (y + lmda);
        z = 0.25 * (z + lmda);
        A = 0.2 * (x + y + 3.0*z);
        delx = (A - x) / A;
        dely = (A - y) / A;
        delz = (A - z) / A;
    } while (Max(Max(fabs(delx), fabs(dely)), fabs(delz)) > r);
    double ea = delx*dely;
    double eb = delz*delz;
    double ec = ea - eb;
    double ed = ea - 6.0*eb;
    double ee = ed + ec + ec;
    return 3.0*sum + fac*( 1.0+ed*(-D1 + D5*ed - D6*delz*ee) + delz*(D2*ee + delz*(-D3*ec + delz*D4*ea)) ) / (A*sqrt(A));
}


// Function to numerically evaluate Carlson's integral of the third kind.
const double C1 = 3.0/14.0;
const double C2 = 1.0/3.0;
const double C3 = 3.0/22.0;
const double C4 = 3.0/26.0;
const double C5 = 0.75*C3;
const double C6 = 1.5*C4;
const double C7 = 0.5*C2;
const double C8 = C3 + C3;
double CarlsonRJ(double x, double y, double z, double p) {
    double rtx, rty, rtz, lmda, delx, dely, delz, delp, A, alpha, beta;
    double r = 0.0015;
    double sum = 0.0;
    double fac = 1.0;
    do {
        rtx = sqrt(x);
        rty = sqrt(y);
        rtz = sqrt(z);
        lmda = rtx*(rty + rtz) + rty*rtz;
        alpha = SQR(p*(rtx + rty + rtz) + rtx*rty*rtz);
        beta = p*SQR(p + lmda);
        sum += fac*CarlsonRC(alpha, beta);
        fac = 0.25*fac;
        x = 0.25 * (x + lmda);
        y = 0.25 * (y + lmda);
        z = 0.25 * (z + lmda);
        p = 0.25 * (p + lmda);
        A = 0.2 * (x + y + z + p + p);
        delx = (A - x) / A;
        dely = (A - y) / A;
        delz = (A - z) / A;
        delp = (A - p) / A;
    } while (Max(Max(Max(fabs(delx), fabs(dely)), fabs(delz)), fabs(delp)) > r);
    double ea = delx*(dely + delz) + dely*delz;
    double eb = delx*dely*delz;
    double ec = delp*delp;
    double ed = ea - 3.0*ec;
    double ee = eb + 2.0*delp*(ea - ec);
    return 3.0*sum + fac*( 1.0+ed*(-C1 + C5*ed - C6*ee) + eb*(C7 + delp*(-C8 + delp*C4))
           + delp*ea*(C2 - delp*C3) - C2*delp*ec ) / (A*sqrt(A));
}



// Function to evaluate the Legendre incomplete integral of the first kind,
// using Carlson integral R_F
double Legendre1st(double phi, double k) {
    double sinPhi = sin(phi);
    double cosPhi = cos(phi);
    return sinPhi*CarlsonRF(cosPhi*cosPhi, (1.0 - sinPhi*k)*(1.0 + sinPhi*k), 1.0);
}

// Overloaded form of Legendre1st to handle the complete integral (only input k)
double Legendre1st(double k) { return CarlsonRF(0.0, 1.0 - k*k, 1.0); }


// Function to evaluate the Legendre incomplete integral of the second kind,
// using Carlson integrals R_F and R_D (degenerate R_F)
double Legendre2nd(double phi, double k) {
    double sinPhi = sin(phi);
    double cos2Phi = 1.0 - sinPhi*sinPhi;
    double p = SQR(k*sinPhi);
    double q = 1 - p;
    return sinPhi * (CarlsonRF(cos2Phi, q, 1.0) - p*CarlsonRD(cos2Phi, q, 1.0)/3.0);
}

// Overloaded form of Legendre2nd to handle the complete integral (only input k)
double Legendre2nd(double k) { return CarlsonRF(0.0, 1.0 - k*k, 1.0) - k*k*CarlsonRD(0.0, 1.0 - k*k, 1.0)/3.0; }
