/* Routines for Numerical Integration
   Author: Callum Marples
*/

#ifndef INTEGRATION_DEF
#define INTEGRATION_DEF

#include <functional>

// Class implementing the trapezium rule.
// Based on 'Trapzd' from Numerical Recipes
class Trapz {
    private:
        std::function<double(double)> Fun;                                      // Function to be integrated ...
        double a;                                                               // and limits (a, b)
        double b;
        int n;                                                                  // Level of refinement
    public:
        double sum;                                                             // Current approximation of integral
        Trapz(std::function<double(double)> FunIn, double aLim, double bLim);
        void Refine();                                                          // Perform one refinement step
        // The actual integration is performed in the 'TrapeziumRule' function, defined outside this class.
        // The function takes F(x) and limits (a, b), creates a Trapz object and refines until meeting a tolerance.
};
double TrapeziumRule(std::function<double(double)> Fun, double a, double b, double eps=1.0e-10, int jMax=20);
double RombergIntegration(std::function<double(double)> Fun, double a, double b);

// Elliptic integrals. Used for surface area calculations.
double CarlsonRF(double x, double y, double z);
double CarlsonRC(double x, double y);
double CarlsonRD(double x, double y, double z);
double CarlsonRJ(double x, double y, double z, double p);

double Legendre1st(double phi, double k);
double Legendre2nd(double phi, double k);
double Legendre1st(double k);
double Legendre2nd(double k);

inline double Max(const double& x, const double& y) { if (x > y) { return x; } else { return y; } }
inline double SQR(const double& x) { return x*x; }

#endif // INTEGRATION_DEF
