/* Ellipsoid Shape Class, "EllipsoidShape"
   Author: Callum Marples
*/

#ifndef SHAPE_DEF
#define SHAPE_DEF

#include <iostream>

const double PI = 3.141592653589793238463;

// This class represents an ellipsoid with semi-axis lengths a, b and c.
// Axes a, b and c are specified such that,
//
// x = a * cos(phi) * sin(theta)
// y = b * sin(phi) * sin(theta)
// z = c * cos(theta)
//
// where theta = [0,pi] and phi = [0,2*pi]
//
// i.e. theta is the angle measured from the c-axis
// and phi is the angle in the x-y plane, starting from x.
class EllipsoidShape {
    public:
        // Data
        double aAxis;
        double bAxis;
        double cAxis;
        // Constructors
        EllipsoidShape();
        EllipsoidShape(const double& r);
        EllipsoidShape(const double& a, const double& c);
        EllipsoidShape(const double& a, const double& b, const double& c);
        EllipsoidShape(const EllipsoidShape& otherEllipsoidShape);
        EllipsoidShape(const EllipsoidShape& otherEllipsoidShape, const double& scale);
        // Output
        friend std::ostream& operator<<(std::ostream& output, const EllipsoidShape& E);
        // Setters
        void SetA(double a) { aAxis = a; }
        void SetB(double b) { bAxis = b; }
        void SetC(double c) { cAxis = c; }
        void SetAxes(double a, double b, double c);
        void SetSkinRadius();
        // Getters
        double GetA() const { return aAxis; }
        double GetB() const { return bAxis; }
        double GetC() const { return cAxis; }
        double GetSkinRadius() { return *pSkinRadius; }
        bool IsSphere() const { return isSphere; }
        // Geometry
        double Volume() { return 4.0*PI*aAxis*bAxis*cAxis/3.0; }
    private:
        // Data
        bool isSphere;
        double* pSkinRadius;
};

#endif // SHAPE_DEF

