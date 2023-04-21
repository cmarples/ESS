/* Class methods for the "EllipsoidShape" class.
   Author: Callum Marples
*/

#include "ellipsoid_shape.hpp"
#include <iostream>
#include <cassert>
#include <cmath>

using namespace std;

//******************************************************************************
// EllipsoidShape Methods
//******************************************************************************

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

// Overridden default constructor.
// Initialises as unit sphere.
EllipsoidShape::EllipsoidShape() {
    aAxis = 1.0;
    bAxis = 1.0;
    cAxis = 1.0;
    isSphere = true;
    SetSkinRadius();
}

// Constructor with 1 value, r, specified.
// Gives sphere of radius r.
EllipsoidShape::EllipsoidShape(const double& r) {
    assert(r > 0);
    aAxis = r;
    bAxis = r;
    cAxis = r;
    isSphere = true;
    SetSkinRadius();
}

// Constructor with 2 values, giving a spheroid with repeated axis a and distinct axis c
EllipsoidShape::EllipsoidShape(const double& a, const double& c) {
    assert(a > 0);
    assert(c > 0);
    aAxis = a;
    bAxis = a;
    cAxis = c;
    if (a == c) { isSphere = true; }
    else { isSphere = false; }
    SetSkinRadius();
}

// Constructor with 3 values, giving a triaxial ellipsoid of axes a, b and c
EllipsoidShape::EllipsoidShape(const double& a, const double& b, const double& c) {
    assert(a > 0 && b > 0 && c > 0);
    aAxis = a;
    bAxis = b;
    cAxis = c;
    if (a == b && a == c) { isSphere = true; }
    else { isSphere = false; }
    SetSkinRadius();
}

// Overridden copy constructor.
// Copies entries of otherEllipsoidShape into a new EllipsoidShape.
EllipsoidShape::EllipsoidShape(const EllipsoidShape& otherEllipsoidShape) {
    aAxis = otherEllipsoidShape.aAxis;
    bAxis = otherEllipsoidShape.bAxis;
    cAxis = otherEllipsoidShape.cAxis;
}

// Constructor with an ellipsoid shape and a double to yield an ellipsoid
// with axes scaled by specified factor.
EllipsoidShape::EllipsoidShape(const EllipsoidShape& otherEllipsoidShape, const double& scale) {
    aAxis = scale * otherEllipsoidShape.aAxis;
    bAxis = scale * otherEllipsoidShape.bAxis;
    cAxis = scale * otherEllipsoidShape.cAxis;
}

// Set axis lengths after construction
void EllipsoidShape::SetAxes(double a, double b, double c) {
    assert(a > 0 && b > 0 && c > 0);
    aAxis = a;
    bAxis = b;
    cAxis = c;
    if (a == b && a == c) { isSphere = true; }
    else { isSphere = false; }
    SetSkinRadius();
}

// Set skin radius (i.e. maximum axis length)
void EllipsoidShape::SetSkinRadius() {
    pSkinRadius = &aAxis;
    if (bAxis > *pSkinRadius) { pSkinRadius = &bAxis; }
    if (cAxis > *pSkinRadius) { pSkinRadius = &cAxis; }
}


// Print ellipsoid axes to screen
std::ostream& operator<<(std::ostream& output, const EllipsoidShape& E) {
    output << "(a, b, c) = (" << E.aAxis << ", " << E.bAxis << ", " << E.cAxis << ")";
    return output;
}
