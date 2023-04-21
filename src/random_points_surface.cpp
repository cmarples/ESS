/* Implementation of random points on the surface of a sphere/ellipsoid
   Author: Callum Marples
*/

#include "random_points_surface.hpp"
#include "rng.hpp"
#include "vectors.hpp"
#include <cmath>
#include <iostream>

using namespace std;
const double PI = 3.141592653589793238463;

//**********************************************************************************************//
// Generate uniformly random points on the surface of a sphere
//**********************************************************************************************//

// Using the cubic rejection method
Vec3 SphereRejection::RandomPoint() {
    Vec3 p;
    double rsq;
    do { // Generate random point in a cube of length 2, accepting it if it lies within a unit sphere
        p[0] = 1.0 - 2.0*rng->RandNum();
        p[1] = 1.0 - 2.0*rng->RandNum();
        p[2] = 1.0 - 2.0*rng->RandNum();
        rsq = p[0]*p[0] + p[1]*p[1] + p[2]*p[2];
        noAttempts += 1;
    } while (rsq > 1.0);
    rsq = 1 / sqrt(rsq); // This is now the norm of point p
    p *= rsq * radius;
    return p;
}


// Using Marsaglia's method
Vec3 SphereMarsaglia::RandomPoint() {
    Vec3 p;
    double rsq, r1, r2;
    do {
        r1 = 1.0 - 2.0*rng->RandNum();
        r2 = 1.0 - 2.0*rng->RandNum();
        rsq = r1*r1 + r2*r2;
        noAttempts += 1;
    } while (rsq >= 1.0);
    p[2] = radius * (1.0 - 2.0*rsq);
    rsq = 2.0*sqrt(1.0 - rsq);
    p[0] = radius * r1 * rsq;
    p[1] = radius * r2 * rsq;
    return p;
}


// Using Cook's method
Vec3 SphereCook::RandomPoint() {
    Vec3 p;
    double v1, v2, v3, v4, s;
    do {
        v1 = -1.0 + 2.0*rng->RandNum();
        v2 = -1.0 + 2.0*rng->RandNum();
        v3 = -1.0 + 2.0*rng->RandNum();
        v4 = -1.0 + 2.0*rng->RandNum();
        s = v1*v1 + v2*v2 + v3*v3 + v4*v4;
        noAttempts += 1;
    } while (s >= 1.0);
    s = 1.0 / s;
    p[0] = radius * 2.0 * s * (v2*v4 + v1*v3);
    p[1] = radius * 2.0 * s * (v3*v4 - v1*v2);
    p[2] = radius * s * (v1*v1 + v4*v4 - v2*v2 - v3*v3);
    return p;
}


// Using Gaussian random numbers
Vec3 SphereGaussian::RandomPoint() {
    Vec3 p;
    p[0] = rng->RandNor();
    p[1] = rng->RandNor();
    p[2] = rng->RandNor();
    Normalise(p);
    noAttempts += 1;
    return (radius * p);
}


// Using the trigonometric method
Vec3 SphereTrig::RandomPoint() {
    Vec3 p;
    double u = 1.0 - 2.0*rng->RandNum();
    double ph = 2.0*PI*rng->RandNum();
    p[2] = radius * u;
    u = radius * sqrt(1.0 - u*u);
    double cosPhi = cos(ph);
    p[0] = u * cosPhi;
    p[1] = u * sqrt(1.0 - cosPhi*cosPhi);
    if (ph > PI) { p[1] *= -1.0; }
    noAttempts += 1;
    return p;
}


// Using the area element rejection method
Vec3 SphereAreaRejection::RandomPoint() {
    Vec3 p;
    bool accept = false;
    double th, ph, sinTh, u;
    ph = 2.0*PI*rng->RandNum();
    do {
        th = PI*rng->RandNum();
        sinTh = sin(th);    // Calculate acceptance probability
        u = rng->RandNum(); // Take random number and decide whether to keep the point, or try again
        if (u <= sinTh) { accept = true; }
        noAttempts += 1;
    } while (accept == false);
    p[2] = radius * sqrt(1.0 - sinTh*sinTh);
    if (th > 0.5*PI) { p[2] *= -1.0; }
    sinTh *= radius;
    double cosPhi = cos(ph);
    p[0] = sinTh * cosPhi;
    p[1] = sinTh * sqrt(1.0 - cosPhi*cosPhi);
    if (ph > PI) { p[1] *= -1.0; }

    return p;
}



//**********************************************************************************************//
// Generate uniformly random points on the surface of an ellipsoid
//**********************************************************************************************//

void SphereMarsaglia::SetParams(double r, shared_ptr<RNG>& rngIn) {
    noAttempts = 0;
    radius = r;
    rng = rngIn;
}

void SphereTrig::SetParams(double r, shared_ptr<RNG>& rngIn) {
    noAttempts = 0;
    radius = r;
    rng = rngIn;
}

Polars EllipsoidSamplerBase::Carts2Pols(Vec3& Carts) {
    Polars p;
    p.ph = atan2(aAxis*Carts[1], bAxis*Carts[0]);
    if (p.ph < 0.0) { p.ph += 2.0*PI; }
    p.th = acos(Carts[2] / cAxis);
    return p;
}

// Naive sphere scaling (gives non-uniform distribution!)
EllipsoidNaive::EllipsoidNaive(double a, double b, double c, shared_ptr<RNG>& rngIn) {
    noAttempts = 0;
    aAxis = a;
    bAxis = b;
    cAxis = c;
    rng = rngIn;
    S.SetParams(1.0, rngIn);
}

// Using the naive spherical scaling method
Vec3 EllipsoidNaive::RandomPoint() {
    Vec3 p = S.RandomPoint();
    p[0] *= aAxis;
    p[1] *= bAxis; // Scale sphere point to ellipsoid (anisotropic scaling)
    p[2] *= cAxis;
    noAttempts += 1;
    return p;
}



// Using the spherical point rejection sampling method
// See Chen and Glotzer 2008
EllipsoidGradRej::EllipsoidGradRej(double a, double b, double c, std::shared_ptr<RNG>& rngIn) {
    noAttempts = 0;
    aAxis = a;
    bAxis = b;
    cAxis = c;
    rng = rngIn;
    S.SetParams(1.0, rngIn);
    // parameters for this algorithm
    gMin = c;
    if (b < gMin) { gMin = b; }
    if (a < gMin) { gMin = a; }
    a4 = 1.0/(a*a*a*a);
    b4 = 1.0/(b*b*b*b);
    c4 = 1.0/(c*c*c*c);

}

Vec3 EllipsoidGradRej::RandomPoint() {
    Vec3 p;
    bool accept = false;
    double g, u;
    do {
        p = S.RandomPoint(); // Generate unit sphere point
        p[0] *= aAxis;
        p[1] *= bAxis; // Scale sphere point to ellipsoid (anisotropic scaling)
        p[2] *= cAxis;
        g = gMin*sqrt(a4*p[0]*p[0] + b4*p[1]*p[1] + c4*p[2]*p[2]); // Calculate acceptance probability
        u = rng->RandNum(); // Take random number and decide whether to keep the point, or try again
        if (u <= g) { accept = true; }
        noAttempts += 1;
    } while (accept == false);
    return p;
}

EllipsoidGradTrig::EllipsoidGradTrig(double a, double b, double c, std::shared_ptr<RNG>& rngIn) {
    noAttempts = 0;
    aAxis = a;
    bAxis = b;
    cAxis = c;
    rng = rngIn;
    S.SetParams(1.0, rngIn);
    // parameters for this algorithm
    gMin = c;
    if (b < gMin) { gMin = b; }
    if (a < gMin) { gMin = a; }
    a4 = 1.0/(a*a*a*a);
    b4 = 1.0/(b*b*b*b);
    c4 = 1.0/(c*c*c*c);

}

Vec3 EllipsoidGradTrig::RandomPoint() {
    Vec3 p;
    bool accept = false;
    double g, u;
    do {
        p = S.RandomPoint(); // Generate unit sphere point
        p[0] *= aAxis;
        p[1] *= bAxis; // Scale sphere point to ellipsoid (anisotropic scaling)
        p[2] *= cAxis;
        g = gMin*sqrt(a4*p[0]*p[0] + b4*p[1]*p[1] + c4*p[2]*p[2]); // Calculate acceptance probability
        u = rng->RandNum(); // Take random number and decide whether to keep the point, or try again
        if (u <= g) { accept = true; }
        noAttempts += 1;
    } while (accept == false);
    return p;
}

// Using area rejection (polar to cartesian)
EllipsoidAreaRej::EllipsoidAreaRej(double a, double b, double c, std::shared_ptr<RNG>& rngIn) {
    noAttempts = 0;
    aAxis = a;
    bAxis = b;
    cAxis = c;
    rng = rngIn;
    // parameters for this algorithm
    bc = b*b*c*c;
    ac = a*a*c*c;
    ab = a*a*b*b;
    double sinThMax = sqrt(-b*b/(2.0*(c*c-b*b)));
    M = sqrt(ac*sinThMax*sinThMax + ab*(1.0 - sinThMax*sinThMax))*sinThMax;
    M = 1.0/M;
}

Vec3 EllipsoidAreaRej::RandomPoint() {
    Vec3 p;
    bool accept = false;
    double m, u;
    double th, ph;
    double sinTh, sinPh, sin2Th, sin2Ph;
    do {
        th = PI*rng->RandNum();
        ph = 2.0*PI*rng->RandNum();
        sinTh = sin(th);
        sin2Th = sinTh * sinTh;
        sinPh = sin(ph);
        sin2Ph = sinPh * sinPh;
        m = M*sqrt(sin2Th*(bc*(1.0 - sin2Ph) + ac*sin2Ph) + ab*(1.0 - sin2Th))*sinTh; // Calculate acceptance probability
        u = rng->RandNum(); // Take random number and decide whether to keep the point, or try again
        if (u <= m) { accept = true; }
        noAttempts += 1;
    } while (accept == false);
    p[0] = aAxis * sinTh * sqrt(1.0 - sin2Ph);
    if (ph > 0.5*PI && ph < 1.5*PI) { p[0] *= -1.0; }
    p[1] = bAxis * sinTh * sinPh;
    p[2] = cAxis * sqrt(1 - sin2Th);
    if (th > 0.5*PI) { p[2] *= -1.0; }
    return p;
}

Polars EllipsoidAreaRej::RandomPointPols() {
    Polars p;
    bool accept = false;
    double m, u;
    double th, ph;
    double sinTh, sinPh, sin2Th, sin2Ph;
    do {
        th = PI*rng->RandNum();
        ph = 2.0*PI*rng->RandNum();
        sinTh = sin(th);
        sin2Th = sinTh * sinTh;
        sinPh = sin(ph);
        sin2Ph = sinPh * sinPh;
        m = M*sqrt(sin2Th*(bc*(1.0 - sin2Ph) + ac*sin2Ph) + ab*(1.0 - sin2Th))*sinTh; // Calculate acceptance probability
        u = rng->RandNum(); // Take random number and decide whether to keep the point, or try again
        if (u <= m) { accept = true; }
        noAttempts += 1;
    } while (accept == false);
    p.th = th;
    p.ph = ph;
    return p;
}


// Using area rejection (polar to cartesian)
SpheroidAreaRej::SpheroidAreaRej(double a, double c, std::shared_ptr<RNG>& rngIn) {
    noAttempts = 0;
    aAxis = a;
    bAxis = a;
    cAxis = c;
    rng = rngIn;
    // parameters for this algorithm
    a4 = a*a*a*a;
    ac = a*a*c*c;
    if (a <= c) { // Prolate (or sphere)
        M = 1.0/(a*c);
    }
    else if (a > c) { // Oblate
        double sinThMax = sqrt(-a*a/(2.0*(c*c-a*a)));
        M = sqrt(ac*sinThMax*sinThMax + a4*(1.0 - sinThMax*sinThMax))*sinThMax;
        M = 1/M;
    }
}

Vec3 SpheroidAreaRej::RandomPoint() {
    Vec3 p;
    bool accept = false;
    double m, u;
    double th;
    double ph = 2.0*PI*rng->RandNum();
    double sinPh = sin(ph);
    double sinTh, sin2Th;
    do {
        th = PI*rng->RandNum();
        sinTh = sin(th);
        sin2Th = sinTh * sinTh;
        m = M*sqrt(sin2Th*ac + a4*(1.0 - sin2Th))*sinTh; // Calculate acceptance probability
        u = rng->RandNum(); // Take random number and decide whether to keep the point, or try again
        if (u <= m) { accept = true; }
        noAttempts += 1;
    } while (accept == false);
    p[0] = aAxis * sinTh * sqrt(1.0 - sinPh*sinPh);
    if (ph > 0.5*PI && ph < 1.5*PI) { p[0] *= -1.0; }
    p[1] = bAxis * sinTh * sinPh;
    p[2] = cAxis * sqrt(1 - sin2Th);
    if (th > 0.5*PI) { p[2] *= -1.0; }
    return p;
}

Polars SpheroidAreaRej::RandomPointPols() {
    Polars p;
    bool accept = false;
    double m, u;
    double th;
    double ph = 2.0*PI*rng->RandNum();
    double sinTh, sin2Th;
    do {
        th = PI*rng->RandNum();
        sinTh = sin(th);
        sin2Th = sinTh * sinTh;
        m = M*sqrt(sin2Th*ac + a4*(1.0 - sin2Th))*sinTh; // Calculate acceptance probability
        u = rng->RandNum(); // Take random number and decide whether to keep the point, or try again
        if (u <= m) { accept = true; }
        noAttempts += 1;
    } while (accept == false);
    p.th = th;
    p.ph = ph;
    return p;
}


// Using area rejection (Mercator coordiantes)
EllipsoidAreaMer::EllipsoidAreaMer(double a, double b, double c, std::shared_ptr<RNG>& rngIn) {
    noAttempts = 0;
    aAxis = a;
    bAxis = b;
    cAxis = c;
    rng = rngIn;
    // parameters for this algorithm
    bc = b*b*c*c;
    ac = a*a*c*c;
    ab = a*a*b*b;
    M = 1.0/(a*c);
}

Vec3 EllipsoidAreaMer::RandomPoint() {
    Vec3 p;
    bool accept = false;
    double m, u, v, ph, sechv, sech2v, sinPh, sin2Ph;
    do {
        v = 2*PI*(-1.0 + 2.0*rng->RandNum());
        ph = 2.0*PI*rng->RandNum();
        sechv = exp(v);
        sechv = 2.0*sechv / (sechv*sechv + 1.0);
        sech2v = sechv * sechv;
        sinPh = sin(ph);
        sin2Ph = sinPh * sinPh;
        m = M*sqrt(sech2v*(bc*(1.0 - sin2Ph) + ac*sin2Ph) + ab*(1.0 - sech2v))*sech2v; // Calculate acceptance probability
        u = rng->RandNum(); // Take random number and decide whether to keep the point, or try again
        if (u <= m) { accept = true; }
        noAttempts += 1;
    } while (accept == false);

    p[0] = aAxis * sechv * sqrt(1.0 - sin2Ph);
    if (ph > 0.5*PI && ph < 1.5*PI) { p[0] *= -1.0; }
    p[1] = bAxis * sechv * sinPh;
    p[2] = cAxis * sqrt(1 - sech2v);
    if (v < 0) { p[2] *= -1.0; }
    return p;
}

// Using area rejection (Mercator coordiantes)
SpheroidAreaMer::SpheroidAreaMer(double a, double c, std::shared_ptr<RNG>& rngIn) {
    noAttempts = 0;
    aAxis = a;
    bAxis = a;
    cAxis = c;
    rng = rngIn;
    // parameters for this algorithm
    a2 = a*a;
    c2 = c*c;
    M = 1.0/(a*c);
}

Vec3 SpheroidAreaMer::RandomPoint() {
    Vec3 p;
    bool accept = false;
    double m, u, v, sechv, sech2v;
    double ph = 2.0*PI*rng->RandNum();
    double sinPh = sin(ph);
    do {
        v = 2*PI*(-1.0 + 2.0*rng->RandNum());
        sechv = exp(v);
        sechv = 2.0*sechv / (sechv*sechv + 1.0);
        sech2v = sechv * sechv;
        m = M*sqrt(a2*(1-sech2v) + c2*sech2v)*aAxis*sech2v; // Calculate acceptance probability
        u = rng->RandNum(); // Take random number and decide whether to keep the point, or try again
        if (u <= m) { accept = true; }
        noAttempts += 1;
    } while (accept == false);
    p[0] = aAxis * sechv * sqrt(1.0 - sinPh*sinPh);
    if (ph > 0.5*PI && ph < 1.5*PI) { p[0] *= -1.0; }
    p[1] = bAxis * sechv * sinPh;
    p[2] = cAxis * sqrt(1 - sech2v);
    if (v < 0) { p[2] *= -1.0; }
    return p;
}

// Generic surface sampler of Detwiler et al 2008
EllipsoidGeneric::EllipsoidGeneric(double a, double b, double c, std::shared_ptr<RNG>& rngIn) {
    noAttempts = 0;
    aAxis = a;
    bAxis = b;
    cAxis = c;
    rng = rngIn;
    S.SetParams(1.0, rngIn);
    // parameters for this algorithm
    R = a;
    if (b > R) { R = b; }
    if (c > R) { R = c; }
    a2 = 1.0/(a*a);
    b2 = 1.0/(b*b);
    c2 = 1.0/(c*c);
}

Vec3 EllipsoidGeneric::RandomPoint() {
    bool accept;
    Vec3 n, u, v, p;
    double t, alpha, beta, gamma, psi, b0, b1, B;
    int j;
    do {
        // Random point on unit sphere i.e. normal vector at point on sphere
        n = S.RandomPoint();
        // Find two orthogonal vectors in the tangent plane
        u[0] = 1.0;
        u[1] = 1.0;
        u[2] = -(n[0] + n[1])/n[2];
        v[0] = n[2] - n[1]*u[2];
        v[1] = n[0]*u[2] - n[2];
        v[2] = n[1] - n[0];
        Normalise(u);
        Normalise(v);
        // Generate point on disk tangent to bounding sphere at point r
        psi = 2.0*PI*rng->RandNum();
        B = R*sqrt(rng->RandNum());
        // B^2 from uniform distribution in (0, R)
        b0 = cos(psi);
        b1 = B * sqrt(1 - b0*b0);
        b0 *= B;
        if (psi > PI) { b1 *= -1; } // If psi > PI, then sin(psi) is negative
        // Convert point in plane to 3D Cartesians
        p = R*n + b0*u + b1*v;
        // Find intersections of line p0 + tp with the ellipsoid
        alpha = a2*n[0]*n[0] + b2*n[1]*n[1] + c2*n[2]*n[2];
        beta = -2.0 * (a2*n[0]*p[0] + b2*n[1]*p[1] + c2*n[2]*p[2]);
        gamma = a2*p[0]*p[0] + b2*p[1]*p[1] + c2*p[2]*p[2] - 1.0;
        t = beta*beta - 4.0*alpha*gamma;
        if (t < 0) { // No intersections - try again
            accept = false;
        }
        else {
            j = 1 + round(rng->RandNum()); // Random integer (1 or 2)
            if (t > 1e-12) { // Two intersections - select one of them at random
                if (j == 1) { t = 0.5 * (-beta + sqrt(t)) / alpha; }
                else { t = 0.5 * (-beta - sqrt(t)) / alpha; }
                accept = true;
            }
            else if (fabs(t) < 1e-12 && j == 1) { // One intersection and accepted
                t = -0.5*beta/alpha;
                accept = true;
            }
        }
        noAttempts += 1;
    } while (accept == false);
    return p - t*n;
}
