/* 3-element vector class, "Vec3"
   Author: Callum Marples
*/

#ifndef VECTORDEF
#define VECTORDEF

#include <cassert>
#include <array>
#include <cmath>

//******************************************************************************
// This class implements 3-element vectors
//******************************************************************************
class Vec3 {
    private:
        std::array<double, 3> mData;
    public:
        // Constructors
        Vec3(double x, double y, double z) : mData{x, y, z} {}
        Vec3() : mData{0.0, 0.0, 0.0} {}
        Vec3(const Vec3& otherVec3) : mData{otherVec3.getVector()} {}
        // Set components
        void setVector(const double& x, const double& y, const double& z) { mData[0] = x, mData[1] = y, mData[2] = z; }
		// Get components
		double& operator[](int i) { assert(i > -1 && i < 3); return mData[i]; }
		std::array<double, 3> getVector() const { return std::array<double, 3> {mData[0], mData[1], mData[2]}; }
        // Arithmetic operators
        Vec3& operator=(const Vec3& otherVec3) { mData = otherVec3.getVector(); return *this; }
        Vec3& operator+=(const Vec3& u) { for (int i=0; i!=3; ++i) { mData[i] += u.mData[i]; } return *this; }
        Vec3& operator-=(const Vec3& u) { for (int i=0; i!=3; ++i) { mData[i] -= u.mData[i]; } return *this; }
        Vec3& operator*=(const double& x) { for (auto& v : mData) { v *= x; } return *this; }
        Vec3& operator/=(const double& x) { assert(x != 0.0); for (auto& v : mData) { v /= x; } return *this; }
};

// Arithmetic operators
inline Vec3 operator+(Vec3 u, Vec3 v) { return u += v; }
inline Vec3 operator-(Vec3 u, Vec3 v) { return u -= v; }
inline Vec3 operator-(Vec3 u) { return {-u[0], -u[1], -u[2]}; }
inline Vec3 operator*(Vec3 v, double x) { return v *= x; }
inline Vec3 operator*(double x, Vec3 v) { return v *= x; }
inline Vec3 operator/(Vec3 v, double x) { assert(x != 0.0); return v /= x; }

// Products
inline double Dot(Vec3& u, Vec3& v) { return (u[0]*v[0] + u[1]*v[1] + u[2]*v[2]); }
inline Vec3 Hadamard(Vec3& u, Vec3& v) { return {u[0]*v[0], u[1]*v[1], u[2]*v[2]}; }
inline Vec3 Cross(Vec3& u, Vec3& v) { return {u[1]*v[2] - u[2]*v[1], u[2]*v[0] - u[0]*v[2], u[0]*v[1] - u[1]*v[0]}; }
inline double NormSquared(Vec3& u) { return (u[0]*u[0] + u[1]*u[1] + u[2]*u[2]); }
inline double Norm(Vec3& u) { return (sqrt(NormSquared(u))); }
inline void Normalise(Vec3& u) { double norm = Norm(u); assert(norm != 0.0); u /= norm; return; }
inline Vec3 Reciprocals(Vec3& u) { assert(u[0] != 0.0 && u[1] != 0.0 && u[2] != 0.0); return {1.0/u[0], 1.0/u[1], 1.0/u[2]}; }


#endif // VECTORDEF
