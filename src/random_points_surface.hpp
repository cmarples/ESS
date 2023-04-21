/* Generation of random points on the surface of a sphere/ellipsoid
   Author: Callum Marples
*/

#ifndef SAMPLINGDEF
#define SAMPLINGDEF

#include "rng.hpp"
#include "vectors.hpp"
#include <memory>

// Sphere surface samplers. These algorithms all generate points on the surface
// of a unit sphere.
class SphereSamplerBase {
    protected:
        double radius;
        int noAttempts;
        std::shared_ptr<RNG> rng;
    public:
        int GetNoAttempts() { return noAttempts; }
};

class SphereRejection : public SphereSamplerBase {
    public:
        SphereRejection(double r, std::shared_ptr<RNG>& rngIn) { noAttempts = 0; radius = r; rng = rngIn; }
        Vec3 RandomPoint();
};

class SphereMarsaglia : public SphereSamplerBase {
    public:
        SphereMarsaglia(double r, std::shared_ptr<RNG>& rngIn) { noAttempts = 0; radius = r; rng = rngIn; }
        SphereMarsaglia() { noAttempts = 0; radius = 1.0; }
        void SetParams(double r, std::shared_ptr<RNG>& rngIn);
        Vec3 RandomPoint();
};

class SphereCook : public SphereSamplerBase {
    public:
        SphereCook(double r, std::shared_ptr<RNG>& rngIn) { noAttempts = 0; radius = r; rng = rngIn; }
        Vec3 RandomPoint();
};

class SphereGaussian : public SphereSamplerBase {
    public:
        SphereGaussian(double r, std::shared_ptr<RNG>& rngIn) { noAttempts = 0; radius = r; rng = rngIn; }
        Vec3 RandomPoint();
};

class SphereTrig : public SphereSamplerBase {
    public:
        SphereTrig(double r, std::shared_ptr<RNG>& rngIn) { noAttempts = 0; radius = r; rng = rngIn; }
        SphereTrig() { noAttempts = 0; radius = 1.0; }
        void SetParams(double r, std::shared_ptr<RNG>& rngIn);
        Vec3 RandomPoint();
};

class SphereAreaRejection : public SphereSamplerBase {
    public:
        SphereAreaRejection(double r, std::shared_ptr<RNG>& rngIn) { noAttempts = 0; radius = r; rng = rngIn; }
        Vec3 RandomPoint();
};




struct Polars {
    double th;
    double ph;
};


// Ellipsoid samplers. These algorithms generate random points on the surface of an ellipsoid
// with axis lengths a, b and c.
class EllipsoidSamplerBase {
    protected:
        double aAxis;
        double bAxis;
        double cAxis;
        int noAttempts;
        std::shared_ptr<RNG> rng;
    public:
        int GetNoAttempts() { return noAttempts; }
        Polars Carts2Pols(Vec3& Carts);
};

// This method generates non-uniform points.
class EllipsoidNaive : public EllipsoidSamplerBase {
    public:
        EllipsoidNaive(double a, double b, double c, std::shared_ptr<RNG>& rngIn);
        SphereMarsaglia S;
        Vec3 RandomPoint();
};

class EllipsoidGradRej : public EllipsoidSamplerBase {
    private:
        double gMin;
        double a4;
        double b4;
        double c4;
    public:
        EllipsoidGradRej(double a, double b, double c, std::shared_ptr<RNG>& rngIn);
        SphereMarsaglia S;
        Vec3 RandomPoint();
};

class EllipsoidGradTrig : public EllipsoidSamplerBase {
    private:
        double gMin;
        double a4;
        double b4;
        double c4;
    public:
        EllipsoidGradTrig(double a, double b, double c, std::shared_ptr<RNG>& rngIn);
        SphereTrig S;
        Vec3 RandomPoint();
};

class EllipsoidAreaRej : public EllipsoidSamplerBase {
    private:
        double bc;
        double ac;
        double ab;
        double M;
    public:
        EllipsoidAreaRej(double a, double b, double c, std::shared_ptr<RNG>& rngIn);
        Vec3 RandomPoint();
        Polars RandomPointPols();
};

class SpheroidAreaRej : public EllipsoidSamplerBase {
    private:
        double ac;
        double a4;
        double M;
    public:
        SpheroidAreaRej(double a, double c, std::shared_ptr<RNG>& rngIn);
        Vec3 RandomPoint();
        Polars RandomPointPols();
};

class EllipsoidAreaMer : public EllipsoidSamplerBase {
    private:
        double bc;
        double ac;
        double ab;
        double M;
    public:
        EllipsoidAreaMer(double a, double b, double c, std::shared_ptr<RNG>& rngIn);
        Vec3 RandomPoint();
};

class SpheroidAreaMer : public EllipsoidSamplerBase {
    private:
        double a2;
        double c2;
        double M;
    public:
        SpheroidAreaMer(double a, double c, std::shared_ptr<RNG>& rngIn);
        Vec3 RandomPoint();
};

class EllipsoidGeneric : public EllipsoidSamplerBase {
    private:
        double R;
        double a2;
        double b2;
        double c2;
    public:
        EllipsoidGeneric(double a, double b, double c, std::shared_ptr<RNG>& rngIn);
        SphereMarsaglia S;
        Vec3 RandomPoint();

};


#endif  // SAMPLINGDEF
