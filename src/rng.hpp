// This file defines a set of wrapper classes for the interface between ESS and the TRNG library.
// Within ESS, one of the classes derived from base class, RNG, is selected based on the input file.
// The base class contains TRNG classes for obtaining uniform random numbers, along with virtual functions to generate them.
// Each derived "RNG_*" class contains a long period random number engine from TRNG and uses it to define RandNum().
//
// The RandNor() function uses Doornik's implementation of the Ziggurat method to generate normally distributed random numbers.
// Within ESS, random numbers of the required type are obtained by calling the RandNum() and RandNor() functions.
//
// Author: Callum Marples

#ifndef RNG_DEF
#define RNG_DEF

#include "trng/yarn5s.hpp"
#include "trng/mt19937_64.hpp"
#include "trng/lagfib4xor.hpp"
#include "trng/uniform01_dist.hpp"
#include "trng/normal_dist.hpp"

class RNG;

// Ziggurat constants and functions
#define ZIGNOR_C 128			       /* number of blocks */
#define ZIGNOR_R 3.442619855899	/* start of the right tail */
				   /* (R * phi(R) + Pr(X>=R)) * sqrt(2\pi) */
#define ZIGNOR_V 9.91256303526217e-3
void zigNorInit(int iC, double dR, double dV);
double DRanNormalTail(double dMin, int iNegative, RNG* rng);
double DRanNormalZig(RNG* rng);

class RNG {
    protected:
        trng::uniform01_dist<> uniform;                     // Uniform random numbers between 0 and 1
    public:
        RNG() { zigNorInit(ZIGNOR_C, ZIGNOR_R, ZIGNOR_V); }
        virtual double RandNum() = 0;                       // For uniform deviates
        double RandNor() { return DRanNormalZig(this); }    // For normal deviates
};

// YARN3
class RNG_yn : public RNG {
    private:
        trng::yarn5s R;
    public:
        RNG_yn(unsigned long& s) {  R.seed(s); zigNorInit(ZIGNOR_C, ZIGNOR_R, ZIGNOR_V); }
        double RandNum() { return uniform(R); }
};

// Mersenne twister
class RNG_mt : public RNG {
    private:
        trng::mt19937_64 R;
    public:
        RNG_mt(unsigned long& s) {  R.seed(s); zigNorInit(ZIGNOR_C, ZIGNOR_R, ZIGNOR_V); }
        double RandNum() { return uniform(R); }
};

// Additive lagged Fibonacci
class RNG_lf : public RNG {
    private:
        trng::lagfib4xor_19937_64 R;
    public:
        RNG_lf(unsigned long& s) {  R.seed(s); zigNorInit(ZIGNOR_C, ZIGNOR_R, ZIGNOR_V); }
        double RandNum() { return uniform(R); }
};


#endif // RNG_DEF
