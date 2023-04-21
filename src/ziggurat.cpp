// Ziggurat Method for generating normally distributed random numbers

#include "rng.hpp"
#include <cmath>

// The code below is a modified version of the Ziggurat method code of Doornik.
// The random number calls use the generators defined in rng.hpp.

/*==========================================================================
 *  This code is Copyright (C) 2005, Jurgen A. Doornik.
 *  Permission to use this code for non-commercial purposes
 *  is hereby given, provided proper reference is made to:
 *		Doornik, J.A. (2005), "An Improved Ziggurat Method to Generate Normal
 *          Random Samples", mimeo, Nuffield College, University of Oxford,
 *			and www.doornik.com/research.
 *		or the published version when available.
 *	This reference is still required when using modified versions of the code.
 *  This notice should be maintained in modified versions of the code.
 *	No warranty is given regarding the correctness of this code.
 *==========================================================================*/

 /*------------------------------ General Ziggurat --------------------------*/
double DRanNormalTail(double dMin, int iNegative, RNG* rng)
{
	double x, y;
	do
	{	x = log(rng->RandNum()) / dMin;
		y = log(rng->RandNum());
	} while (-2 * y < x * x);
	return iNegative ? x - dMin : dMin - x;
}

// Constants already defined in rng.hpp
//#define ZIGNOR_C 128			       /* number of blocks */
//#define ZIGNOR_R 3.442619855899	/* start of the right tail */
//				   /* (R * phi(R) + Pr(X>=R)) * sqrt(2\pi) */
//#define ZIGNOR_V 9.91256303526217e-3

/* s_adZigX holds coordinates, such that each rectangle has*/
/* same area; s_adZigR holds s_adZigX[i + 1] / s_adZigX[i] */
static double s_adZigX[ZIGNOR_C + 1], s_adZigR[ZIGNOR_C];

void zigNorInit(int iC, double dR, double dV)
{
	int i;	double f;

	f = exp(-0.5 * dR * dR);
	s_adZigX[0] = dV / f; /* [0] is bottom block: V / f(R) */
	s_adZigX[1] = dR;
	s_adZigX[iC] = 0;

	for (i = 2; i < iC; ++i)
	{
		s_adZigX[i] = sqrt(-2 * log(dV / s_adZigX[i - 1] + f));
		f = exp(-0.5 * s_adZigX[i] * s_adZigX[i]);
	}
	for (i = 0; i < iC; ++i)
		s_adZigR[i] = s_adZigX[i + 1] / s_adZigX[i];
}
double DRanNormalZig(RNG* rng)
{
	unsigned int i;
	double x, u, f0, f1;

	for (;;)
	{
		u = 2 * rng->RandNum() - 1;
		i = std::floor(rng->RandNum() * 128);
		/* first try the rectangular boxes */
		if (fabs(u) < s_adZigR[i])
			return u * s_adZigX[i];
		/* bottom box: sample from the tail */
		if (i == 0)
			return DRanNormalTail(ZIGNOR_R, u < 0, rng);
		/* is this a sample from the wedges? */
		x = u * s_adZigX[i];
		f0 = exp(-0.5 * (s_adZigX[i] * s_adZigX[i] - x * x) );
		f1 = exp(-0.5 * (s_adZigX[i+1] * s_adZigX[i+1] - x * x) );
      	if (f1 + rng->RandNum() * (f0 - f1) < 1.0)
			return x;
	}
}

