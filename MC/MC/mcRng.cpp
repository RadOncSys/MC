#include "mcRng.h"
#include <exception>

mcRng::mcRng(void) : p96Gen_(nullptr), pIGen_(nullptr), pJGen_(nullptr), decrement_(0)
{
	//init(33, 97);
}

mcRng::mcRng(int ijSeed, int klSeed)
{
	init(ijSeed, klSeed);
}

mcRng::~mcRng(void)
{
}

// Initializes the RANMAR random number generator (see below). The first
// seed, ijSeed, is used to initialize the variables i and j. The second
// seed, klSeed, is used to initialize the variables k and l. The seeds
// must obey the following conditions:
//
//         0 <= ijSeed <= 31328   (where 31328 = 177 * 177 - 1)
//         0 <= klSeed <= 30081   (where 30081 = 169 * 178 - 1)
//
// A small change from the EGS4 version of RANMAR is that (0, 0) is here
// accepted as a valid pair of seeds. In EGS4, (0, 0) would be changed to
// a default of (1802, 9373), which has significance only as a standard test.
// If the RANMAR algorithm has been coded correctly and those two seeds are
// chosen, then the initial values of (i, j, k, l) are (12, 34, 56, 78), and 
// multiplying numbers 20001 through 20006 of the random number sequence by
// 16777216.0 (=2^24) yields the integers
//
//                 6533892.0,  14220222.0,   7275067.0,
//                 6172232.0,   8354498.0,  10633180.0.
//
// Inasmuch as there are over 3 x 10^4 possible values of each of the two
// seeds, RANMAR can generate over 9 x 10^8 different random number sequences,
// each with a length of about 10^30. The existence of two seeds facilitates
// cooperative efforts: if each cooperating party is assigned a first seed
// and varies only the second, then statistical independence of all the runs
// is assured.

void mcRng::init(int ijSeed, int klSeed)
{
	// Check the validity of the two seeds:
	if (ijSeed < 0)
		throw std::exception("mcRng::init: negative ijSeed");

	if (ijSeed > 31328)
		throw std::exception("mcRng::init: ijSeed exceeds 31328");

	if (klSeed < 0)
		throw std::exception("mcRng::init: negative klSeed");

	if (klSeed > 30081)
		throw std::exception("mcRng::init: klSeed exceeds 30081");

	decrement_ = 362436.0 / 16777216.0;

	// Use the seeds to compute the initial values of i, j, k, l:
	int i = ((ijSeed / 177) % 177) + 2;
	int j = (ijSeed % 177) + 2;
	int k = ((klSeed / 169) % 178) + 1;
	int l = klSeed % 169;

	// Compute 97 generator values to 24 binary digits (bits) each:
	for (int iGen = 0; iGen < 97; iGen++)
	{
		genArray_[iGen] = 0.0;
		double powerOfTwo = 1.0;
		for (int iBit = 0; iBit < 24; iBit++)
		{
			// Set the decimal number represented by the bit:
			powerOfTwo *= 0.5;
			// By adding or not adding the decimal number, set the bit to one or zero:
			int m = (((i * j) % 179) * k) % 179;
			i = j;
			j = k;
			k = m;
			l = (53 * l + 1) % 169;
			if ((l * m) % 64 >= 32)
				genArray_[iGen] += powerOfTwo;
		}
	}

	/* Set the pointers to their initial values: */
	pIGen_ = p96Gen_ = genArray_ + 96;
	pJGen_ = genArray_ + 32;
}


// Implementation of the RANMAR random number generator. 
// References:
//     F. James, Comp. Phys. Comm. 60 (1990), 329-44.
//     F. James, "A review of pseudo-random number generators,"
//         CERN Report SOFTWR 88-20 (1988).
//     G. Marsaglia and A. Zaman, "Toward a universal random number
//         generator," Florida State University Report FSU-SCRI-87-50 (1987).

double mcRng::rnd()
{
	const double decrementOfDecrement = 7654321.0 / 16777216.0;
	const double incrementOfDecrement = 16777213.0 / 16777216.0;

	// Compute a new generator:
	double gen = *pIGen_ - *pJGen_;
	if (gen < 0.0) gen += 1.0;

	// Recompute the decrement:
	decrement_ -= decrementOfDecrement;
	if (decrement_ < 0.0) decrement_ += incrementOfDecrement;

	// Decrement the generator to obtain the random number:
	double rnd = gen - decrement_;
	if (rnd < 0.0) rnd += 1.0;

	// Save the generator and move the pointers for use during
	// the next function call:
	*pIGen_ = gen;
	if (--pIGen_ < genArray_) pIGen_ = p96Gen_;
	if (--pJGen_ < genArray_) pJGen_ = p96Gen_;

	return rnd;
}
