#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "generators.h"

static unsigned int jz, jsr = 123456789;

#define SHR3 (jz=jsr, jsr^=(jsr<<13), jsr^=(jsr>>17), jsr^=(jsr<<5),jz+jsr)
#define UNI (.5 + (signed) SHR3 * .2328306e-9)

void setSeed(unsigned int seed)
{
    jsr = seed;
}

float uniform(void)
{
    // Implementation of Marsaglia's SHR3 algorithm for generating random uniformly distributed numbers
    return UNI;
}

long long checkBit(long long num, int bit)
{
    // Checks whether or not a bit at index bit is set in num
    return (num & (1 << bit)) > 0;
}

int arrayContains(int *array, int length, int number)
{
    /*
        array => array to search
        length => length of the array to search
        number => number to search for
        returns => 1 if found, else 0
    */
    for (int i=0; i<length; i++)
    {
        if(array[i] == number)
            return 1;
    }
    return 0;
}

void lfsr(struct genState *state)
{
    long long A_next = 0, A_upper, A;
    A = state->current & state->mask;
    A_upper = checkBit(A, state->N - 1);
	for (int i=0; i < state->N - 1; i++)
    {
        //if (arrayContains(state->taps, state->tapsLen, i))
        //    A_next += (checkBit(state->current, i) ^ A_upper) << ((i + 1) % state->N);
		if (i == state->current_tap)
		{
			A_next += (checkBit(A, i) ^ A_upper) << ((i + 1) % state->N);
			state->current_tap = state->taps[++state->_tapCnt];
		}
        else
            A_next += checkBit(A, i) << ((i + 1) % state->N);
		// if state->N is a power of 2 then a faster way would be ((i + 1) & (state->N - 1))
    }
    A_next += A_upper;

    state->current = A_next;
	state->current_tap = state->taps[0];
	state->_tapCnt = 0;
}

float generate(struct genState *state)
{
    // Implementation of linear feedback shift register random number generator
    lfsr(state);
	return 	((state->current & state->shifted) +
			((state->current & (state->shifted<<state->M))>>state->M) +
			((state->current & (state->shifted<<state->twoM))>>state->twoM) +
			((state->current & (state->shifted<<state->threeM))>>state->threeM)
				- state->mean) / state->std;
}

void initializeGenerator(struct genState *state, int N, int M, int *taps, float mean, float std)
{
	// taps needs to be a sorted array (ascending) for example: [0 1 5 31]
    // maybe include sorting in initialization
    state->N = N;
    state->M = M;
    state->taps = taps;
	state->current_tap = taps[0];
	state->_tapCnt = 0;
    state->twoM = 2*M;
    state->threeM = 3*M;
    state->mask = ((long long)1<<N) - 1;
    state->shifted = ((long long)1<<M) - 1;
	state->mean = mean;
	state->std = std;
}

void seedGenerator(struct genState *state, long long seed)
{
    state->current = seed;
}

void initializeGauss(struct gaussGenState *state)
{
    state->has_gauss = 0;
    state->gauss = 0.0;
}

float gauss(struct gaussGenState *state)
{
    /*
        Implementation of Marsaglia's polar method for calculating normally distributed gaussian variables
        seeding of rand needs to be done outside of function (with setSeed())
    */
    if (state->has_gauss){
        const float temp = state->gauss;
        state->has_gauss = 0;
        state->gauss = 0.0;
        return temp;
    }   
    else {
        float f, x1, x2, r2;

        do {
            x1 = UNI * 2 - 1;
            x2 = UNI * 2 - 1;
            r2 = x1 * x1 + x2 * x2;
        } while (r2 >= 1.0);

        f = sqrt((-2.0/r2) * log(r2));

        state->gauss = f * x1;
        state->has_gauss = 1;
        return f * x2;
    }
}

float gaussbm(struct gaussGenState *state)
{
    /*
        Implementation of Box-Muller method for calculating normally distributed gaussian variables
        seeding of rand needs to be done outside of function (with setSeed())
    */
    if (state->has_gauss){
        const float temp = state->gauss;
        state->has_gauss = 0;
        state->gauss = 0.0;
        return temp;
    }
    else {
        float u, v, f1, f2;

        u = UNI;
        v = UNI;

        f1 = sqrt(-2.0 * log(u));
        f2 = 2*M_PI*v;

        state->gauss = f1 * cos(f2);
        state->has_gauss = 1;
        return f1 * sin(f2);
    }
}

float gaussInv(void)
{
    /*
        Implementation of Inverse cumulative distribution method for calculating normally distributed gaussian variables
        Approximation relative error less than 1.15 x 10e-9 in the entire region.
        seeding of rand needs to be done outside of function (with setSeed())
    */
    float p, q, r;

    float a[6] = {-3.969683028665376e+01,  2.209460984245205e+02,
                    -2.759285104469687e+02,  1.383577518672690e+02,
                    -3.066479806614716e+01,  2.506628277459239e+00};
    float b[5] = {-5.447609879822406e+01,  1.615858368580409e+02,
                    -1.556989798598866e+02,  6.680131188771972e+01,
                    -1.328068155288572e+01};
    float c[6] = {-7.784894002430293e-03, -3.223964580411365e-01,
                    -2.400758277161838e+00, -2.549732539343734e+00,
                    4.374664141464968e+00,  2.938163982698783e+00};
    float d[4] = { 7.784695709041462e-03,  3.224671290700398e-01,
                    2.445134137142996e+00,  3.754408661907416e+00};

    p = UNI;

    if (p < 0.02425){
        q = sqrt(-2*log(p));
        return ((((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
               ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1));
    }

    if ((1-0.02425) < p){
        q = sqrt(-2*log(1-p));
        return -((((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / 
                ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1));
    }

    q = p - 0.5;
    r = q*q;
    return ((((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q /
           (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1));
}

// ------------------- ZIGGURAT -------------------

static unsigned int iz, kn[128];
static int hz;
static float wn[128], fn[128];

#define RNOR (hz=SHR3, iz=hz&127, (abs(hz)<kn[iz])? hz * wn[iz] : nfix())

float nfix(void)
{
    const float r = 3.442619855899;
    static float x, y;
    for (;;)
    {
        x = hz * wn[iz];
        if (iz == 0)
        {
            do { x = -log(UNI) * 0.2904764516; y = -log(UNI); }
            while (y + y < x * x);
            return (hz > 0) ? r + x : -r - x;
        }
        if (fn[iz] + UNI * (fn[iz - 1] - fn[iz]) < exp(-.5 * x * x)) return x;

        hz = SHR3;
        iz = hz & 127;
        if (abs(hz) < kn[iz]) return (hz * wn[iz]);
    }

}

void zigset(void)
{
    // Setup all necessary tables for generating random normally distributed variables
    const float m1 = 2147483648.0;
    float dn = 3.442619855899, tn = dn, vn = 9.91256303526217e-3, q;
    int i;

    q = vn / exp(-.5 * dn * dn);
    kn[0] = (dn / q) * m1;
    kn[1] = 0;

    wn[0] = q / m1;
    wn[127] = dn / m1;

    fn[0] = 1.;
    fn[127] = exp(-.5 * dn * dn);

    for (i = 126; i >= 1; i--)
    {
        dn = sqrt(-2. * log(vn / dn + exp(-.5 * dn * dn)));
        kn[i + 1] = (dn / tn) * m1;
        tn = dn;
        fn[i] = exp(-.5 * dn * dn);
        wn[i] = dn / m1;
    }
}

float ziggurat(void)
{
    /*
        Implementation of Marsaglia's Ziggurat algorithm for generating random normally distributed numbers
        Before using, setSeed() and zigset() need to be called, otherwise all 0's returned
    */
    return RNOR;
}

// ------------------- END ZIGGURAT -------------------