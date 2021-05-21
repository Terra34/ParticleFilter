#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "generators.h"

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

void lfsr(struct genHelper *helper)
{
    long long A_next = 0, A_upper;
    helper->current = helper->current & helper->mask;
    for (int i=0; i < helper->N - 1; i++)
    {
        A_upper = checkBit(helper->current, helper->N - 1);
        //if (arrayContains(helper->taps, helper->tapsLen, i))
        //    A_next += (checkBit(helper->current, i) ^ A_upper) << ((i + 1) % helper->N);
		if (i == helper->current_tap)
		{
			A_next += (checkBit(helper->current, i) ^ A_upper) << ((i + 1) % helper->N);
			helper->current_tap = helper->taps[++helper->_tapCnt];
		}
        else
            A_next += checkBit(helper->current, i) << ((i + 1) % helper->N);
		// if helper->N is a power of 2 then a faster way would be ((i + 1) & (helper->N - 1))
    }
    A_next += A_upper;

    helper->current = A_next;
	helper->current_tap = helper->taps[0];
	helper->_tapCnt = 0;
}

float generate(struct genHelper *helper)
{
    lfsr(helper);
    return (((( helper->current & helper->shifted ) + 
            ( (helper->current & (helper->shifted<<helper->M))>>helper->M ) +
            ( (helper->current & (helper->shifted<<helper->twoM))>>helper->twoM ) +
            ( (helper->current & (helper->shifted<<helper->threeM))>>helper->threeM )
            ) & helper->shifted) - helper->mean) / helper->std;
}

void initializeGenerator(struct genHelper *helper, int N, int M, int *taps, float mean, float std)
{
	// taps needs to be a sorted array (ascending) for example: [0 1 5 31]
    // maybe include sorting in initialization
    helper->N = N;
    helper->M = M;
    helper->taps = taps;
	helper->current_tap = taps[0];
	helper->_tapCnt = 0;
    helper->twoM = 2*M;
    helper->threeM = 3*M;
    helper->mask = ((long long)1<<N) - 1;
    helper->shifted = ((long long)1<<M) - 1;
	helper->mean = mean;
	helper->std = std;
}

void seedGenerator(struct genHelper *helper, long long seed)
{
    helper->current = seed;
}

void initializeGauss(struct gaussGenHelper *state)
{
    state->has_gauss = 0;
    state->gauss = 0.0;
}

float gauss(struct gaussGenHelper *state)
{
    if (state->has_gauss){
        const float temp = state->gauss;
        state->has_gauss = 0;
        state->gauss = 0.0;
        return temp;
    }
    else {
        float f, x1, x2, r2;

        do {
            x1 = ((float)rand() / RAND_MAX) * 2 - 1;
            x2 = ((float)rand() / RAND_MAX) * 2 - 1;
            r2 = x1 * x1 + x2 * x2;
        } while (r2 >= 1.0 || r2 == 0.0);

        f = sqrt(-2.0 * log(r2) / r2);

        state->gauss = f * x1;
        state->has_gauss = 1;
        return f * x2;
    }
}