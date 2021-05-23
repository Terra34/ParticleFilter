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
