/*
 * Copyright:
 * ----------------------------------------------------------------
 * This confidential and proprietary software may be used only as
 * authorised by a licensing agreement from ARM Limited
 *   (C) COPYRIGHT 2014, 2016 ARM Limited
 *       ALL RIGHTS RESERVED
 * The entire notice above must be reproduced on all authorised
 * copies and copies may only be made to the extent permitted
 * by a licensing agreement from ARM Limited.
 * ----------------------------------------------------------------
 * File:     main.c
 * Release Information : Cortex-M1 DesignStart-r0p1-00rel0
 * ----------------------------------------------------------------
 *
 */

/*
 * --------Included Headers--------
 */

#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <float.h>
#include <stdbool.h>

// Xilinx specific headers
#include "xparameters.h"
#include "m1_for_arty.h" // Project specific header
#include "xil_printf.h"

/*******************************************************************/

#define STCTRL      (*( (volatile uint32_t *) 0xE000E010 ))
#define STRELOAD    (*( (volatile uint32_t *) 0xE000E014 ))
#define STCURR      (*( (volatile uint32_t *) 0xE000E018 ))  
#define PI 3.14159265359
#define TWO_PI 6.28318530718

#define SYSTEM_CLOCK	(100000000UL) /*	HCLK frequency	*/


void wait(int num)
{
	for (int i=0; i < num; i++)
		i++;
}

void newLine(){
	print("\r\n");
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

struct genHelper {
    int *taps;
	int current_tap;
	int _tapCnt;
    int N;
    int M;
    int twoM;
    int threeM;
    int rightShift;
	float mean;
	float std;
    long long current;
    long long shifted;
    long long mask;
};

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

void time_generator()
{
    int iterations = 1000;
	char myString[256];
    int taps[3] = {0, 3, 31};
    struct genHelper _generator;
    initializeGenerator(&_generator, 32, 8, taps, 127.63402f, 73.727015f);
    seedGenerator(&_generator, ((long long)1<<32) - 1);
    float data;
    uint32_t start, end;
	
	STCTRL = (1<<0) | (1<<2);
	wait(1000);
    start = STCURR;

    for (int i=0; i<iterations; i++)
    {
        data = generate(&_generator);
    }
	newLine();
	
    end = STCURR;
	STCTRL = (0<<0);

    sprintf(myString, "LFSR generator clocks elapsed: %fms", (float)(start-end)/SYSTEM_CLOCK*1000);
	print(myString);
	newLine();
}

struct gaussGenHelper 
{
    float gauss;
	bool has_gauss;
};

float gauss(struct gaussGenHelper *state)
{
    if (state->has_gauss){
        const float temp = state->gauss;
        state->has_gauss = false;
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
        state->has_gauss = true;
        return f * x2;
    }
}

void time_gauss()
{
	char myString[256];
    int iterations = 1000;
    struct gaussGenHelper _normState;
    float data;
    uint32_t start, stop;

    _normState.has_gauss = false;
    _normState.gauss = 0.0;

    srand(100);

	
	STCTRL = (1<<0) | (1<<2);
	wait(1000);
    start = STCURR;

    for (int i=0; i<iterations; i++)
    {
        data = gauss(&_normState);
    }

    stop = STCURR;
	STCTRL = (0<<0);

    sprintf(myString, "Gauss generator clocks elapsed: %fms", (float)(start-stop)/SYSTEM_CLOCK*1000);
	print(myString);
	newLine();
}

int main(void){
	SystemInit();
	STRELOAD = 16000000;
	
	time_gauss();
	time_generator();
}
