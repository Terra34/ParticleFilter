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

// Xilinx specific headers
#include "xparameters.h"
#include "m1_for_arty.h" // Project specific header
#include "xil_printf.h"
#include "generators.h"
#include "helpers.h"
#include "particleFilter.h"

/*******************************************************************/

#define STCTRL      (*( (volatile uint32_t *) 0xE000E010 ))
#define STRELOAD    (*( (volatile uint32_t *) 0xE000E014 ))
#define STCURR      (*( (volatile uint32_t *) 0xE000E018 ))  

#define SYSTEM_CLOCK	(100000000UL) /*	HCLK frequency	*/

void time_lfsr(int iterations)
{
	char myString[256];
    int taps[3] = {0, 3, 31};
    struct genState state;
    initializeGenerator(&state, 32, 8, taps, 517.21916f, 150.88634811f);
    seedGenerator(&state, ((long long)1<<32) - 1);
    float data;
    uint32_t start, end;
	
	STCTRL = (1<<0) | (1<<2);
	wait(10);
    start = STCURR;

    for (int i=0; i<iterations; i++)
    {
        data = generate(&state);
    }
	newLine();
	
    end = STCURR;
	STCTRL = (0<<0);

    sprintf(myString, "LFSR generator time elapsed: %fms", (float)(start-end)/SYSTEM_CLOCK*1000);
	print(myString);
	newLine();
}

void test_lfsr(int iterations, int seed)
{
	char myString[32];
    float data;
	int taps[3] = {0, 3, 31};
    struct genState state;
    initializeGenerator(&state, 32, 8, taps, 517.21916f, 150.88634811f);
    seedGenerator(&state, ((long long)1<<32) - 1);
    

    for (int i=0; i<iterations; i++)
    {
        data = generate(&state);
		sprintf(myString, "%f", data);
		print(myString);
		newLine();
    }	
}

void time_uniform(int iterations)
{
	char myString[256];
	setSeed(100);
    float data;
    uint32_t start, end;
	
	STCTRL = (1<<0) | (1<<2);
	wait(10);
    start = STCURR;

    for (int i=0; i<iterations; i++)
    {
        data = uniform();
    }
	newLine();
	
    end = STCURR;
	STCTRL = (0<<0);

    sprintf(myString, "Uniform(SHR3) generator time elapsed: %fms", (float)(start-end)/SYSTEM_CLOCK*1000);
	print(myString);
	newLine();
}

void test_uniform(int iterations, int seed)
{
	char myString[32];
	setSeed(seed);
    float data;

    for (int i=0; i<iterations; i++)
    {
        data = uniform();
		sprintf(myString, "%f", data);
		print(myString);
		newLine();
    }	
}

void time_boxmuller(int iterations)
{
	char myString[256];
	float data;
    uint32_t start, end;
	struct gaussGenState state;
	setSeed(100);
	initializeGauss(&state);
	
	STCTRL = (1<<0) | (1<<2);
	wait(10);
    start = STCURR;

    for (int i=0; i<iterations; i++)
    {
        data = gaussbm(&state);
    }
	newLine();
	
    end = STCURR;
	STCTRL = (0<<0);

    sprintf(myString, "Box-Muller gauss generator time elapsed: %fms", (float)(start-end)/SYSTEM_CLOCK*1000);
	print(myString);
	newLine();
}

void test_boxmuller(int iterations, int seed)
{
	char myString[32];
	setSeed(seed);
    float data;
	struct gaussGenState state;
	initializeGauss(&state);

    for (int i=0; i<iterations; i++)
    {
        data = gaussbm(&state);
		sprintf(myString, "%f", data);
		print(myString);
		newLine();
    }	
}

void time_marsaglia(int iterations)
{
    char myString[256];
	float data;
    uint32_t start, end;
	struct gaussGenState state;
	setSeed(100);
	initializeGauss(&state);
	
	STCTRL = (1<<0) | (1<<2);
	wait(10);
    start = STCURR;

    for (int i=0; i<iterations; i++)
    {
        data = gauss(&state);
    }
	newLine();
	
    end = STCURR;
	STCTRL = (0<<0);

    sprintf(myString, "Marsaglia polar gauss generator time elapsed: %fms", (float)(start-end)/SYSTEM_CLOCK*1000);
	print(myString);
	newLine();
}

void test_marsaglia(int iterations, int seed)
{
	char myString[32];
	setSeed(seed);
    float data;
	struct gaussGenState state;
	initializeGauss(&state);

    for (int i=0; i<iterations; i++)
    {
        data = gauss(&state);
		sprintf(myString, "%f", data);
		print(myString);
		newLine();
    }	
}

void time_inverse(int iterations)
{
	char myString[256];
	float data;
    uint32_t start, end;
	setSeed(100);
	
	STCTRL = (1<<0) | (1<<2);
	wait(10);
    start = STCURR;

    for (int i=0; i<iterations; i++)
    {
        data = gaussInv();
    }
	newLine();
	
    end = STCURR;
	STCTRL = (0<<0);

    sprintf(myString, "Inverse cumulative distribution approximation gauss generator time elapsed: %fms", (float)(start-end)/SYSTEM_CLOCK*1000);
	print(myString);
	newLine();
}

void test_inverse(int iterations, int seed)
{
	char myString[32];
	setSeed(seed);
    float data;

    for (int i=0; i<iterations; i++)
    {
        data = gaussInv();
		sprintf(myString, "%f", data);
		print(myString);
		newLine();
    }	
}

void time_ziggurat(int iterations)
{
	char myString[256];
	float data;
    uint32_t start, end;
	setSeed(100);
	zigset();
	
	STCTRL = (1<<0) | (1<<2);
	wait(10);
    start = STCURR;

    for (int i=0; i<iterations; i++)
    {
        data = ziggurat();
    }
	newLine();
	
    end = STCURR;
	STCTRL = (0<<0);

    sprintf(myString, "Ziggurat method gauss generator time elapsed: %fms", (float)(start-end)/SYSTEM_CLOCK*1000);
	print(myString);
	newLine();
}

void test_ziggurat(int iterations, int seed)
{
	char myString[32];
    float data;
	setSeed(seed);
	zigset();

    for (int i=0; i<iterations; i++)
    {
        data = ziggurat();
		sprintf(myString, "%f", data);
		print(myString);
		newLine();
    }	
}

int main(void){
	int iter = 500;
	int seed = 110;
	SystemInit();
	STRELOAD = 0x00FFFFFF;
	
	time_uniform(iter);
	time_boxmuller(iter);
	time_marsaglia(iter);
	time_inverse(iter);
	time_ziggurat(iter);
	time_lfsr(iter);
	
	
	//test_uniform(iter, seed); // PASS
	//test_boxmuller(iter, seed); // PASS
	//test_marsaglia(iter, seed); // PASS
	//test_inverse(iter, seed); // PASS
	//test_ziggurat(iter, seed); // PASS
	//test_lfsr(iter, seed); // PASS
	
    return 0;
}
