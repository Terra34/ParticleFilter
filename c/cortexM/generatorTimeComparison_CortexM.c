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

// Xilinx specific headers
#include "xparameters.h"
#include "m1_for_arty.h" // Project specific header
#include "xil_printf.h"
#include "generators.h"

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

void newLine()
{
	print("\r\n");
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

void time_gauss()
{
	char myString[256];
    int iterations = 1000;
    struct gaussGenHelper _normState;
    float data;
    uint32_t start, stop;

    initializeGauss(&_normState);

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
