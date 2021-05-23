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
#include "xuartlite.h"
#include "generators.h"
#include "helpers.h"
#include "particleFilter.h"

/*******************************************************************/

#define STCTRL      (*( (volatile uint32_t *) 0xE000E010 ))
#define STRELOAD    (*( (volatile uint32_t *) 0xE000E014 ))
#define STCURR      (*( (volatile uint32_t *) 0xE000E018 ))  

int main(void){
    // filter parameters
    float param_t = 0.1f;
    int numParticles = 128;

    // sensor parameters
    float predict_std = 10.f;
   
    // filter variables
    float filterState[2];
    float _particles[numParticles*2];
    float _copyParticles[numParticles*2];
    float weights[numParticles];
    float t = param_t;
    float *pitchRoll;
    float yaw;

    // helper variables
    float *particles, *copyParticles;
    float sum_weights = 1.f;
    float resetWeights = 1.f / numParticles;
    struct genState _State;
	int taps[3] = {0, 3, 31};
	//struct gaussGenState _State;
	
    char myString[256];
    SystemInit();
	STRELOAD = 0x00FFFFFF;

	float meas[9];
    uint32_t start, end;

    // Filter initialization
    srand(101);
    initializeFilter(&particles, &copyParticles, _particles, _copyParticles, weights, numParticles, resetWeights);
	initializeGenerator(&_State, 32, 8, taps, 517.21916f, 150.88634811f);
	seedGenerator(&_State, ((long long)1<<32) - 1);
	
	//initializeGauss(&_State);

    pitchRoll = getPitchRoll(meas);
    
    STCTRL = (1<<0) | (1<<2);
	wait(100);
    start = STCURR;
    
    // Predict + Update step
    predictUpdate(particles, weights, meas, pitchRoll, &sum_weights, numParticles, predict_std, t, &_State);	
    
    end = STCURR;
	STCTRL = (0<<0);

    sprintf(myString, "Predict + Update time elapsed: %fms", (float)(start-end)/SYSTEM_CLOCK*1000);
	print(myString);
	newLine();

    STCTRL = (1<<0) | (1<<2);
	wait(100);
    start = STCURR;

    // Resample + Estimate step
    resampleEstimate(&particles, &copyParticles, weights, sum_weights, numParticles, resetWeights, filterState);
    yaw = getYaw(filterState, meas+6);

    end = STCURR;
	STCTRL = (0<<0);

    sprintf(myString, "Resample + Estimate time elapsed: %fms", (float)(start-end)/SYSTEM_CLOCK*1000);
	print(myString);
	newLine();
}
