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

/*******************************************************************/

#define PI 3.14159265359
#define TWO_PI 6.28318530718
#define BUFFER_SIZE		36

void BufferFullHandler(u8 *buff, int bufflen, unsigned int *received, float *meas, int measLen);

int searchSorted(float* arr, float r, int numParticles){
    /*
        arr => Array that needs to be searched
        r => Value to search by
    */
    int L = 0;
    int R = numParticles - 1;
    int m;

    while (R >= L)
    {
        m = (int) ((L+R) / 2);
        if (r > *(arr + m))
        {
            if (r <= *(arr + m + 1) || (m + 1 > numParticles - 1))
                return m+1;
            L = m + 1;
        }
        else // if (r <= *(arr + m))
        {
            if (r > *(arr + m - 1) || (m - 1 < 0))
                return m;
            R = m - 1;
        }
    }
    return -1;
}

float * getPitchRoll(float* meas){
    // meas is acc_data => [accx accy accz]
    static float pitchRoll[2];
    // roll
    pitchRoll[0] = atan2(-1. * *(meas + 1), sqrt(pow2(*meas) + pow2(*(meas + 2)))) * 180 / PI;

    // pitch
    pitchRoll[1] = atan2(*meas, sqrt(pow2(*(meas + 1)) + pow2(*(meas + 2)))) * 180 / PI;

    return pitchRoll;
}

float getYaw(float* fS, float* meas){
    /*  
        fS => filtered state - assumes degrees!
        meas => magnetometer data [magx magy magz]
    */
   float sina, sinb, cosa, cosb, XH, YH;

   sina = sin(*fS * PI / 180);
   sinb = sin(*(fS + 1) * PI / 180);
   cosa = cos(*fS * PI / 180);
   cosb = cos(*(fS + 1) * PI / 180);
   XH = *meas * cosb + *(meas+1) * sinb * sina + *(meas+2) * sinb * cosa;
   YH = *(meas+1) * cosa + *(meas+2) * sina;
   return atan2(-YH, XH) * 180 / PI;
}


int main(void){
    // filter parameters
    float param_t = 0.1f;
    int numParticles = 128;

    // sensor parameters
    float predict_std = 10.f;
   
    // filter variables
    float filterState[2];
    float _particles[numParticles][2];
    float _copyParticles[numParticles][2];
    float weights[numParticles];
    float t = param_t;
    float *pitchRoll;
    float yaw;

    // helper variables
    float (*particles)[2], (*copyParticles)[2], (*temp)[2];
    float sum_weights = 1.f;
    float resetWeights = 1.f / numParticles;
    float x;
    int k;
    float residual;
    float sum_residual;
    float residuals[numParticles];
    float w;
    int num_copies;
    int indexes[numParticles];
    struct genHelper _State;
	int taps[3] = {0, 3, 31};
	//struct gaussGenHelper _State;
	
	// Uart read data variables
	char myString[256];
	int Status;
	unsigned int ReceivedCount = 0;
	XUartLite uart;
	u8 RecvBuffer[BUFFER_SIZE];
	float meas[9];
	
	SystemInit();
	
	// Initialize UART
	Status = XUartLite_Initialize(&uart, XPAR_AXI_UARTLITE_0_DEVICE_ID);
	if (Status != XST_SUCCESS) {
		return XST_FAILURE;
	}
	
	Status = XUartLite_SelfTest(&uart);
	if (Status != XST_SUCCESS) {
		return XST_FAILURE;
	}

    // Filter initialization
    srand(101);
    particles = _particles;
    copyParticles = _copyParticles;
    for (int i=0; i<numParticles; i++){
        weights[i] = resetWeights;
        particles[i][0] = (float)rand() / RAND_MAX * 360 - 180;
        particles[i][1] = (float)rand() / RAND_MAX * 360 - 180;
    }
	initializeGenerator(&_State, 32, 8, taps, 127.63402f, 73.727015f);
	seedGenerator(&_State, ((long long)1<<32) - 1);
	
	//init_gauss(&_State);
	
	while (1){
		// wait for data from UART
		ReceivedCount += XUartLite_Recv(&uart,
							RecvBuffer + ReceivedCount,
							BUFFER_SIZE - ReceivedCount);
		
		if (ReceivedCount == BUFFER_SIZE)
		{
			BufferFullHandler(RecvBuffer, BUFFER_SIZE, &ReceivedCount, meas, 9);
			pitchRoll = getPitchRoll(meas);
			
			// Predict + Update step
			sum_weights = 0.f;
			for (int i = 0; i < numParticles; i++){
				particles[i][0] += meas[3] * t + generate(&_State) * predict_std;
				particles[i][1] += meas[4] * t + generate(&_State) * predict_std;
				
				//particles[i][0] += meas[3] * t + gauss(&_State) * predict_std;
				//particles[i][1] += meas[4] * t + gauss(&_State) * predict_std;

				if (particles[i][0] > 180)
					particles[i][0] -= 360;
				else if (particles[i][0] < -180)
					particles[i][0] += 360;

				if (particles[i][1] > 180)
					particles[i][1] -= 360;
				else if (particles[i][1] < -180)
					particles[i][1] += 360;

				weights[i] *= 1. / (pow2(particles[i][0] - pitchRoll[0]) 
							+ pow2(particles[i][1] - pitchRoll[1]));
				weights[i] += FLT_MIN;

				sum_weights += weights[i];
			}	
			
			// Resample + Estimate step
			k = 0;

			w = numParticles * weights[0] / sum_weights;
			num_copies = (int)w;
			residual = w - num_copies;
			residuals[0] = residual;
			sum_residual = residual;

			for (int j=0; j<num_copies; j++){
				indexes[k++] = 0;
			}
			for (int i = 0; i < numParticles; i++){
				w = numParticles * weights[i] / sum_weights;
				num_copies = (int)w;
				residual = w - num_copies;
				residuals[i] = residual + residuals[i-1];
				sum_residual += residual;

				for (int j=0; j<num_copies; j++){
					indexes[k++] = i;
				}
			}

			while (k < numParticles){
				x = (float)rand() * sum_residual / RAND_MAX;
				indexes[k++] = searchSorted(residuals, x, numParticles);
			}
			
			filterState[0] = 0.f;
			filterState[1] = 0.f;
			for (int i=0; i<numParticles; i++){
				copyParticles[i][0] = particles[indexes[i]][0];
				copyParticles[i][1] = particles[indexes[i]][1];
				// resample happens after every step, which means that the weights will always be 
				// 1/numParticles at this step (because we just resampled right before estimating)
				filterState[0] += copyParticles[i][0];
				filterState[1] += copyParticles[i][1];
				weights[i] = resetWeights;
			}
			// swap array pointers
			temp = particles;
			particles = copyParticles;
			copyParticles = temp;
			
			filterState[0] /= numParticles;
			filterState[1] /= numParticles;

			yaw = getYaw(filterState, meas+6);
			
			
			sprintf(myString, "%f,%f,%f", filterState[0], filterState[1], yaw);
			print(myString);
			newLine();
		}
	}
}

void BufferFullHandler(u8 *buff, int bufflen, unsigned int *received, float *meas, int measLen)
{
	char myString[128];
	int iterations;
	
	if (bufflen%4)
	{
		print("Buffer length is not a multiple of 4! No bytes handled.");
		newLine();
		*received = 0;
		return;
	}
	
	iterations = bufflen/4;
	
	if (iterations != measLen)
	{
		sprintf(myString, "Need to receive exactly %d bytes - %d received! No bytes handled.", measLen, iterations); 
		print(myString);
		newLine();
		*received = 0;
		return;
	}
	
	for (int i=0; i<iterations; i++){
		meas[i] = BytesToFloat(buff);
		buff += 4;
	}
	
	*received = 0;
}