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

#define BUFFER_SIZE		36

void BufferFullHandler(u8 *buff, int bufflen, unsigned int *received, float *meas, int measLen);

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
    initializeFilter(&particles, &copyParticles, _particles, _copyParticles, weights, numParticles, resetWeights);
	initializeGenerator(&_State, 32, 8, taps, 517.21916f, 150.88634811f);
	seedGenerator(&_State, ((long long)1<<32) - 1);
	
	//initializeGauss(&_State);
	
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
			predictUpdate(particles, weights, meas, pitchRoll, &sum_weights, numParticles, predict_std, t, &_State);	
			
			// Resample + Estimate step
			resampleEstimate(&particles, &copyParticles, weights, sum_weights, numParticles, resetWeights, filterState);
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
