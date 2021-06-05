/***************************** Include Files *******************************/

#include <stdio.h>
#include <stdlib.h>
#include "xparameters.h"
#include "m1_for_arty.h" // Project specific header
#include "xil_printf.h"
#include "particleFilter.h"

/*************************** Global Variables ******************************/

#include "PmodNAV.h"
PmodNAV nav;
int shouldSample;

/***************************** Function Definitions ************************/

#define STCTRL (*( ( volatile unsigned long *) 0xE000E010 ))
#define STRELOAD (*( ( volatile unsigned long *) 0xE000E014 ))
#define STCURR (*( ( volatile unsigned long *) 0xE000E018 ))
#define SBIT_ENABLE 0
#define SBIT_TICKINT 1
#define SBIT_CLKSOURCE 2
#define RELOAD_VALUE 9999999 

void SysTick_Handler(void)
{
	shouldSample=1;
}

void readNavData(PmodNAV *nav, float *state)
{
	state[0] = nav->acclData.X;
	state[1] = nav->acclData.Y;
	state[2] = nav->acclData.Z;
	state[3] = nav->gyroData.X;
	state[4] = nav->gyroData.Y;
	state[5] = nav->gyroData.Z;
	state[6] = nav->magData.X;
	state[7] = nav->magData.Y;
	state[8] = nav->magData.Z;
}

int main(void) {
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
	
	float meas[9];
	char myString[256];
	NAV_begin(&nav,XPAR_PMODNAV_0_AXI_LITE_GPIO_BASEADDR,XPAR_PMODNAV_0_AXI_LITE_SPI_BASEADDR);
	NAV_Init(&nav);
	
	STRELOAD = RELOAD_VALUE;
	STCTRL = (1<<SBIT_ENABLE) | (1<<SBIT_TICKINT) | (1<<SBIT_CLKSOURCE);
	
	srand(101);
    setSeed(101);
    initializeFilter(&particles, &copyParticles, _particles, _copyParticles, weights, numParticles, resetWeights);
    zigset();
	
	while(1){
		if (shouldSample==1)
		{
			NAV_GetData(&nav);
			readNavData(&nav, meas);
			pitchRoll = getPitchRoll(meas);
			predictUpdateZiggurat(particles, weights, meas, pitchRoll, &sum_weights, numParticles, predict_std, t);
			resampleEstimate(&particles, &copyParticles, weights, sum_weights, numParticles, resetWeights, filterState);
			yaw = getYaw(filterState, meas+6);
			sprintf(myString, "r%frp%fpy%fy\n", filterState[0], filterState[1], yaw);
			print(myString);
			shouldSample=0;
		}
	}
	return 0;
}