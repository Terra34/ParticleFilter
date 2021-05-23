#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "../headers/generators.h"
#include "../headers/particleFilter.h"

#define PI 3.14159265359
#define TWO_PI 6.28318530718


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
	//struct gaussGenHelper _State;
	
	float meas[9];

    // Filter initialization
    srand(100);
    initializeFilter(&particles, &copyParticles, _particles, _copyParticles, weights, numParticles, resetWeights);
	initializeGenerator(&_State, 32, 8, taps, 127.63402f, 73.727015f);
	seedGenerator(&_State, ((long long)1<<32) - 1);
	
	//init_gauss(&_State);
	
    pitchRoll = getPitchRoll(meas);
    
    // Predict + Update step
    predictUpdate(particles, weights, meas, pitchRoll, &sum_weights, numParticles, predict_std, t, &_State);	
    
    // Resample + Estimate step
    resampleEstimate(&particles, &copyParticles, weights, sum_weights, numParticles, resetWeights, filterState);

    yaw = getYaw(filterState, meas+6);
    
    
    printf("%f,%f,%f\n", filterState[0], filterState[1], yaw);
}
