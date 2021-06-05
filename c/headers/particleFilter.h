#ifndef PARTICLEFILTER_H_
#define PARTICLEFILTER_H_

#include "generators.h"
#define PI 3.14159265359
#define TWO_PI 6.28318530718

int searchSorted(float* arr, float r, int numParticles);
float * getPitchRoll(float* meas);
float getYaw(float* fS, float* meas);
void initializeFilter(	float **particles, float **copyParticles, float *_particles, float *_copyParticles, 
						float *weights, int numParticles, float resetWeights);
void predictUpdate(	float *particles, float *weights, float *meas, float *pitchRoll,
					float *sum_weights, int numParticles, float predict_std, float t, 
					struct genState *state);
void predictUpdateGauss(float *particles, float *weights, float *meas, float *pitchRoll,
						float *sum_weights, int numParticles, float predict_std, float t, 
						struct gaussGenState *state);
void predictUpdateZiggurat(float *particles, float *weights, float *meas, float *pitchRoll,
						float *sum_weights, int numParticles, float predict_std, float t);
void resampleEstimate(	float **particles, float **copyParticles, float *weights, 
						float sum_weights, int numParticles, float resetWeights, 
						float *filterState);

#endif
