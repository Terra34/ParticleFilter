#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdint.h>

static unsigned int jz, jsr = 123456789;

#define PI 3.14159265359

#define SHR3 (jz=jsr, jsr^=(jsr<<13), jsr^=(jsr>>17), jsr^=(jsr<<5),jz+jsr)
#define LOW_MASK (0x0000FFFF)
#define FIXED_ACCURACY (0.005493164063)
#define LFSR_RIGHT_SHIFT (3)

int16_t uniform16(void){
    /*
        Generate uniform distribution (for particle generation)
    */
    return (int16_t)(LOW_MASK&SHR3);
}

int16_t lfsr2(void){
    /*
        Gaussov RNG - basically CLT with 4 uniform variates
        right shift determines sigma
    */
    unsigned int n1, n2;
    n1 = SHR3;
    n2 = SHR3;

    return ((int16_t)((((n1 >> 16) + (LOW_MASK&n1) + (n2 >> 16) + (LOW_MASK&n2)) >> 2) - 0x00007FFD)) >> LFSR_RIGHT_SHIFT;
}

int searchSorted(int* arr, int r, int numParticles){
    /*
        arr => Array that needs to be searched
        r => Value to search by
    */
    int L = 0;
    int R = numParticles - 1;
    int m;

    while (R >= L)
    {
        m = ((L+R) >> 1);
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

float pow2f(float num){
    return num*num;
}

int pow2(int num){
    return (num*num);
}

int16_t float_to_fixed(float input){
    // assumes [-180, 180] data
    return (int16_t)(round(input / FIXED_ACCURACY));
}

float fixed_to_float(int16_t input){
    return ((float)input * FIXED_ACCURACY);
}

int16_t * getPitchRoll(float* meas){
    // meas is acc_data => [accx accy accz]
    static int16_t pitchRoll[2];
    float temp;

    // roll
    temp = atan2(-1. * *(meas + 1), sqrt(pow2f(*meas) + pow2f(*(meas + 2)))) * 180 / PI;
    pitchRoll[0] = float_to_fixed(temp);

    // pitch
    temp = atan2(*meas, sqrt(pow2f(*(meas + 1)) + pow2f(*(meas + 2)))) * 180 / PI;
    pitchRoll[1] = float_to_fixed(temp);

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

void initializeFilter(	int16_t **particles, int16_t **copyParticles, int16_t *_particles, int16_t *_copyParticles, 
						uint16_t *weights, int numParticles){
	*particles = _particles;
    *copyParticles = _copyParticles;
    for (int i=0; i<numParticles; i++){
        _particles[i*2] = uniform16();
        _particles[i*2 + 1] = uniform16();
    }
}

void predictUpdate(	int16_t *particles, uint16_t *weights, int16_t *gyro, int16_t *pitchRoll,
					unsigned int *sum_weights, int numParticles){
    /*
        gyro is gyroscope data scaled by sampling interval (for example 0.1s)
    */
    int temp;
	*sum_weights = 0;
	for (int i = 0; i < numParticles; i++){
		particles[i*2] += gyro[0] + lfsr2();
		particles[i*2 + 1] += gyro[1] + lfsr2();

        //std::cout << particles[i*2] << ", " << particles[i*2 + 1];

        // all ones divided by upper 16 bits of the result!
        temp = (pow2f((int)particles[i*2] - pitchRoll[0]) 
					            + pow2f((int)particles[i*2 + 1] - pitchRoll[1])) >> 16;
		// as only the upper 16 bits are considered, a lot of the particles are treated as equally similar even though they might not be (weights diversification)
		// However, if this effect needs to be further boosted, the lower bits of the result can be nulled
		// temp = temp & 0xFFFC; 
        if(temp == 0) temp++;
		weights[i] = 0xFFFF / temp;

		*sum_weights += weights[i];
	}
}

void resampleEstimate( int16_t **particles, int16_t **copyParticles, uint16_t *weights,
                       unsigned int sum_weights, int numParticles,
                       float *filterState){
    /*
        filterState is an int because in it we accumulate all states so it needs to be
        of bigger span than the particle representation.
    */
    int k, sum_residual, pitch, roll;
	int residuals[numParticles], indexes[numParticles];
    uint16_t resetWeight;
	int16_t *temp;

    resetWeight = (uint16_t)(sum_weights/numParticles);
    
    k = 0;

    while (weights[0] > resetWeight)
    {
        weights[0] -= resetWeight;
        indexes[k++] = 0;
    }
    sum_residual = weights[0];
    residuals[0] = weights[0];

    for (int i = 1; i < numParticles; i++){
        while (weights[i] > resetWeight)
        {
            weights[i] -= resetWeight;
            indexes[k++] = i;
        }
        sum_residual += weights[i];
        residuals[i] = weights[i] + residuals[i-1];
	}

	while (k < numParticles){
		indexes[k++] = searchSorted(residuals, SHR3%(sum_residual+1), numParticles);
	}

    pitch = 0;
	roll = 0;
	for (int i=0; i<numParticles; i++){
		(*copyParticles)[i*2] = (*particles)[indexes[i]*2];
		(*copyParticles)[i*2 + 1] = (*particles)[indexes[i]*2 + 1];
		roll += (*copyParticles)[i*2];
		pitch += (*copyParticles)[i*2 + 1];
	}
	// swap array pointers
	temp = *particles;
	*particles = *copyParticles;
	*copyParticles = temp;

	filterState[0] = ((float)roll / numParticles) * FIXED_ACCURACY;
	filterState[1] = ((float)pitch / numParticles) * FIXED_ACCURACY;
}

int main(void)
{
    // measurement variables
    float _meas[9] = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};

    // filter parameters
    float param_t = 0.1;
    int numParticles = 512;
   
    // filter variables
    float filterState[2];
    int16_t _particles[numParticles*2];
    int16_t _copyParticles[numParticles*2];
    uint16_t weights[numParticles];
    float t = param_t;
    int16_t *pitchRoll;
    float yaw;

    // helper variables
    unsigned int sum_weights;
    int16_t *particles, *copyParticles;
    int16_t gyroMeas[2];


    initializeFilter(&particles, &copyParticles, _particles, _copyParticles, weights, numParticles);
    gyroMeas[0] = float_to_fixed(param_t * _meas[3]);
    gyroMeas[1] = float_to_fixed(param_t * _meas[4]);
    pitchRoll = getPitchRoll(_meas);
    predictUpdate(particles, weights, gyroMeas, pitchRoll, &sum_weights, numParticles);
    resampleEstimate(&particles, &copyParticles, weights, sum_weights, numParticles, filterState);
    yaw = getYaw(filterState, _meas + 6);
    
    return 0;
}