#ifndef GENERATORS_H_
#define GENERATORS_H_

struct genState {
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

struct gaussGenState 
{
    float gauss;
	int has_gauss;
};

long long checkBit(long long num, int bit);
int arrayContains(int *array, int length, int number);
void lfsr(struct genState *state);
float generate(struct genState *state);
void initializeGenerator(struct genState *state, int N, int M, int *taps, float mean, float std);
void seedGenerator(struct genState *state, long long seed);
float gauss(struct gaussGenState *state);
void initializeGauss(struct gaussGenState *state);


#endif
