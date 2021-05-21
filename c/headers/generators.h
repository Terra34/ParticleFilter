#ifndef GENERATORS_H_
#define GENERATORS_H_

struct genHelper {
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

struct gaussGenHelper 
{
    float gauss;
	int has_gauss;
};

long long checkBit(long long num, int bit);
int arrayContains(int *array, int length, int number);
void lfsr(struct genHelper *helper);
float generate(struct genHelper *helper);
void initializeGenerator(struct genHelper *helper, int N, int M, int *taps, float mean, float std);
void seedGenerator(struct genHelper *helper, long long seed);
float gauss(struct gaussGenHelper *state);
void initializeGauss(struct gaussGenHelper *state);


#endif