#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "../headers/generators.h"

void time_generator()
{
    int iterations = 50000;
    int taps[3] = {0, 3, 31};
    struct genHelper _generator;
    initializeGenerator(&_generator, 32, 8, taps, 127.63402f, 73.727015f);
    seedGenerator(&_generator, ((long long)1<<32) - 1);
    float data;
    clock_t start, end;
	
    start = clock();

    for (int i=0; i<iterations; i++)
    {
        data = generate(&_generator);
    }
	
    end = clock();

    printf("LFSR generator clocks elapsed: %fms\n", (float)(end-start)/CLOCKS_PER_SEC);
    return;
}

void time_gauss()
{
    int iterations = 50000;
    struct gaussGenHelper _normState;
    float data;
    clock_t start, end;

    initializeGauss(&_normState);

    srand(100);

    start = clock();
    for (int i=0; i<iterations; i++)
    {
        data = gauss(&_normState);
    }
    end = clock();

    printf("Gauss generator clocks elapsed: %fms\n", (float)(end-start)/CLOCKS_PER_SEC);
}

int main(void){
	time_gauss();
	time_generator();

    return 0;
}
