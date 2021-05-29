#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "../headers/generators.h"

int main(void)
{
    float data;
    int niter = 1000000;
    struct gaussGenState state;
    int taps[3] = {0, 3, 31};
    struct genState _generator;
    initializeGenerator(&_generator, 32, 8, taps, 127.63402f, 73.727015f);
    seedGenerator(&_generator, ((long long)1<<32) - 1);
    clock_t start, end;

    initializeGauss(&state);

    start = clock();
    for(int i=0; i<niter; i++)
    {
        data = gauss(&state);
    }
    end = clock();

    printf("Marsaglia time elapsed: %fms\n", (float)(end-start)/CLOCKS_PER_SEC);

    start = clock();
    for(int i=0; i<niter; i++)
    {
        data = gaussbm(&state);
    }
    end = clock();

    printf("BoxMuller time elapsed: %fms\n", (float)(end-start)/CLOCKS_PER_SEC);
    

    start = clock();
    for(int i=0; i<niter; i++)
    {
        data = gaussInv();
    }
    end = clock();

    printf("Inverse time elapsed: %fms\n", (float)(end-start)/CLOCKS_PER_SEC);
    
    zigset();
    start = clock();
    for(int i=0; i<1000000; i++)
    {
        data = ziggurat();
    }
    end = clock();

    printf("Ziggurat time elapsed: %fms\n", (float)(end-start)/CLOCKS_PER_SEC);

    start = clock();
    for(int i=0; i<niter; i++)
    {
        data = uniform();
    }
    end = clock();

    printf("Uniform time elapsed: %fms\n", (float)(end-start)/CLOCKS_PER_SEC);

    start = clock();
    for(int i=0; i<niter; i++)
    {
        data = generate(&_generator);
    }
    end = clock();

    printf("LFSR time elapsed: %fms\n", (float)(end-start)/CLOCKS_PER_SEC);
    
    return 0;
}