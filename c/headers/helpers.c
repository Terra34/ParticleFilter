#include "xparameters.h"
#include "m1_for_arty.h"

void wait(int num){
	for (int i=0; i < num; i++)
		i++;
}

void newLine(){
	print("\r\n");
}

float pow2(float x)
{
	return x*x;
}