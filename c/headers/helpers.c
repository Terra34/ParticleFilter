#include "xil_printf.h"
#include "helpers.h"

void wait(int num){
	for (int i=0; i < num; i++)
		i++;
}

void newLine(){
	print("\r\n");
}

float BytesToFloat(u8 *buff){
	// Assumes little endian in buff
	uint32_t temp;
	
	temp = ((buff[3] << 24) |
			(buff[2] << 16) |
			(buff[1] << 8)  |
			buff[0]);
	
	return *((float *) &temp);
}
