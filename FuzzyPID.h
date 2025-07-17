#ifndef FUZZYPID_H
#define FUZZYPID_H

#include <stdint.h>


//#include "pid.h"

//PID
typedef struct
{
 float setVaule;  
 float deta_kp; 
 float date_ki;  
 float date_kd;  
 float lasterror; 
 float preerror;

 float maximum; 
 float minimum;  
 
 float qKp;    
 float qKi;      
 float qKd;    
}FUZZYPID;


extern FUZZYPID FPID;

void Fuzzytrans(float _Set_Vaule,float _Measure_Vaule,float pre_Measure_Vaule);
void Calculatetimer(float cal_time, uint8_t *Array, uint16_t Length);
void DWT_Init(void);
uint32_t DWT_GetUs(void);
	
#endif

