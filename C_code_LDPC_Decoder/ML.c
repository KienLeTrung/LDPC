#include "mex.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
double *r;
int nBit, no_punctured_columns, Z;
short *p_channel;

short square,ratio;
int temp,i1,i2;

 r = mxGetPr(prhs[0]); 
// nBit = mxGetScalar(prhs[1]);
 nBit = mxGetN(prhs[0]);
 no_punctured_columns = mxGetScalar(prhs[1]);
 Z = mxGetScalar(prhs[2]);
 
 nBit= nBit + (2 + no_punctured_columns)*Z;
 
 plhs[0] = mxCreateNumericMatrix(1,nBit,mxINT16_CLASS,mxREAL);
 p_channel = (short int *)mxGetData(plhs[0]);
 
  
//  temp = (int) (std*128) * (std*128);
//  temp = temp + (1<<6);
//  temp = temp >> 7;
//  square=temp;
//  //mexPrintf("%d\n",square);
//  
//  temp = (int) (1.41421356*128*128);
//  //mexPrintf("%d\n",temp);
//  temp=temp + square/2; 
//  temp=temp/square;
//  ratio=temp;
//  //mexPrintf("%d\n",ratio);
 
for (i1=0, i2=0; i1<nBit; i1++)
{
    if ( (i1 < 2*Z) || ( (i1 >= nBit-no_punctured_columns*Z) && (i1 < nBit) ) )
        p_channel[i1]=0;
    else if ( (i1 >= 2*Z) && (i1 < nBit-no_punctured_columns*Z) )
    {
//         temp = (int) (r[i2]*128) * ratio;
//         temp = temp + (1<<6);
//         temp = temp >> 7;
        p_channel[i1] = (int) (r[i2]*128);
        i2++;
    }    
}
    
 return;
 }