#include "mex.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <windows.h>

 void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
     //input
     int *c,B; //input sequence and length   
     int BG;    // base graph
     int Zc;    // lifting size
     int *No;   // number of shifting values for an element in base graph of generator matrix
     int *R;    // shift values
     int *p;    // initial positions of shift values
     int Kb;    // information bits
     
     //output
     int *v;        
     
     //variable
     int ncols, nrows; 
     int i1,i2,i3,i4,t;
     
  int i;
   long long int Elapsed;  
     // get input
     c=(int *)mxGetData(prhs[0]);    
     B=mxGetN(prhs[0]);  
     BG=mxGetScalar(prhs[1]);
     Zc=mxGetScalar(prhs[2]);
     No=(int *)mxGetData(prhs[3]);
     R=(int *)mxGetData(prhs[4]);
     p=(int *)mxGetData(prhs[5]);
     Kb=mxGetScalar(prhs[6]);
          
     if (BG==1)
     {
         nrows=46;
         ncols=22;
     }
     else if (BG==2)
     {
         nrows=42;
         ncols=10;
     }
     
   plhs[0] = mxCreateNumericMatrix(1,(Kb+nrows) * Zc,mxINT32_CLASS,mxREAL);
   v = (int *)mxGetData(plhs[0]);
   //v=mxMalloc(sizeof(int) * (Kb+nrows) * Zc);
   memset(v,0,sizeof(int) *(Kb+nrows) * Zc);  
   
//    Elapsed = GetTickCount();
// for (i=0;i<20;i++)
// {
   //parity check part of codeword
     for (i1=0,t=Kb*Zc; i1 < nrows; i1++)
     {
        for (i2=0; i2 < Zc; i2++)
        {           
            for (i3=0; i3 < Kb; i3++)
            {
                for (i4=0; i4 < No[i1 * ncols + i3]; i4++)
                {
                    v[t] = v[t] + c[ i3*Zc + (R[ p[i1 * ncols + i3]+i4 ] + i2 + Zc) % Zc ];   
                                                 //start pointer   %element %shift due to Z
                }
            }
            //v[t]=v[t]%2;
            v[t]=v[t]&1;
            t++;
        }
     }   

   //information part of codeword
     memcpy(v,c,B*sizeof(int));              
//}

// mexPrintf("That took %lld milliseconds", GetTickCount() - Elapsed);
       return;
 }