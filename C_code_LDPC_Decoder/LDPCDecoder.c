#include "mex.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

 short sign(short x) {
     return (x > 0) - (x < 0);
 } 

 void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
 
 // input
 short No_iteration; 
 int *shift_value, *Col_position, *no_one_element, Z, BG,Kb; 
 double rate;   
 short *msgChannel; 
 
 //output
 char *output_estimate; 
 
 //variables
 int nrows, ncols, nEdge_base_graph; 
 int irow, iShift, iZ,p1,p2,t1,temp_row,temp_col,iCheck,iBit,iEdge; 
 int nEdge,nCheck,nBit;
 int rows_total_La1=0,rows_total_La2=0,rows_total_La3=0,*no_rows_La1, *no_rows_La2, *no_rows_La3;
 int i1, i2, i3, sum, sum1, sum2, no_punctured_columns, layer, iLayer, n=1, flag=0;
 int *idxBit, idxCheck,*degBit, *degCheck,*degBit_base_graph, *pointerCheck, *pointerBit, *pointerCheck_temp, *pointerBit_temp;
 int *idxEdge2Bit,*idxEdge2Check;
 int idxCheck_iEdge, idxBit_iEdge;
 int *no_rows;
 int rows_total;
 short *msgBit2Check, *msgCheck2Bit, *msgBit;  //variables for message passing algorithm
 char *v_estimate;      //estimated codeword  
 short min;
 short sgn;
 int isEqual;
 char syn;      //syndrome
 
//get input
 No_iteration=mxGetScalar(prhs[0]);
 shift_value=(int *)mxGetData(prhs[1]);
 nEdge_base_graph=mxGetN(prhs[1]);
 Col_position=(int *)mxGetData(prhs[2]);
 no_one_element=(int *)mxGetData(prhs[3]);
 nrows=mxGetN(prhs[3]);
 Z=mxGetScalar(prhs[4]);
 BG=mxGetScalar(prhs[5]);
 Kb=mxGetScalar(prhs[6]);
 rate=mxGetScalar(prhs[7]);
 msgChannel=(short int *)mxGetData(prhs[8]);

 nEdge=nEdge_base_graph*Z; 
 nCheck=nrows*Z;
 if (BG==1)
     ncols=68;
 else if (BG==2)
     ncols=52; 
 nBit=ncols*Z;
 
 //initial positions of pointers to check nodes
 //degrees of check nodes
 pointerCheck = mxMalloc(sizeof(int) * nCheck); 
 degCheck=mxMalloc(sizeof(int) * nCheck);
 for (i1=0,temp_row=0; i1 < nrows; i1++)
 {
     for (i2=0; i2 < Z; i2++)
     {
          degCheck[i1*Z+i2] = no_one_element[i1];   //degree equals number of 1 elements in a row
          pointerCheck[i1*Z+i2] = temp_row;         
          temp_row = temp_row + no_one_element[i1];
     }
 } 
 
 //initial positions of pointers to bit nodes
 //degrees of bit nodes
 pointerBit = mxMalloc(sizeof(int) * nBit);
 degBit=mxMalloc(sizeof(int) * nBit);
 degBit_base_graph=mxMalloc(sizeof(int) * ncols);
 memset(degBit_base_graph,0,sizeof(int) * ncols);
 for (iBit=0; iBit < nEdge_base_graph; ++iBit)
 {
    ++degBit_base_graph[Col_position[iBit]];    //number of 1 elements in a columns in base graph
 }
 for (i1=0,temp_col=0; i1 < ncols; i1++)
 {
     for (i2=0; i2 < Z; i2++)
     {
         degBit[i1*Z+i2] = degBit_base_graph[i1];   //degree equals number of 1 elements in a column
         pointerBit[i1*Z+i2] = temp_col;        
        temp_col = temp_col + degBit_base_graph[i1];
     }
 } 
  
//indice and degrees of check nodes and bit nodes 
//divide layer for message passing algorithm
 idxBit=mxMalloc(sizeof(int) * nEdge);
// idxCheck=mxMalloc(sizeof(int) * nEdge); 
 
 idxEdge2Bit=mxMalloc(sizeof(int) * nEdge); 
 idxEdge2Check=mxMalloc(sizeof(int) * nEdge);
 pointerCheck_temp = mxMalloc(sizeof(int) * nCheck);
 memcpy(pointerCheck_temp,pointerCheck,sizeof(int) * nCheck);
 pointerBit_temp = mxMalloc(sizeof(int) * nBit);
 memcpy(pointerBit_temp,pointerBit,sizeof(int) * nBit);
 
 no_rows_La1=mxMalloc(sizeof(int) * nBit);
 no_rows_La2=mxMalloc(sizeof(int) * nBit);
 no_rows_La3=mxMalloc(sizeof(int) * nBit); 
 
 for (irow=0,p1=0,t1=0,iEdge=0; irow < nrows; ++irow)   //loop for rows in base graph
 {
     temp_row=irow*Z;
     sum=0;sum1=0;sum2=0;
     
     for (iShift=0; iShift < no_one_element[irow]; ++iShift)  //loop for 1 elements in one row of base graph   
     {
         temp_col=Col_position[p1]*Z;   
         if ( ((rate==0.2)&&(BG==2)&&(Kb==10)) || ((rate==0.33)&&(BG==1)&&(Kb==22)) )   //layer , no rate matching
         {
             layer=2;
             if (Col_position[p1]==0||Col_position[p1]==1)
                 sum++;
         }
         
         else if ( ( (BG==1) && (rate>=0.33) && (rate<=0.89) )||( (BG==2) && (rate>=0.2) && (rate<=0.67) ) ) //layer, rate matching
         { 
             layer=3;
             no_punctured_columns=ceil(nBit/Z-2-Kb/rate);
             if (Col_position[p1]==0 || Col_position[p1]==1)               
             sum1++;
             if ( (Col_position[p1] >= ncols-no_punctured_columns) && (Col_position[p1] < ncols) )
             sum2++;                    
         }
         
         for (iZ=0,p2=0; iZ < Z; iZ++)       //loop for lift size  
         {
            idxBit[t1] = (shift_value[p1]+p2)%Z + temp_col;     //column positions            
            p2++;         
            idxCheck = temp_row + iZ;                       // row positions
            
            idxEdge2Check[ pointerCheck_temp[idxCheck] ] = iEdge;  //label and store the edges connecting to check nodes        
            ++pointerCheck_temp[idxCheck];          
        
            idxEdge2Bit[ pointerBit_temp[idxBit[t1]] ] = iEdge;   //label and store the edges connecting to bit nodes    
            ++pointerBit_temp[idxBit[t1]];        
        
             iEdge++;
             t1++;
         }
         
         p1++;         
     }
 
     if ( ((rate==0.2)&&(BG==2)&&(Kb==10)) || ((rate==0.33)&&(BG==1)&&(Kb==22)) )
     {     
     if (sum==1||sum==0)
         {
             for (i1=0; i1 < Z; i1++)
             {
                no_rows_La1[rows_total_La1]=irow*Z + i1;               
                rows_total_La1++;
             }
         }
      else
         {
             for (i1=0; i1 < Z; i1++)
             {
                no_rows_La2[rows_total_La2]=irow*Z + i1;               
                rows_total_La2++;
             }
         }
       }
      else if (( (BG==1) && (rate>=0.33) && (rate<=0.89) )||( (BG==2) && (rate>=0.2) && (rate<=0.67) ))
      { 
         sum=sum1+sum2;
         if (sum==1)
         {
             for (i1=0; i1 < Z; i1++)
             {
                 no_rows_La1[rows_total_La1]=irow*Z + i1;                 
                 rows_total_La1++;
             }
         }
          else if (sum==2)
         {
             if (sum1==1)
             {
                 for (i1=0; i1 < Z; i1++)
                {
                     no_rows_La2[rows_total_La2]=irow*Z + i1;                  
                     rows_total_La2++;
                }
             }
             else if (sum1==2)
             {
                  for (i1=0; i1 < Z; i1++)
                {
                     no_rows_La3[rows_total_La3]=irow*Z + i1;                  
                     rows_total_La3++;
                }
             }
          }
          else if (sum==3)
         {
              for (i1=0; i1 < Z; i1++)
                {
                     no_rows_La2[rows_total_La2]=irow*Z + i1;                   
                     rows_total_La2++;
                }
         }
         else if (sum==0)
         {
               for (i1=0; i1 < Z; i1++)
                {
                     no_rows_La3[rows_total_La3]=irow*Z + i1;                    
                     rows_total_La3++;
                }
         }
      }
 }  

 // allocate memory for message passing algorithm 
 msgBit2Check=mxMalloc(sizeof(short) * nEdge);
 msgCheck2Bit=mxMalloc(sizeof(short) * nEdge);        
 msgBit=mxMalloc(sizeof(short) * nBit);
 memset(msgCheck2Bit,0,sizeof(short) * nEdge);
 
 // initial values of LLR of bit nodes
 memcpy(msgBit,msgChannel,nBit*sizeof(short));  
 
 //message passing algorithm 
 while (n<=No_iteration)
 {
    for (iLayer=1;iLayer<=layer;iLayer++)    
     {
         if (iLayer==1)
         {
            no_rows=no_rows_La1;
            rows_total=rows_total_La1;     
         }
         else if (iLayer==2)
         {
            no_rows=no_rows_La2;
            rows_total=rows_total_La2;            
         }
         else if (iLayer==3)
         {
             no_rows=no_rows_La3;
             rows_total=rows_total_La3;
         }
       
     //message from bit nodes to check nodes
     for(i1 = 0; i1 < rows_total; ++i1)
 	{         
         for(i2 = 0; i2 < degCheck[no_rows[i1]]; ++i2)
         {            
msgBit2Check[ idxEdge2Check[ pointerCheck[no_rows[i1]] + i2] ] = msgBit[ idxBit[ idxEdge2Check[ pointerCheck[no_rows[i1]] + i2] ] ]
                                                             -msgCheck2Bit[ idxEdge2Check[ pointerCheck[no_rows[i1]] + i2] ];
       
         }
 	}
         
    //message from check nodes to bit nodes
	for(i1 = 0; i1 < rows_total; ++i1)
    {        
         for(i2 = 0; i2 < degCheck[no_rows[i1]]; ++i2)
         {
            min=32640;
            sgn=1;
            for(i3 = 0; i3 < degCheck[no_rows[i1]]; ++i3)
            {
                if (idxEdge2Check[ pointerCheck[no_rows[i1]] + i3] != idxEdge2Check[ pointerCheck[no_rows[i1]] + i2 ])
                {                  
                    sgn *=sign(msgBit2Check[ idxEdge2Check[ pointerCheck[no_rows[i1]] + i3] ]);
                    if (abs(msgBit2Check[ idxEdge2Check[ pointerCheck[no_rows[i1]] + i3] ]) < min)
                        min=abs(msgBit2Check[ idxEdge2Check[ pointerCheck[no_rows[i1]] + i3] ]);
                }
           }
           msgCheck2Bit[ idxEdge2Check[ pointerCheck[no_rows[i1]] + i2 ] ]=sgn * min;
           
          //atanh(1)=19.07, 19.07 is converted to fixed-point 9_7= 19.07*2^7=2441
           if ( msgCheck2Bit[ idxEdge2Check[ pointerCheck[no_rows[i1]] + i2 ] ] > 2441)  
                msgCheck2Bit[ idxEdge2Check[ pointerCheck[no_rows[i1]] + i2 ] ] = 2441;
            if (msgCheck2Bit[ idxEdge2Check[ pointerCheck[no_rows[i1]] + i2 ] ] < -2441)
                msgCheck2Bit[ idxEdge2Check[ pointerCheck[no_rows[i1]] + i2 ] ] = -2441;
       }
    } 
            // LLR
         v_estimate=mxMalloc(sizeof(char) * nBit);
    for(iBit = 0; iBit < nBit; ++iBit)
    {
        msgBit[iBit]=msgChannel[iBit];
        for(i1 = 0; i1 < degBit[iBit]; ++i1)
            msgBit[iBit] += msgCheck2Bit[idxEdge2Bit[ pointerBit[iBit]+i1 ] ];   
        if (msgBit[iBit]>=0)
            v_estimate[iBit]=0;
        else
            v_estimate[iBit]=1;  
    } 
//     // LLR
//     for(iBit = 0; iBit < nBit; ++iBit)
//     {
//         msgBit[iBit]=msgChannel[iBit];
//         for(i1 = 0; i1 < degBit[iBit]; ++i1)
//             msgBit[iBit] += msgCheck2Bit[idxEdge2Bit[ pointerBit[iBit]+i1 ] ];       
//     }
//          
//     //estimate codeword  
//     v_estimate=mxMalloc(sizeof(char) * nBit);
//     for(iBit = 0; iBit < nBit; ++iBit)
//     {
//         if (msgBit[iBit]>=0)
//             v_estimate[iBit]=0;
//         else
//             v_estimate[iBit]=1;      
//     }
    
    //check syndrome=0
    for (i1=0,isEqual=1; i1 < nCheck; ++i1)
    {
        syn=0;
        for(i2 = 0; i2 < degCheck[i1]; ++i2)
        {            
            //sum of bits of estimated codeword in the positions where there are 1 elements in graph
            syn += v_estimate[ idxBit[ idxEdge2Check[ pointerCheck[ i1 ] + i2  ] ] ];         
        }
        syn=syn%2;   
        if (syn != 0)  
        {
            isEqual=0;
            break;
        }
    }       
   
    if (isEqual==1)   //syndrome=0, break
    {
        flag=1;
        break;
    }
    
   }
    
     if (flag==1)
     {
         break;
     }
    
    n++;    
 }
 
 plhs[0] = mxCreateNumericMatrix(1,nBit,mxINT8_CLASS,mxREAL);
 output_estimate = (signed char *)mxGetData(plhs[0]);
 memcpy(output_estimate,v_estimate,nBit*sizeof(char));
 
 //free memory
 mxFree(no_rows_La1);
mxFree(no_rows_La2);
mxFree(no_rows_La3);
mxFree(idxBit);
//mxFree(idxCheck);
mxFree(degBit);
mxFree(degCheck);
mxFree(pointerCheck);
mxFree(pointerBit);
mxFree(pointerCheck_temp);
mxFree(pointerBit_temp);
mxFree(idxEdge2Bit);
mxFree(idxEdge2Check);
mxFree(msgBit2Check);
mxFree( msgCheck2Bit);
mxFree(msgBit);
mxFree(v_estimate);
                
     return;
 }