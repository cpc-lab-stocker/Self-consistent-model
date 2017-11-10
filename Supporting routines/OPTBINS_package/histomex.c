/* function [counts, binwidth, centers] = histo_mex[Y,M]*/



#include "mex.h"		/*ncludes mex handlers*/


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{ 

/*VARIABLE DECLARATIONS***********/

int        i = 0, j = 0, k = 0, l = 0;
int        y_length;
int        m_length;

int        m_max = 0;
int        m     =0;

int        index =0;

int        y_row;
int        y_col;
int        edges_num_el;  
int        m_dim;         
int        y_dim;



double*    y_data;
double*    m_data;

mxArray*   ymin_array;
mxArray*   ymax_array;
mxArray*   edges_array;
mxArray*   binwidth_array;
mxArray*   counts_array;
mxArray*   factors_array;
mxArray*   bin_array;
mxArray*   centers_array;
mxArray*   temp_array[3];

double*    ymin_data;
double*    ymax_data;
double*    edges_data;
double*    binwidth_data;
double*    factors_data;
double*    bin_data;
double*    counts_data;
double*    temp_data[3];

int*      counts_ptr;
int       counts_dim[20];

/*CHECK PARAMETER ARGUMENTS***********/
if(nrhs != 2)
  mexErrMsgTxt("Requires two input arguments.\n");
  
if(nlhs != 3)
  mexErrMsgTxt("Requires three ouput arguments.\n");
  

  
/*GET mxARRAY INFORMATION************************/

y_length = mxGetNumberOfElements(prhs[0]);
m_length = mxGetNumberOfElements(prhs[1]);


y_row = mxGetM(prhs[0]);
y_col = mxGetN(prhs[0]);

if(m_length != y_row)
   mexErrMsgTxt("The dimensions of the data and the bin number array disagree.\n");

y_data = mxGetPr(prhs[0]);
m_data = mxGetPr(prhs[1]);

y_dim = mxGetNumberOfDimensions(prhs[0]);
m_dim = mxGetNumberOfDimensions(prhs[1]);


y_row = mxGetM(prhs[0]);
y_col = mxGetN(prhs[0]);


/*FIND THE EXTREMA ******************************/
ymin_array = mxCreateDoubleMatrix(1,y_row, mxREAL);
ymax_array = mxCreateDoubleMatrix(1,y_row, mxREAL);

ymin_data = mxGetPr(ymin_array);
ymax_data = mxGetPr(ymax_array);


for(i = 0; i < y_row; i++)
    {
    ymin_data[i] = y_data[j];
    for(; j < y_length; j += y_row)
       {
       if(y_data[j] < ymin_data[i])
          ymin_data[i] = y_data[j];
          
       if(y_data[j] > ymax_data[i])
          ymax_data[i] = y_data[j];
       }
       j = i + 1;

     }


/*DEfINE BIN EDGES*************************/
/*the extreme bins go off to inifinity*/

for(i = 0;i < m_length; i++)
    if(m_data[i] > m_max)
       m_max = m_data[i];

binwidth_array  = mxCreateDoubleMatrix(1,y_row, mxREAL);
binwidth_data   = mxGetPr(binwidth_array);

edges_array = mxCreateDoubleMatrix(y_row,m_max, mxREAL);
edges_data = mxGetPr(edges_array);


for(i = 0; i < y_row; i++)
    binwidth_data[i] = (ymax_data[i] - ymin_data[i])/m_data[i];
 

 
for(i = 0; i < y_row; i++ )
   {
    m = m_data[i] - 1;
    for(j= 0; j < m ; j++ )
        edges_data[i + j*y_row] = (binwidth_data[i])*(j+1) + ymin_data[i];
    edges_data[i + j*y_row ] = mxGetInf();
       
   }

/*obtain the counts ******************************/

for(i = 0; i <m_length; i++)
    counts_dim[i] = (int)m_data[i];
    
counts_ptr = counts_dim;


if (y_row == 1)
   {
   counts_dim[1] = (int)m_data[0];
   counts_dim[0] = 1;
   counts_array = mxCreateNumericArray(2, counts_ptr, mxDOUBLE_CLASS, mxREAL);
   }
else
  counts_array =  mxCreateNumericArray(m_length,counts_ptr, mxDOUBLE_CLASS, mxREAL);
  
  counts_data = mxGetPr(counts_array);
   
   
   
factors_array = mxCreateNumericMatrix( 1, y_row, mxDOUBLE_CLASS, mxREAL);
factors_data = mxGetPr(factors_array);

m_max = 0;

for(i = 0; i < y_row; i++)
    factors_data[i] = 1.0;
    

for(i = 0; i < (y_row - 1); i++)
  for(j = i + 1; j < y_row; j++)
      factors_data[j] = factors_data[j]*m_data[i]; 
 
              l = 0;
for(i = 0; i < y_col ; i++)
{
 
 bin_array = mxCreateDoubleMatrix(1, y_row, mxREAL);
 bin_data  = mxGetPr(bin_array);
 
 for(j = 0; j < y_row; j++ )
 {
    
  
    for(m = 0; m < m_data[j]; m++)
      {     
          
   if( y_data[j + i*y_row] <= edges_data[j + y_row*m] )
     {
      bin_data[j] = m +1;
      break;
     }
     
     
   }  

 }
 for(k = 0; k <  y_row ; k++)
    index =  index + factors_data[k]*(bin_data[k] - 1);
    index++;
     
 counts_data[index-1] = counts_data[index-1] +1;
 index = 0; 
 mxDestroyArray(bin_array);     
}

/*define the bin centers*/
centers_array = mxCreateCellMatrix(1,m_length);

for( i = 0; i < y_row; i++)
{
temp_array[i] = mxCreateDoubleMatrix(1,(int)m_data[i],mxREAL);
temp_data[i] = mxGetPr(temp_array[i]);

for(j = 0; j < (int)m_data[i]; j++)
   temp_data[i][j] = binwidth_data[i]*(j) + ymin_data[i] + binwidth_data[i]/2;
   
mxSetCell(centers_array,i,temp_array[i]);
}

/*return the computed arrays*/
      
plhs[0] = counts_array;
plhs[1] = binwidth_array;
plhs[2] = centers_array;

mxDestroyArray(ymin_array);
mxDestroyArray(ymax_array);
mxDestroyArray(edges_array);
mxDestroyArray(factors_array);




}
/*END MEX FUNCTION*/
