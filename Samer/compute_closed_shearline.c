/* $Revision: 1.8.6.4 $ */
/*=========================================================
 * convec.c
 * example for illustrating how to use pass complex data 
 * from MATLAB to C and back again
 *
 * convolves  two complex input vectors
 *
 * This is a MEX-file for MATLAB.
 * Copyright 1984-2011 The MathWorks, Inc.
 *=======================================================*/
#include "mex.h"
#include "matrix.h"

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector_types.h>
#include <vector_functions.h>

/* computational subroutine */
void convec( double *xr, double *xi, size_t nx,
             double *yr, double *yi, size_t ny,
             double *zr, double *zi)
{
    mwSize i,j;
  
    zr[0]=0.0;
    zi[0]=0.0;
    /* perform the convolution of the complex vectors */
    for(i=0; i<nx; i++) {
        for(j=0; j<ny; j++) {
            *(zr+i+j) = *(zr+i+j) + *(xr+i) * *(yr+j) - *(xi+i) * *(yi+j);
            *(zi+i+j) = *(zi+i+j) + *(xr+i) * *(yi+j) + *(xi+i) * *(yr+j);
        }
    }
}

/* The gateway routine. */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
	
	
	double* timespan = mxGetPr(prhs[0]);
	double* domain = mxGetPr(prhs[1]);
	UINT64_T* resolution = (UINT64_T*) mxGetData(prhs[2]);
	double* eigvals = mxGetPr(prhs[3]);
	double* eigvecs = mxGetPr(prhs[4]);
	int xDim, yDim;

	// get dimensions of any array
	//xDim = (int) mxGetM(prhs[2]);
	//yDim = (int) mxGetN(prhs[2]);
	//mexPrintf("x Dimensions = %d.\n",xDim);
	//mexPrintf("y Dimensions = %d.\n",yDim);
	
	// test type of array data
	//mexPrintf("test = %d.\n", mxIsDouble(prhs[4]));
   
	// get number of parameters
	//mexPrintf("Number of parameters is %d\n", nrhs);
	
	
	
	mexPrintf("%lu %lu\n", resolution[0], resolution[1]);
	//mexPrintf("%d %d\n", resolution[0], resolution[1]);
    /*double  *xr, *xi, *yr, *yi, *zr, *zi;
    size_t rows, cols;
    size_t nx, ny;
    
    // check for the proper number of arguments 
    if(nrhs != 2)
      mexErrMsgIdAndTxt( "MATLAB:convec:invalidNumInputs",
              "Two inputs required.");
    if(nlhs > 1)
      mexErrMsgIdAndTxt( "MATLAB:convec:maxlhs",
              "Too many output arguments.");
    //Check that both inputs are row vectors
    if( mxGetM(prhs[0]) != 1 || mxGetM(prhs[1]) != 1 )
      mexErrMsgIdAndTxt( "MATLAB:convec:inputsNotVectors",
              "Both inputs must be row vectors.");
    rows = 1; 
    // Check that both inputs are complex
    if( !mxIsComplex(prhs[0]) || !mxIsComplex(prhs[1]) )
      mexErrMsgIdAndTxt( "MATLAB:convec:inputsNotComplex",
              "Inputs must be complex.\n");
  
    // get the length of each input vector 
    nx = mxGetN(prhs[0]);
    ny = mxGetN(prhs[1]);


    // get pointers to the real and imaginary parts of the inputs 
    xr = mxGetPr(prhs[0]);
    xi = mxGetPi(prhs[0]);
    yr = mxGetPr(prhs[1]);
    yi = mxGetPi(prhs[1]);
  
    // create a new array and set the output pointer to it 
    cols = nx + ny - 1;
    plhs[0] = mxCreateDoubleMatrix( (mwSize)rows, (mwSize)cols, mxCOMPLEX);
    zr = mxGetPr(plhs[0]);
    zi = mxGetPi(plhs[0]);

    // call the C subroutine 
    convec(xr, xi, nx, yr, yi, ny, zr, zi);*/

    return;
}
