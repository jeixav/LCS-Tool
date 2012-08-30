
1) Enable the use of OpenMP with Matlab
in C:\Users\sbarakat\AppData\Roaming\MathWorks\MATLAB\R2012a\mexopts.bat
set COMPFLAGS=/openmp $COMPFLAGS

2) Select the system dataset by uncommenting in the function GetDoubleGyreFlowV
in file RegularGrid.cpp

3) Review the parameters at the beginning compute_closed_shearline.cpp

4) Compile and run the code:

clear;
setenv('OMP_NUM_THREADS','16')
maxNumCompThreads(16)
cd E:\Chaos2Rost\Projects\Haller\lcs_tool\
addpath('flow_templates')
doubleGyre = double_gyre;
doubleGyre = set_flow_resolution(uint64([2 1]*400),doubleGyre);
doubleGyre.flow.resolution = doubleGyre.resolution;
%doubleGyre = shear_lcs_script(doubleGyre);

lapacklib = fullfile(matlabroot, 'extern', 'lib', 'win64', 'microsoft', 'libmwlapack.lib');
blaslib = fullfile(matlabroot, 'extern', 'lib', 'win64', 'microsoft', 'libmwblas.lib');
mex('-v', '-largeArrayDims', 'ClosedIntegralLines\compute_closed_shearline.cpp', 'ClosedIntegralLines\RegularGrid.cpp', blaslib, lapacklib, ...
	'-IC:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v4.0/include/', ...
	'-ID:/NVIDIA_GPU_Computing_SDK_40_64/C/common/inc/', ...
	'-ID:/NVIDIA_GPU_Computing_SDK_40_64/OpenCL/common/inc/', ...
	'-ID:/NVIDIA_GPU_Computing_SDK_40_64/shared/inc/', ...
	'-ID:/Libraries/gsl-1.8-src/src/gsl/1.8/gsl-1.8/', ...
	'-ID:/Libraries/teem/include/', ...
	'-IC:/Program Files (x86)/boost/boost_1_47/', ...
	'D:/Libraries/gsl-1.8-src/src/gsl/1.8/gsl-1.8/VC8/libgslcblas/Release-StaticLib_x64/libgslcblas.lib', ...
	'D:/Libraries/gsl-1.8-src/src/gsl/1.8/gsl-1.8/VC8/libgsl/Release-StaticLib_x64/libgsl.lib', ...
	'D:/Libraries/teem/bin/Release/teem.lib', ...
	'D:/Libraries/teem/arch/win32/lib/Release/bz2.lib', ...
	'D:/Libraries/teem/arch/win32/lib/Release/png.lib', ...
	'D:/Libraries/teem/arch/win32/lib/Release/z.lib');
doubleGyre = shear_lcs_script_closed(doubleGyre);







