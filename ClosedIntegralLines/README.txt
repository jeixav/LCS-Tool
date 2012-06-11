cd Y:\sbarakat\Projects\Haller\lcs_tool\flow_templates\

doubleGyre = double_gyre;

cd ..

doubleGyre = shear_lcs_script(doubleGyre);

lapacklib = fullfile(matlabroot, 'extern', 'lib', 'win64', 'microsoft', 'libmwlapack.lib');
blaslib = fullfile(matlabroot, 'extern', 'lib', 'win64', 'microsoft', 'libmwblas.lib');
mex('-v', '-largeArrayDims', 'ClosedIntegralLines\compute_closed_shearline.cpp', 'ClosedIntegralLines\RegularGrid.cpp', blaslib, lapacklib, ...
	'-IC:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v4.0/include/', ...
	'-ID:/NVIDIA_GPU_Computing_SDK_40_64/C/common/inc/', ...
	'-ID:/NVIDIA_GPU_Computing_SDK_40_64/OpenCL/common/inc/', ...
	'-ID:/NVIDIA_GPU_Computing_SDK_40_64/shared/inc/', ...
	'-ID:/Libraries/gsl-1.8-src/src/gsl/1.8/gsl-1.8/', ...
	'D:/Libraries/gsl-1.8-src/src/gsl/1.8/gsl-1.8/VC8/libgslcblas/Release-StaticLib_x64/libgslcblas.lib', ...
	'D:/Libraries/gsl-1.8-src/src/gsl/1.8/gsl-1.8/VC8/libgsl/Release-StaticLib_x64/libgsl.lib');
doubleGyre = shear_lcs_script_closed(doubleGyre);











cd flow_templates\
doubleGyre = double_gyre;
cd ..
doubleGyre = shear_lcs_script(doubleGyre);