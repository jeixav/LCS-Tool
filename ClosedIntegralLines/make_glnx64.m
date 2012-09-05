userHome = getenv('HOME');

mex('-v','-largeArrayDims','CC=gcc-4.4','CXXFLAGS=\$CFLAGS -fopenmp',...
    'compute_closed_shearline.cpp','RegularGrid.cpp','-lgsl','-lgslcblas')
