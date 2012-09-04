userHome = getenv('HOME');

mex('-v','-largeArrayDims','CC=gcc-4.4','compute_closed_shearline.cpp',...
    'RegularGrid.cpp','-lgsl','-lgslcblas')
