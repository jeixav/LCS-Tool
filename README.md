LCS Tool
========

LCS Tool performs computations for the analysis of Lagrangian coherent structures. It is based on work of <a href="http://www.georgehaller.com/">George Haller’s Nonlinear Dynamics Group.</a>

Demo scripts
------------

MATLAB scripts demonstrating the use of LCS Tool are in the folder named demo. To run these scripts, start MATLAB in the LCS Tool folder, then at the MATLAB prompt type, for example:

	cd demo/double_gyre/
	addpath ../..
	matlabpool open
	flow_animation

The time required to execute some of the demo scripts is given below. These times were measured on a laptop with an Intel Core i5-3360M processor with 16 gigabytes of memory:

- double\_gyre/hyperbolic\_shear\_lcs_details: 11 minutes
- bickley\_jet/stretchlines: 26 minutes
- ocean\_dataset/demo\_script_ocean: 18 minutes

Examples of images produced by these scripts are shown below.

![Output of double_gyre/hyperbolic_shear_lcs_details, attracting LCSs.](https://raw.github.com/jeixav/LCS-Tool/master/demo/double_gyre/hyperbolic_shear_lcs_details_attracting.png "Double gyre forward time LCS analysis.")

![Output of bickley_jet/stretchlines.m, forward time stretchlines.](https://raw.github.com/jeixav/LCS-Tool/master/demo/bickley_jet/stretchlines_forward.png "Bickley jet forward time stretchlines.")

![Output of ocean_dataset/hyperbolic_shear_lcs.m, backward-time LCSs.](https://raw.github.com/jeixav/LCS-Tool/master/demo/ocean_dataset/hyperbolic_shear_lcs_backward.png "Ocean dataset forward time LCS analysis.")
