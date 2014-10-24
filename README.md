LCS Tool
========

LCS Tool performs computations for the analysis of Lagrangian coherent structures. It is developed by the Nonlinear Dynamical Systems Group at [ETH Zurich](http://ETHZ.CH), led by [Prof. George Haller](http://GeorgeHaller.COM).

Demo scripts
------------

MATLAB scripts demonstrating the use of LCS Tool are in the folder named demo. To run these scripts, start MATLAB in the LCS Tool folder, then at the MATLAB prompt type, for example:

	cd demo/double_gyre/
	addpath ../..
	flow_animation

The time required to execute some of the demo scripts is given below. These times were measured on a laptop with an Intel Core i5-3360M processor with 16 gigabytes of memory:

- double\_gyre/hyperbolic\_shear\_lcs: 10 minutes
- bickley\_jet/stretchlines: 30 minutes
- ocean\_dataset/hyperbolic\_shear\_lcs\_details: 20 minutes

References
----------

Article on LCS Tool:

- K. Onu, F. Huhn and G. Haller, "LCS Tool : A Computational Platform for Lagrangian Coherent Structures". arXiv: [1406.3257 \[nlin.CD\]](http://arxiv.org/abs/1406.3527).

The algorithms used in LCS Tool are based on methods from the following publications:

- G. Haller and F. J. Beron-Vera, "Coherent Lagrangian vortices: the black holes of turbulence". _Journal of Fluid Mechanics_ 731 (2013), DOI: [10.1017/jfm.2013.391](http://dx.doi.org/10.1017/jfm.2013.391).
- M. M. Farazmand and G. Haller, "Attracting and repelling Lagrangian coherent structures from a single computation", _Chaos_ 23 (2013), DOI: [10.1063/1.4800210](http://dx.doi.org/10.1063/1.4800210).
- G. Haller and F. J. Beron-Vera, "Geodesic Theory of Transport Barriers in Two-Dimensional Flows", _Physica D: Nonlinear Phenomena_ 241 (2012), DOI: [10.1016/j.physd.2012.06.012](http://dx.doi.org/10.1016/j.physd.2012.06.012).
- M. M. Farazmand and G. Haller, "Computing Lagrangian coherent structures from their variational theory", _Chaos_ 22 (2012), DOI: [10.1063/1.3690153](http://dx.doi.org/10.1063/1.3690153).
