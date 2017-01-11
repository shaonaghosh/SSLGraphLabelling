SSLGraphLabelling
Semi Supervised graph labelling with various methods such as Harmonic Energy Minimization, linear programming for max-flow/mon-cut and quadratic optimization for Graph Laplacian methods.

Dependency : Matlab. Uses toolbox 'nearestneighbour' by Richard Brown.

Runs  on both simulated and UCI Irvine data sets such as USPS and MNIST. One data (MNIST and USPS) accompanies this repository.
mainfile2.m is the script that calls the optimizer. The optimizer internally calls different methods of graph labelling, given a few labels such as:

1. Harmonic Inverse from the paper 'Supervised Learning using Gaussian Fields and Harmonic Functions' by Zhu and Ghahramani ICML 2003
http://www.aaai.org/Papers/ICML/2003/ICML03-118.pdf
2. Harmonic Energy Minimization using Quadratic Programming from the paper 'Predicting the Labelling of a graph by p Semi Norm interpolation' Herbster and Lever, COLT, 2009 for {p = 2} http://discovery.ucl.ac.uk/1311163/1/fomarch11.pdf#
3. Linear Programming for max flow/min cut based graph labelling investigating open question {p=1}
4. Constrained Quadratic Formulation for {p=1} and {p=1}



