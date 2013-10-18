Lux
===

Test library for MCMC computations and figuring out the whole
C++/cmake/Rcpp/R mess.  This code is written as a library for
an Rcpp/R package to allow efficient posterior simulation from
a first order random walk with t-distributed steps to add the
tools to R necessary for inference on heavy-tailed time-series
and 1-dimensional movement models.  This code only simulates
states from the hidden markov model conditional on the parameter,
inference on the parameters conditional on the states can be
carried out in R with any LM/GLM/GLMM/GAM package as part of a
Gibbs sampler.  The key is making the start-up/update overhead
low once the problem is set up, which is true for Lux but
more difficult for some of the other R packages.

Notes about the code quality: 
1) For now, it is not always correct so get in touch before 
   you choose to use it.  
2) This was my chance to learn C++ in a numerical context where
   I knew from previous attempts in R that speed was an issue.
   I've learned STL/templates along the way and the approach I
   take has changed from class to class.  It's going to be
   messy until I have the chancen to go back and rewrite/test.


