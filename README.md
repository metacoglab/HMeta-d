HMeta-d
===

**Hierarchical meta-d' model (HMeta-d)**

Steve Fleming
stephen.fleming@ucl.ac.uk 

This MATLAB toolbox implements the meta-d’ model (Maniscalco & Lau, 2012) in a hierarchical Bayesian framework using Matlab and JAGS, a program for conducting MCMC inference on arbitrary Bayesian models. A paper with more details on the method and the advantages of estimating meta-d’ in a hierarchal Bayesian framework is available here https://academic.oup.com/nc/article/doi/10.1093/nc/nix007/3748261/HMeta-d-hierarchical-Bayesian-estimation-of.

For a more general introduction to Bayesian models of cognition see Lee & Wagenmakers, Bayesian Cognitive Modeling: A Practical Course http://bayesmodels.com/

The model builds on work by Michael Lee on Bayesian estimation of Type 1 SDT parameters: https://link.springer.com/article/10.3758/BRM.40.2.450 

The code is designed to work “out of the box” without much coding on the part of the user, and it receives data in the same format as Maniscalco & Lau’s toolbox, allowing easy switching and comparison between the two.

1) To get started, you need to first ensure JAGS (an MCMC language similar to BUGS) is installed on your machine. See here for further details:

http://mcmc-jags.sourceforge.net/

**Note that there are re compatibility issues between matjags and JAGS 4.X To run the MATLAB code you will need to install JAGS 3.4.0 rather than the latest version.** The model files work fine with JAGS 4.X when called from R with rjags.

2) The main functions are fit_meta_d_mcmc (for fitting individual subject data) and fit_meta_d_mcmc_group (for hierarchical fits of group data). More information is contained in the help of these two functions and in the wiki https://github.com/smfleming/HMM/wiki/HMeta-d-tutorial. To get started try running exampleFit or exampleFit_group. 

A walkthrough of the model and intuitions behind different usages can be found in Olivia Faull's step-by-step tutorial developed for the Zurich Computational Psychiatry course: https://github.com/metacoglab/HMeta-d/blob/master/CPC_metacog_tutorial/cpc_metacog_tutorial.m

Please get in touch with your experiences with using the toolbox, and any bug reports or issues to me at stephen.fleming@ucl.ac.uk 

**License**

This code is being released with a permissive open-source license. You should feel free to use or adapt the utility code as long as you follow the terms of the license, which are enumerated below. If you use the toolbox in a publication we ask that you cite the following paper:

Fleming, S.M. (2017) HMeta-d: hierarchical Bayesian estimation of metacognitive efficiency from confidence ratings, Neuroscience of Consciousness, 3(1) nix007, https://doi.org/10.1093/nc/nix007

Copyright (c) 2017, Stephen Fleming

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
