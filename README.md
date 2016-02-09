HMM
===

Hierarchical meta-d' model
***this is a beta version under active development***

Steve Fleming
stephen.fleming@ucl.ac.uk 

This toolbox implements the meta-d’ model (Maniscalco & Lau, 2012) in a hierarchical Bayesian framework using Matlab and JAGS, a program for conducting MCMC inference on arbitrary Bayesian models. A paper on the advantages of estimating meta-d’ in a hierarchal Bayesian framework is in preparation. For a more general introduction to Bayesian models of cognition see Lee & Wagenmakers, Bayesian Cognitive Modeling: A Practical Course http://bayesmodels.com/

The code is designed to work “out of the box” without much coding on the part of the user, and it receives data in the same format as Maniscalco & Lau’s toolbox, allowing easy switching and comparison between the two.

1) To get started, you need to first ensure JAGS (an MCMC language similar to BUGS) is installed on your machine. See here for further details:

http://mcmc-jags.sourceforge.net/

2) The main functions are fit_meta_d_mcmc (for fitting individual subject data) and fit_meta_d_mcmc_group (for hierarchical fits of group data). More information is contained in the help of these two functions. To get started try running exampleFit or exampleFit_group.

Please get in touch with your experiences with using the toolbox, and any bug reports or issues to me at stephen.fleming@ucl.ac.uk 
