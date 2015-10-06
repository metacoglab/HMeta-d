load dic
model in "/Users/sfleming/Documents/HMM/Bayes_metad.txt"
data in jagsdata.R
compile, nchains(1)
parameters in jagsinit3.R
initialize
update 1000
monitor set meta_d, thin(1)
monitor set cS1, thin(1)
monitor set cS2, thin(1)
monitor deviance
update 10000
coda *, stem('CODA3')
