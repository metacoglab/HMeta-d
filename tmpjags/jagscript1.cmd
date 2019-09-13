load dic
model in "/Users/sfleming/Dropbox/Utils/HMeta-d/Bayes_metad.txt"
data in jagsdata.R
compile, nchains(1)
parameters in jagsinit1.R
initialize
update 1000
monitor set meta_d, thin(1)
monitor set d1, thin(1)
monitor set c1, thin(1)
monitor set cS1, thin(1)
monitor set cS2, thin(1)
monitor deviance
update 10000
coda *, stem('CODA1')
