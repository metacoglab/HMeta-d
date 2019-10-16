load dic
model in "/Users/sfleming/Dropbox/Utils/HMeta-d/Matlab/Bayes_metad_group_nodp.txt"
data in jagsdata.R
compile, nchains(1)
parameters in jagsinit1.R
initialize
update 1000
monitor set d1, thin(1)
monitor set c1, thin(1)
monitor set mu_logMratio, thin(1)
monitor set sigma_logMratio, thin(1)
monitor set mu_c2, thin(1)
monitor set sigma_c2, thin(1)
monitor set Mratio, thin(1)
monitor set cS1, thin(1)
monitor set cS2, thin(1)
monitor deviance
update 10000
coda *, stem('CODA1')
