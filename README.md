# Instruction to set up env
## Log into lxplus
```
ssh -Y <username>@lxplus.cern.ch -L localhost:XXXX:localhost:XXXX
```
or if using desy remote server, do
```
ssh <username>@@naf-cms.desy.de -L localhost:XXXX:localhost:XXXX
```


## Set up CMSSW environment
```
source /cvmfs/cms.cern.ch/cmsset_default.sh  
scram p cmssw CMSSW_13_2_5  
cd CMSSW_13_2_5/src  
cmsenv  
```

## Clone the repository
```
git clone https://:@gitlab.cern.ch:8443/cms-podas23/dpg/hcal-das.git
```
or from github if you don't have gitlab permission
```
git clone https://github.com/lwang046/PODAS-HCAL.git
```


## Open jupyter notebook
```
jupyter notebook --no-browser --port=XXXX  
```
