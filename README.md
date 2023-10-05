# Instruction to set up env
## Log into lxplus
```
ssh -Y @lxplus.cern.ch -L localhost:XXXX:localhost:XXXX
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

## Open jupyter notebook
```
jupyter notebook --no-browser --port=XXXX  
```