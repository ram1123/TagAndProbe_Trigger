# TagAndProbe_Trigger

###### https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2#Working_points_for_92X_samples_R

## Package Setup

```bash
cmsrel CMSSW_9_4_13 
cd CMSSW_9_4_13/src  
cmsenv  
git clone git@github.com:ram1123/TagAndProbe_Trigger.git
scram b -j 4  
```

## For Test Run 

```bash
cd $CMSSW_BASE/src/TagAndProbe_Trigger/NtupleProducer/test   
cmsRun runNtupler.py  
```

## To Submit Crab Job

```bash
# set crab environment
source /cvmfs/cms.cern.ch/crab3/crab.csh
# Submit the job
python multicrab_Ele.py
```

## Get Efficiency

```bash
cp TnP_ntuple.root $CMSSW_BASE/src/TagAndProbe_Trigger/TagAndProbeMacros  (one example file is already present in the directory)  
cd $CMSSW_BASE/src/TagAndProbe_Trigger/TagAndProbeMacros/  
./compileNrun_tnp.sh TagAndProbe_Ele.C 
#(OR TagAndProbe_Mu.C, This script is self explaining when run it)   

./runPlotting.sh histNames_Ele.txt 
#(OR histNames_Mu.txt, this will generate final results in the "results" directory)   
```

### Condor job submission

```bash
voms-proxy-init --voms cms --valid 168:00
python Submit_lpc_CondorJob_WV.py
```
