The ntuple making part works for 2021 Winter20 samples in this branch.
For 2018 we need to use CMSSW_10_2_X and for Summer19 one should use EGMHLT_run3 branch
Although the macros works for all years.

# Code setup

```bash
cmsrel CMSSW_11_2_4
cd CMSSW_11_2_4/src
cmsenv
git clone git@github.com:arunhep/TagAndProbe_Trigger.git
git checkout EGMHLT_run3_110X
scramv1 b
```

# For Test Run (check the input file and triggers)

```bash
cd $CMSSW_BASE/src/TagAndProbe_Trigger/NtupleProducer/test
cmsRun runNtupler.py
```

# For crab job submission (check the paths etc in the file)

```bash
voms-proxy-init --voms cms --valid 168:00
source /cvmfs/cms.cern.ch/common/crab-setup.sh
python multicrab_Ele.py
```

### Macros to produce the plots

```bash
cd TagAndProbe_Trigger/TagAndProbeMacros
./run.sh
```
