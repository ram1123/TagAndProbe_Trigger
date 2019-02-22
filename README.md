# TagAndProbe_Trigger

###### https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2#Working_points_for_92X_samples_R
cmsrel CMSSW_10_2_11  
cd CMSSW_10_2_11/src  
cmsenv  
git cms-init  
git cms-merge-topic cms-egamma:EgammaID_1023 ##if you want the V2 IDs, otherwise skip
git clone git@github.com:arunhep/TagAndProbe_Trigger.git
scram b -j8  

## For Test Run 
cd $CMSSW_BASE/src/TagAndProbe_Trigger/NtupleProducer/test   
<br>  
cmsRun runNtupler.py  
<br>  
cp TnP_ntuple.root $CMSSW_BASE/src/TagAndProbe_Trigger/TagAndProbeMacros  (one example file is already present in the directory)  
<br>  
cd $CMSSW_BASE/src/TagAndProbe_Trigger/TagAndProbeMacros/  
<br>  
./compileNrun_tnp.sh TagAndProbe_Ele.C 
(OR TagAndProbe_Mu.C, This script is self explaining when run it)   
<br>   
./runPlotting.sh histNames_Ele.txt 
(OR histNames_Mu.txt, this will generate final results in the "results" directory)   
