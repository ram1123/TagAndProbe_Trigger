### compile the code will generate executable "tnp"

g++ Run_TnP.cxx -o tnp  -std=c++0x `root-config --libs --cflags` -include TagAndProbe_Ele.C or TagAndProbe_Mu.C

## To dump all root file name from eos use below command

```bash
eosls /store/user/<PATH_of_Ntuples> | awk '{print "root://cmseos.fnal.gov//store/user/<PATH_of_Ntuples>/"$1}' > RootFiles_EOS.txt
```

for example:
```bash
eosls /store/user/rasharma/aQGC_Ntuples/TriggerEfficiency/ElectronNtuple_HWW_2017_fixed_v2_92X/SingleElectron/Run2017B/190409_221045/0000/ | awk '{print "root://cmseos.fnal.gov//store/user/rasharma/aQGC_Ntuples/TriggerEfficiency/ElectronNtuple_HWW_2017_fixed_v2_92X/SingleElectron/Run2017B/190409_221045/0000/"$1}' > RootFiles_EOS.txt
```

