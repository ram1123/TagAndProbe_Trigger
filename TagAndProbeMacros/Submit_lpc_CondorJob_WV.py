#! /usr/bin/env python
import os
import glob
import math
from array import array
import sys
import time
import subprocess
import tarfile
import datetime
import commands

OUTDIR = 'SingleElectron_Run2017B'
changes = raw_input("\n\nWrite change summary: ")

print "==> ",changes


currentDir = os.getcwd();
CMSSWDir =  currentDir+"/../";
print "PWD = ",currentDir

TestRun = 0

if OUTDIR == "":
	JobName = "job"
else:
	JobName = OUTDIR
# Get date and time for output directory
## ADD "test" IN OUTPUT FOLDER IF YOU ARE TESTING SO THAT LATER YOU REMEMBER TO WHICH DIRECTORY YOU HAVE TO REMOVE FROM EOS
os.system('xrdfs root://cmseos.fnal.gov/ mkdir /store/user/rasharma/aQGC_Ntuples/TriggerEfficiency/TriggerScaleFactors/' + OUTDIR)
if TestRun:
	outputFolder = "/store/user/rasharma/aQGC_Ntuples/TriggerEfficiency/TriggerScaleFactors/"+OUTDIR+'/'+datetime.datetime.now().strftime('%Y_%m_%d_%Hh%M')+"_TEST/";
	OutputLogPath = "Logs/"+OUTDIR+'/' + datetime.datetime.now().strftime('%Y_%m_%d_%Hh%M') + "_TEST";
else:
	outputFolder = "/store/user/rasharma/aQGC_Ntuples/TriggerEfficiency/TriggerScaleFactors/"+OUTDIR+'/'+datetime.datetime.now().strftime('%Y_%m_%d_%Hh%M');
	OutputLogPath = "Logs/"+OUTDIR+'/' + datetime.datetime.now().strftime('%Y_%m_%d_%Hh%M');


print "Name of output dir: ",outputFolder
print "Name of Log dir: ",OutputLogPath
# create a directory on eos
os.system('xrdfs root://cmseos.fnal.gov/ mkdir ' + outputFolder)
# create directory in pwd for log files
os.system('mkdir -p ' + OutputLogPath)

def exclude_function(filename):
    if filename.endswith('.root'):
            return True
    else:
            return False

## Function to create a tar file
def make_tarfile(output_filename, source_dir):
    with tarfile.open(output_filename, "w:gz") as tar:
            tar.add(source_dir, arcname=os.path.basename(source_dir),  exclude=exclude_function)

# Get CMSSW directory path and name
cmsswDirPath = commands.getstatusoutput('echo ${CMSSW_BASE}')
CMSSWRel = os.path.basename(cmsswDirPath[1])

print "CMSSW release used : ",CMSSWRel

# create tarball of present working CMSSW base directory
os.system('rm CMSSW*.tgz')
make_tarfile(CMSSWRel+".tgz", cmsswDirPath[1])

# send the created tarball to eos
os.system('xrdcp -f ' + CMSSWRel+".tgz" + ' root://cmseos.fnal.gov/'+outputFolder+"/" + CMSSWRel+".tgz")
os.system('git diff > main.patch')
os.system('git diff > mypatch.patch')
os.system("sed -i '1s/^/Changes Summay : "+changes+"\\n/' mypatch.patch")
os.system('git log -1 --format="%H" >> mypatch.patch ')
os.system('xrdcp -f mypatch.patch root://cmseos.fnal.gov/'+outputFolder+'/mypatch.patch')
os.system('xrdcp -f main.patch root://cmseos.fnal.gov/'+outputFolder+'/main.patch')



outJDL = open("runstep2condor_WV.jdl","w");


outJDL.write("Executable = runstep2condor_WV.sh\n");
outJDL.write("Universe = vanilla\n");
#outJDL.write("Requirements =FileSystemDomain==\"fnal.gov\" && Arch==\"X86_64\"");
outJDL.write("Notification = ERROR\n");
outJDL.write("Should_Transfer_Files = YES\n");
outJDL.write("WhenToTransferOutput = ON_EXIT\n");
#outJDL.write("include : list-infiles.sh |\n");
#outJDL.write("Transfer_Input_Files = "+inputlist+"\n");
outJDL.write("x509userproxy = $ENV(X509_USER_PROXY)\n");

outJDL.write("Output = "+OutputLogPath+"/"+JobName+".stdout\n");
outJDL.write("Error  = "+OutputLogPath+"/"+JobName+".stdout\n");
outJDL.write("Log = "+OutputLogPath+"/"+JobName+".log\n");
outJDL.write("Queue\n");
	    
outJDL.close();

command1 = "g++ Run_TnP.cxx -o  tnp_Ele -std=c++0x `root-config --libs --cflags` -include  TagAndProbe_Ele.C"
command2 = "./tnp_Ele RootFiles_Run2017B.txt" 

outScript = open("runstep2condor_WV.sh","w");
outScript.write('#!/bin/bash');
outScript.write("\n"+'source /uscmst1/prod/sw/cms/bashrc prod');
outScript.write("\n"+'echo "Starting job on " `date`');
outScript.write("\n"+'echo "Running on: `uname -a`"');
outScript.write("\n"+'echo "System software: `cat /etc/redhat-release`"');
outScript.write("\n"+'source /cvmfs/cms.cern.ch/cmsset_default.sh');
outScript.write("\n"+'### copy the input root files if they are needed only if you require local reading');
outScript.write("\n"+'xrdcp -s root://cmseos.fnal.gov/'+outputFolder+"/" + CMSSWRel+".tgz  .");
outScript.write("\n"+'tar -xf '+ CMSSWRel +'.tgz' );
outScript.write("\n"+'rm '+ CMSSWRel +'.tgz' );
outScript.write("\n"+'cd ' + CMSSWRel + '/src/TagAndProbe_Trigger/TagAndProbeMacros/' );
outScript.write("\n"+'echo "============================================" ');
outScript.write("\n"+'rm *.root');
outScript.write("\n"+'echo "====> List output files : " ');
outScript.write("\n"+'ls');
outScript.write("\n"+'echo "============================================" ');
outScript.write("\n"+'echo "pwd = $PWD"');
outScript.write("\n"+'echo "===========\t Project Rename 	=====================" ');
outScript.write("\n"+'scramv1 b ProjectRename');
outScript.write("\n"+'echo "===========\t Load CMS environment	=====================" ');
outScript.write("\n"+'eval `scram runtime -sh`');
outScript.write("\n"+'echo "===========\t scram clean and scram b ===================" ');
outScript.write("\n"+'scramv1 b clean; scramv1 b');
outScript.write("\n"+'echo "============================================" ');
outScript.write("\n"+'echo "====> List output files : " ');
outScript.write("\n"+'ls');
outScript.write("\n"+'echo "===========\t Run main script	=====================" ');
outScript.write("\n"+command1);
outScript.write("\n"+command2);
outScript.write("\n"+'echo "====> List output files : " ');
outScript.write("\n"+'ls -ltrh');
outScript.write("\n"+'echo "xrdcp output for condor"');
outScript.write("\n"+'for FILE in *.root');
outScript.write("\n"+'do');
outScript.write("\n"+'\techo "xrdcp -r -f ${FILE} root://cmseos.fnal.gov/'+outputFolder+'/"');
outScript.write("\n"+'\txrdcp -r -f ${FILE} root://cmseos.fnal.gov/'+outputFolder+'/ 2>&1');
outScript.write("\n"+'done');
outScript.write("\n"+'cd ${_CONDOR_SCRATCH_DIR}');
outScript.write("\n"+'rm -rf ' + CMSSWRel);
outScript.write("\n");
outScript.close();
os.system("chmod 777 runstep2condor_WV.sh");

print "===> Set Proxy Using:";
print "\tvoms-proxy-init --voms cms --valid 168:00";
print "\"condor_submit runstep2condor_WV.jdl\" to submit";
os.system("condor_submit runstep2condor_WV.jdl")
