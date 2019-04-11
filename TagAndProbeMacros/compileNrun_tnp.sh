#!/bin/bash

CMACRO=$1  # TagAndProbe_Ele.C or TagAndProbe_Mu.C
EXEC="tnp"

if [[ "$#" -lt 1 ]];
then
    echo " "
    echo "====>>  no argument is passed. Pass either TagAndProbe_Ele.C OR TagAndProbe_Mu.C"
    echo " "
    exit
fi


if [[ "$CMACRO" =~ "Ele" ]]; then
    echo "compiling for electron : $CMACRO" 
    EXEC="tnp_Ele"
fi
if [[ "$CMACRO" =~ "Mu" ]]; then
    echo " "
    echo "compiling for muon : $CMACRO"
    echo " "
    EXEC="tnp_Mu"
fi
echo "ram: "$EXEC
echo "ram: "$CMACRO
g++ Run_TnP.cxx -o $EXEC  -std=c++0x `root-config --libs --cflags` -include $CMACRO


############  Enter the name of text file in which root files are entered

echo " "
echo "1) Run Interactively : Enter : 0"
echo "2) Run in lxplus batch : Enter : 1"
echo " "

#read run
run=0
if [[ $run == 0 ]]; then
ls -ltrh *.txt
echo "You are running in interactive mode!!!"
echo "Enter the name of text file having all the root files..."
#read textFile
#textFile="RootFiles_EOS.txt"
textFile="RootFiles.txt"
FirstFile=`head -1 $textFile`
#if [[ "$FirstFile" =~ "root://" ]];   # If your files are present not at cern then access files through xrootd services. That needs to generate a valid proxy.
#then
#        echo " "
#        echo "==== need to generate proxy, since it is going to use xrootd services"
#        echo "==== running voms-proxy-init --voms cms first ===="
#        echo "==== Enter password to generate proxy ===="
#        echo " "
#        voms-proxy-init --voms cms --valid 168:00
#       echo "-------------------------------"
#       echo "---"
#       echo "---   now running ./$EXEC $textFile  "
#       echo "---"
#       echo "-------------------------------"
#       ./$EXEC $textFile
#else
#       echo "-------------------------------"
#       echo "---"
#       echo "---   now running ./$EXEC $textFile  "
#       echo "---"
#       echo "-------------------------------"
#     ./$EXEC $textFile
#  fi
fi

if [[ $run == 1 ]]; then 
QUEUE="1nd"
ls -ltrh *.txt
echo "You are going to submit jobs on lxplus batch!!!"
echo "Enter the name of text file having all the root files"
read textFile
echo "-------------------------------"
echo "---"
echo "---   now running ./$EXEC TEXT_File_WITH_ALL_ROOT_FILES  "
echo "---"
echo "-------------------------------"
FirstFile=`head -1 $textFile`
echo $FirstFile
if [[ "$FirstFile" =~ "root://" ]];   # If your files are present not at cern then access files through xrootd services. That needs to generate a valid proxy.
then
    if [[ `ls` =~ "x509up" ]]; then  
        TimeLeft=`voms-proxy-info | cut -d ":"  -f 2 | sed 's/ Digital .*//' | tail -c 5`
        echo $TimeLeft
        if [[ $TimeLeft -lt 10 ]]; then
        echo "Having a vaild proxy with less than 10 hours left. Better generate a new proxy again."
        echo "==== running voms-proxy-init --voms cms first ===="
        echo "==== Enter password to generate proxy ===="
        voms-proxy-init --voms cms --valid 168:00
        cp `voms-proxy-info --path` .  ### if ntuple is not at cern eos then proxy is needed
        fi
    else
        echo "==== running voms-proxy-init --voms cms first ===="
        echo "==== Enter password to generate proxy ===="
        voms-proxy-init --voms cms --valid 168:00
        echo "copying the generated proxy file from /tmp to $PWD"
        cp `voms-proxy-info --path` .  ### if ntuple is not at cern eos then proxy is needed
    fi
     echo "Enter the name of File which you want to submit on lxplus batch queue : "
     read batchFileName
     if [[ `ls` =~ $batchFileName ]]; then
         echo "Deleting the existing $batchFileName  and creating a fresh one"
         rm $batchFileName 
     fi
     touch $batchFileName
     echo "#!/bin/bash" >> $batchFileName
     echo "cd $CMSSW_BASE/src/TagAndProbe_Trigger/TagAndProbeMacros" >> $batchFileName
     echo "export X509_USER_PROXY=\`ls x509up*\`" >> $batchFileName
     echo "./$EXEC $textFile" >> $batchFileName
     bsub -q $QUEUE < $batchFileName
else
    echo "All files listed in text file either present in local or at cern itself"
    echo "Do you want to run interactively then press y or Y else job will be submit to batch queues"
    #read Decision
		Decision=Y
    if [[ $Decision == y ]] || [[ $Decision == Y ]]; then
        ./$EXEC $textFile
    else
     echo "Enter the name of File which you want to submit on lxplus batch queue : "
     read batchFileName
     if [[ `ls` =~ $batchFileName ]]; then
         echo "Deleting the existing $batchFileName  and creating a fresh one"
         rm $batchFileName 
     fi  
     touch $batchFileName
     echo "#!/bin/bash" >> $batchFileName
     echo "cd $CMSSW_BASE/src/TagAndProbe_Trigger/TagAndProbeMacros" >> $batchFileName
     echo "./$EXEC $textFile" >> $batchFileName
		 ########################################################
		 ########################################################
     bsub -q $QUEUE < $batchFileName
    fi

 fi  ## check name of root files

fi   ## Run = 1







