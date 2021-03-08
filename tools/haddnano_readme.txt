############################################################
#### Documentation on haddnano.py usage for nanoAODplus ####
############ Nuha for nanoAODplus team #####################
############################################################

=========
Problem?
=========
Cannot merge nanoAODplus Ntuples that has different triggers in different files into one file using standard hadd command (see man hadd)

===================
Proposed solution?
===================
using haddnano.py (tools recommended by CMS nanoAOD team)

=====
How?
=====
One example with CMSSW dependent:
Tested on naf-cms20, SLC6
Procedure:
1. go to cmssw (I tried on CMSSW_10_2_14, which used ROOT version 6.12/07), and setup the environment (cmsenv)
2. go to haddnano.py directory and do: python haddnano.py <output.root> <directory_of_input_file*.root>

---------- printout message ------------------------
Adding file /pnfs/desy.de/cms/tier2/store/user/nujomhar/Data2010/CRAB_UserFiles/Data10MuOniaRunAKKbeauty/200803_133802/0000/Data10_Mu_KKbeauty_13.root
Merging Events
missing: ['HLT_EcalOnly_SumEt160', 'HLT_DoublePhoton15_L1R', 'HLT_SingleLooseIsoTau20_Trk5', 'HLT_SingleLooseIsoTau25_Trk5', 'HLT_Ele15_SW_EleId_L1R', 'HLT_DoubleEle10_SW_L1R', 'HLT_Jet70U', 'HLT_PixelTracks_Multiplicity40', 'HLT_Ele15_SW_L1R', 'HLT_Mu0_TkMu0_Jpsi_NoCharge', 'HLT_MinBiasBSC'] 
 Additional: []
Merging h_trackpt
Merging h_trackptlow
Merging h_tracketa
Merging h_d0pt
Merging h_dstarpt
--------- printout message end -----------------------

======================================
Can haddnano.py be cmssw independent?
======================================

YES!! At the end, this is what we want! (Assuming your input root files are not corrupted!)

1. Tested in own laptop (direct access to the input files e.g in hardisk etc., however not practical)
        - Works fine with ROOT6 version 6.14/04 installed in laptop

2. Tested on several naf server and several ROOT version WITHOUT cmssw setup environment (cmsenv) 
        - Does not work on several root versions in naf unfortunately. Check ROOT version with this command: $root --version
        - Tested on several naf server in SLC6: ssh -X zulaiha@naf-desy-cms.desy.de 
             - by default when login naf: ROOT version 5.34/38 (not work, gave output file but *** Break *** segmentation violation at the end. Be careful because it produced corrupted output file!!)
             - ROOT5 version 5.34/23 (not work, gave this error: Traceback (most recent call last):........, no output is produced)
             - ROOT6 version 6.02/05 (not work, gave same error as above: Traceback (most recent call last):........, no output is produced)
             - ROOT6 version 6.22/00 (WORK! simply source this latest ROOT version for SLC6 in naf: $source /cvmfs/sft.cern.ch/lcg/views/LCG_98/x86_64-slc6-gcc8-opt/setup.sh, produced output)
        - Tested also for el7: ssh -X zulaiha@naf-desy-cms-el7.desy.de
             - ROOT6 version 6.22/00 (WORK! simply source this latest ROOT version for el7 in naf: $source /cvmfs/sft.cern.ch/lcg/views/LCG_98/x86_64-centos7-gcc8-opt/setup.sh, produced output)
             
   ***Notes*** 
   - The best is to use root that is compatible with the machine that you are using, e.g slc-gcc8 root for SLC6 machine and centos7-gcc8 for el7 machine (see above). 
   - This is to avoid complication (eg. takes long time to open TBrowser when use root compatible with SLC6 but in el7 machine)
   - Files and the command that is used for these testing: python haddnano.py testing.root /pnfs/desy.de/cms/tier2/store/user/nujomhar/Data2010/CRAB_UserFiles/Data10MuRunAKKbeauty/200803_190836/0000/Data10_Mu_KKbeauty_11*.root     

==================================
Run haddnano.py in naf background
==================================
nohup sh -c "date; ./haddnano.py <outputfile.root> <dir+inputfile>*.root; date" &> <printout_file>.txt &

where:
        sh = shell. Theoretically zsh and bash can use this command but need to see.

        -c = string

        "<commands>" = all the commands you want to execute, separated by semicolon. In this case, 3 commands will be execute; 
        1. date+time when you submit the job, 2. run your codes (in this case run haddnano.py), 3. date+time after your job finish

        &> <printout_file>.txt = dump all the information (the dates, cout (if you have in your .cc code), etc) in the text file

        & = run in background

e.g
nohup sh -c "date; ./haddnano.py Data10_MuMonitorRB_KKbeauty.root /pnfs/desy.de/cms/tier2/store/user/nujomhar/Data2010/CRAB_UserFiles/Data10ERunBKKbeauty/200805_191353/0000/*.root; date" &> test.txt &

For this example, it takes ~1 hour to merge 999 files

   ***Note***
   - Please use condor to submit huge nanoAOD files and NOT run naf in background. For 2010, it is still affordable!

Any questions, comments and suggestions, please let me know! Thanks!


