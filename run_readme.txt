The following are instructions on how to run the code on an el7 cms workgroup 
server on lxplus or at DESY.
To set up and compile, follow instructions on 

compile_4_2_8.readme.txt   for  2010 0.9 and 7 TeV pp data and MC, and
                                2011 2.76 TeV MC
compile_4_4_7.readme.txt   for  2011 2.76 TeV pp data
compile_5_3_32.readme.txt  for  2011 7 TeV data and MC, and 
                                2012 8 TeV data and MC
compile_7_6_7.readme.txt   for  2015 13 TeV AOD or miniAOD data and MC, on slc6
compile_7_6_4.readme.txt   for  2015 13 TeV AOD or miniAOD data and MC, on el7
compile_10_6_4.readme.txt  for  2016 13 TeV UL AOD data and MC (not yet public),
                                for cross-validation with nanoAODv8

To run, either on the el7 shell, or on the singularity shell, whichever 
works (depending on the setup sometimes one will work, sometimes the other, 
and sometimes both):

   At DESY, if not yet done, load the cmssw environment (automatic on lxplus):
   module use -a /afs/desy.de/group/cms/modulefiles
module load cmssw

   set the cmssw environment:
cmsenv

   run the configuration: (needs to match compilation for AOD or miniAOD)
cmsRun jobname.py

The following example configurations are provided so far:
nanoanalyzer_cfg_2010Data.py  for 2010 7 TeV pp data (similar for 0.9 TeV)
nanoanalyzer_cfg_2010MC.py    for 2010 7 TeV pp MC (similar for 0.9/2.76 TeV)
nanoanalyzer_cfg_2011Datapp2.76.py  for 2011 2.76 TeV pp data 
nanoanalyzer_cfg_2011Data.py  for 2011 7 TeV pp data (similar for 2012 8 TeV)
nanoanalyzer_cfg_2011MC.py    for 2011 7 TeV pp MC (similar for 2012 8 TeV)
nanoanalyzer_cfg_2015Data.py  for 2015 13 TeV pp data
nanoanalyzer_cfg_2015miniMC.py  for 2015 13 TeV pp miniAOD MC

The configurations for the 2015 data have been technically validated,
but not yet been validated for content 

Run 2 UL data have not yet been released.
Official CMS nanoAOD will be the main reference for these. 
