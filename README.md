NanoAODRun1 Producer Tool
----------------------------------

- [Setting up in a CMSSW Docker container (Open Data Users)](#docker)

- [Setting up with CERN credentials](#cern)

### <a name="docker">Setting up in a CMSSW Docker container (Open Data Users)</a>

The following are instructions on how to set up and compile the code on a CMSSW docker container, with `CMSSW_X_X_X`. (Select the container version according to the recommendation for each year.)

Start a docker container following the [instruction](https://opendata.cern.ch/docs/cms-guide-docker).

##### WITHIN the docker container

Once (first time only):

Clone the GitHub repository to `NanoAOD/NanoAnalyzer/`:
```
git clone https://github.com/cms-opendata-analyses/NanoAODRun1ProducerTool.git NanoAOD/NanoAnalyzer/
```

`NanoAOD/NanoAnalyzer/` is the name of the folder. The code won’t compile if it is not cloned to this folder.

If you are already in the `CMSSW_X_X_X/src` directory,

Every time:
```
cmsenv
```

If not,
```
cd CMSSW_X_X_X/src
cmsenv
```

Every time
```
cd NanoAOD/NanoAnalyzer
```

Here you will now find the content of the repository.

Technically, the NanoAnalyzer is a so-called EDAnalyzer
(originally set up with mkdanlzr).

To compile:
```
scram b
```

For a successful compilation, you should see something like:
```
>> Local Products Rules ..... started
>> Local Products Rules ..... done
>> Building CMSSW version CMSSW_4_2_8 ----
>> Entering Package NanoAnalyzer/src
>> Leaving Package NanoAnalyzer/src
>> Package NanoAnalyzer/src built
>> Entering Package NanoAnalyzer/pytest
>> Leaving Package NanoAnalyzer/pytest
>> Package NanoAnalyzer/pytest built
>> Entering Package NanoAnalyzer/tools
>> Leaving Package NanoAnalyzer/tools
>> Package NanoAnalyzer/tools built
>> Entering Package NanoAnalyzer/files
>> Leaving Package NanoAnalyzer/files
>> Package NanoAnalyzer/files built
>> Subsystem NanoAnalyzer built
>> Local Products Rules ..... started
>> Local Products Rules ..... done
gmake[1]: Entering directory `/home/cmsusr/CMSSW_4_2_8'
>> All python modules compiled
>> Pluging of all type refreshed.
gmake[1]: Leaving directory `/home/cmsusr/CMSSW_4_2_8'
```

#### Whenever you want to update the code to the latest version (not first time)

In the `NanoAOD/NanoAnalyzer/` directory, 

get the updates from github
```
git pull
```

#### Run NanoAODAnalyzer
	
Example configuration files for different years of data and simulation are found in the files called `nanoanalyzer_cfg_20*.py`. Choose an appropriate example for your container and data format. Run it using: 
```
cmsRun nanoanalyzer_cfg_<year><Data or MC>.py
```

e.g. `nanoanalyzer_cfg_2010Data.py` for 2010 7 TeV pp AOD data (similar for 0.9 TeV) or `nanoanalyzer_cfg_2015miniMC.py` for 2015 13 TeV pp miniAOD MC.

Note:

The configurations for the 2015 data have been technically validated, but not yet been validated for content. 

Run2 UL data have not yet been released.

Official CMS NanoAOD will be the main reference for these.

<br>

### <a name="cern">Setting up with CERN credentials</a>

For instructions how to load and compile the code with 
   CMSSW 5_3_32 (for 2011 7 TeV and 2012 8 TeV pp data), see
compile_5_3_32.readme.txt

For instructions how to load and compile the code with 
   CMSSW 4_2_8 (for 2010 7 and 0.9 TeV pp data, and 2011 2.79 TeV pp data), see
compile_4_2_8.readme.txt

For instructions how to load and compile the code with 
   CMSSW 7_6_4/el7 (for 2015 5 and 13 TeV pp AOD data), see
compile_7_6_4.readme.txt

For instructions how to load and compile the code with 
   CMSSW 7_6_7/slc6 (for 2015 5 and 13 TeV pp AOD data), see
compile_7_6_7.readme.txt

For instructions how to load and compile the code with 
   CMSSW 10_6_4 (for validation against 2017 and 2018 UL AOD data), see
compile_10_6_4.readme.txt

For instructions how to run in various configurations, see
run_readme.txt
