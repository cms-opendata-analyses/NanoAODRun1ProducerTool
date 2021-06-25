##################################################################################
# Nanoanalyzer configuration file for all years                                  #
# Use for HT Condor and VM                                                       #
# Uncomment & comment relevant lines before you run it                           #
##################################################################################

import sys
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.Types as CfgTypes
import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.Utilities.FileUtils as FileUtils
from RecoMuon.TrackingTools.MuonServiceProxy_cff import *
from FWCore.ParameterSet.VarParsing import VarParsing

process = cms.Process("Nano")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
#process.load("Configuration.Geometry.GeometryIdeal_cff") # for 2011
process.load("Configuration.StandardSequences.Geometry_cff") # for 2010
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.load("Configuration.Geometry.GeometryRecoDB_cff") # for Run2
#process.load('Configuration.StandardSequences.EndOfProcess_cff') # for Run2

################################### VARPASSING ###################################
# Use VarParsing to specify your input directory and output file
# Comment Varparsing if you want to submit job using CRAB
# because CRAB does not support input directory as we put input dataset directly
#options = VarParsing('analysis')
#options.register(
#    "inputDir",
#    "",
#    VarParsing.multiplicity.singleton,
#    VarParsing.varType.string,
#    "Input directory with inputs"
#    )
#options.register(
#    "outputName",
#    "",
#    VarParsing.multiplicity.singleton,
#    VarParsing.varType.string,
#    "Output name"
#    )
#options.parseArguments()

#if (options.inputDir == ""):
#    sys.exit("Directory to find input file where???")
#else:
    # 2010 VM
#    InputDir = "/home/cms-opendata/CMSSW_4_2_8/src/NanoAOD/NanoAnalyzer/" + options.inputDir
    # 2011 NAF/HTC
    # InputDir = "/nfs/dust/cms/user/zulaiha/PhD/CMSSW_5_3_32/src/NanoAOD/NanoAnalyzer/" + options.inputDir
    # 2015 NAF/HTC
    #InputDir = "/nfs/dust/cms/user/zulaiha/PhD/CMSSW_7_6_1/src/NanoAOD/NanoAnalyzer/" + options.inputDir

##################################################################################

###################################### GLOBAL TAG ################################
# Change the global tag accordingly

# "standard" 2010 Data (recheck for particular dataset used)
#process.GlobalTag.globaltag = 'FT_R_42_V10A::All'
# "standard" 2010 MC (recheck for particular MC set used)  
process.GlobalTag.globaltag = 'START42_V17B::All'

# "standard" 2011 Data (recheck for particular dataset used)
#process.GlobalTag.globaltag = 'FT_53_LV5_AN1::All'
# "standard" 2011 MC (recheck for particular dataset used)
#process.GlobalTag.globaltag = 'START53_LV6::All'

##################################################################################

# Intialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('Nano')
# steers how often intermediate output is given
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
# to get summary
process.options   = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
# To avoid crash when some data base input not found (for debug only)
#process.options   = cms.untracked.PSet(SkipEvent = cms.untracked.vstring('ProductNotFound'))
                                        
# Set the maximum number of events to be processed here (-1=all)
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))

#####################################  JSON FILE #################################
# Change the directory and JSON event filter file accordingly
# Only uncomment if you run in Data

# "standard" 2010 Golden JSON (recheck for your particular dataset)
#goodJSON = 'files/Cert_136033-149442_7TeV_Apr21ReReco_Collisions10_JSON_v2.txt'

# "standard" 2011 Golden JSON (recheck for your particular dataset)
#goodJSON = 'files/Cert_160404-180252_7TeV_ReRecoNov08_Collisions11_JSON.txt'

# "standard" 2012 Golden JSON (recheck for your particular dataset)
#goodJSON = 'files/Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt'

##################################################################################

# Get the luminosity list from the JSON
# Only uncomment if you run in Data
#myLumis = LumiList.LumiList(filename = goodJSON).getCMSSWString().split(',')

##################################################################################

# Load jet correction services for all jet algoritms
process.load("JetMETCorrections.Configuration.JetCorrectionServicesAllAlgos_cff")
#

#################################### INPUT FILE ##################################
# To test locally or submit batch job using condor, use this:
#fileinPut = FileUtils.loadListFromFile (InputDir)

#process.source = cms.Source("PoolSource",
#                            fileNames = cms.untracked.vstring(*fileinPut)
#)

# To submit batch job using CRAB or test locally (2nd option), use this:
process.source = cms.Source("PoolSource",
                            # 2010 Data (test on one file)
                            #fileNames = cms.untracked.vstring('root://eospublic.cern.ch//eos/opendata/cms/Run2010B/Mu/AOD/Apr21ReReco-v1/0000/00459D48-EB70-E011-AF09-90E6BA19A252.root')
                            # 2010 MC (test on one file)
                            fileNames = cms.untracked.vstring('root://eospublic.cern.ch//eos/opendata/cms/MonteCarlo2010/Summer12/DYToMuMu_M-20_TuneZ2Star_HFshowerLibrary_7TeV_pythia6/AODSIM/LowPU2010_DR42_PU_S0_START42_V17B-v1/10000/0A78309F-AF5B-E211-B615-003048FFCB8C.root')
                            # 2011 Data (test on one file)
                            #fileNames = cms.untracked.vstring('root://eospublic.cern.ch//eos/opendata/cms/Run2011A/SingleMu/AOD/12Oct2013-v1/10000/00209631-9C37-E311-BACA-002590494FDE.root')
                            # 2011 MC (test on one file)
                            #fileNames = cms.untracked.vstring('root://eospublic.cern.ch//eos/opendata/cms/MonteCarlo2011/Summer11LegDR/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/AODSIM/PU_S13_START53_LV6-v1/00000/0019AB30-B9B7-E311-9E28-003048FF86CA.root')
)


##################################################################################

# Process the lumi
# Only uncomment if you run on Data
#process.source.lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange())
#process.source.lumisToProcess.extend(myLumis)

# *************************************************
# number of events to be skipped (0 by default)   *
# *************************************************
process.source.skipEvents = cms.untracked.uint32(0)

#example for parton filter on MC
#process.partonfilter = cms.EDFilter("PythiaFilter",
#                                    ParticleID = cms.untracked.int32(4)
#)
#process.charmfilter = cms.Path(process.partonfilter)

# Process the analyzer 
# ************************************************************
# define your output path and file here (default is test.root)
# ************************************************************
process.nano = cms.EDAnalyzer('NanoAnalyzer',
                              outFile = cms.string('test.root'),

                              # Change this:
                              # If MC:
                              isData = cms.bool(False)
                              # If Data:
                              #isData = cms.bool(True)
                          )
#process.p = cms.EndPath(process.nano)
#process.schedule = cms.Schedule(process.charmfilter, process.p)

process.p = cms.Path(process.nano)
