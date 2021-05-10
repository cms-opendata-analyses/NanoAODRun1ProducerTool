##############################################
Documentation on edm dump event content
Collection of input tag in AOD/miniAOD
Some description: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideRecoDataTable
##############################################

How to list the available collections from an EDM rootfile?
===========================================================

- Go to any latest CMSSW version
- Do cmsenv
- Do edmDumpEventContent --all root://xrootd-cms.infn.it/<one_rootfile> > <any_textfile>
  e.g.: edmDumpEventContent --all root://xrootd-cms.infn.it//store/data/Run2016G/DoubleMuon/AOD/23Sep2016-v1/80000/5AF8FEE3-3C8C-E611-818C-008CFAF0842A.root > evtcontent_AOD.txt


How to use the collection?
===========================
In the collection, there are column type, module, label, process, full name
           
E.g.: 
    Type                                   module                 label             process             full name
vector<reco::Track>                   "generalTracks"             ""                "RECO"         recoTracks_generalTracks__RECO

