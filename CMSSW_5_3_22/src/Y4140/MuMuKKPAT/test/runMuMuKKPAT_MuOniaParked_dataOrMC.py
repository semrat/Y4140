import FWCore.ParameterSet.Config as cms
process = cms.Process('NTUPLE')
process.options = cms.untracked.PSet(
        wantSummary = cms.untracked.bool(True)
)



# Import of standard configurations
process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.suppressInfo = cms.untracked.vstring( "mkcands" )
process.MessageLogger.suppressWarning = cms.untracked.vstring( "mkcands" )
process.MessageLogger.cerr.FwkReport.reportEvery = 1



MC = False #True
if MC :
        official = True #official = False
MCMotherId = 511 # 511 B0, 531 Bs0
if MCMotherId == 511 :
    MCExclusiveDecay = True	
elif MCMotherId == 531 :
    MCExclusiveDecay = False		



# Input source
process.source = cms.Source("PoolSource",
                            skipEvents = cms.untracked.uint32( 0 ),
                            fileNames = cms.untracked.vstring()
)
if (not MC) :
    sourceFiles = cms.untracked.vstring( # 'root://cms-xrd-global.cern.ch/' prefix could help sometimes
            #'root://cmsxrootd.hep.wisc.edu//store/data/Run2012B/MuOniaParked/AOD/22Jan2013-v1/20003/CE2C94D8-036F-E211-A034-00215E21D8EE.root'
            '/store/data/Run2012B/MuOniaParked/AOD/22Jan2013-v1/20001/EE35A843-3E69-E211-9292-00215E2226AC.root'
    )
elif MC :
        if MCMotherId == 511 :
                if (not official) :
                        sourceFiles = cms.untracked.vstring(
                                #private
                                'file:/lustre/cms/store/group/cristella/Bd2Psi2SKpi-PHSP/MC_generation/141028_153606/merge/MC_Bd2Psi2SKpi.root'
                        )
		else :
                        sourceFiles = cms.untracked.vstring(
                                # offcial MC
                                'file:/lustre/cms/store/mc/Summer12DR53X/BdToPsi2SKPi_MSEL5_TuneZ2star_8TeV-pythia6/AODSIM/PU_RD2_START53_V19F-v1/10000/A8203BEC-BFC3-E411-BCE7-00266CFFC948.root'
           )
	elif MCMotherId == 531 :
                sourceFiles = cms.untracked.vstring(
		'/store/mc/Summer12_DR53X/BsToPsiMuMu_2MuPtEtaFilter_8TeV-pythia6-evtgen/AODSIM/PU_S10_START53_V7A-v1/0000/005DE3B0-FDDC-E111-9812-00266CFFC198.root',
	    )
process.PoolSource.fileNames = sourceFiles ;


process.source.inputCommands = cms.untracked.vstring(
        "keep *",
        "drop L1GlobalTriggerObjectMapRecord_hltL1GtObjectMap__RECO",
        "drop *_MEtoEDMConverter_*_*"
	)
process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32( 50 ) 
	#input = cms.untracked.int32( -1 ) 
	)
process.load('Configuration.Geometry.GeometryIdeal_cff') # 53x
process.load("Configuration.StandardSequences.GeometryExtended_cff") # from Lucia
process.load("Configuration.StandardSequences.Reconstruction_cff") # from Lucia
process.load("Configuration.StandardSequences.MagneticField_cff") # for using TransientTrackBuilder 
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff") # for using TransientTrackBuilder
process.GlobalTag.globaltag = 'FT_R_53_V18::All' #Global tag for 2012B data
#process.GlobalTag.globaltag = 'FT53_V21A_AN6::All' #Global tag for 2012C data
#process.GlobalTag.globaltag = '' #Global tag for 2012D data
process.load('Configuration/EventContent/EventContent_cff')
process.load('L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskAlgoTrigConfig_cff')
process.load('L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskTechTrigConfig_cff')
process.load('HLTrigger/HLTfilters/hltLevel1GTSeed_cfi')



#############################################################################################
##################################good collisions############################################
    
# 53x                                    
pvSelection = cms.PSet(
        minNdof = cms.double( 4. )
        , maxZ    = cms.double( 24. )
        , maxRho  = cms.double( 2. )
)
process.goodOfflinePrimaryVertices = cms.EDFilter("PrimaryVertexObjectFilter",      # checks for fake PVs automatically
                                                  filterParams = pvSelection,
                                                  filter       = cms.bool( False ), # use only as producer
                                                  src          = cms.InputTag( 'offlinePrimaryVertices' )
                                          )
process.primaryVertexFilter = process.goodOfflinePrimaryVertices.clone( filter = True )
process.noscraping = cms.EDFilter("FilterOutScraping",
                                  applyfilter = cms.untracked.bool(True),
                                  debugOn = cms.untracked.bool(False),
                                  #debugOn = cms.untracked.bool(True),
                                  numtrack = cms.untracked.uint32(10),
                                  thresh = cms.untracked.double(0.25)
                          )
process.load("PhysicsTools.PatAlgos.patSequences_cff")
process.load("PhysicsTools.PatAlgos.cleaningLayer1.genericTrackCleaner_cfi")
process.cleanPatTracks.checkOverlaps.muons.requireNoOverlaps = cms.bool(False)
process.cleanPatTracks.checkOverlaps.electrons.requireNoOverlaps = cms.bool(False)
from PhysicsTools.PatAlgos.producersLayer1.muonProducer_cfi import *
patMuons.embedTrack = cms.bool(True)
patMuons.embedPickyMuon = cms.bool(False)
patMuons.embedTpfmsMuon = cms.bool(False)

# Prune generated particles to muons and their parents
process.genMuons = cms.EDProducer("GenParticlePruner",
                                  src = cms.InputTag("genParticles"),
                                  select = cms.vstring(
                                          "drop  *  ",                     # this is the default
                                          "++keep abs(pdgId) = 13",        # keep muons and their parents
                                          "drop pdgId == 21 && status = 2" # remove intermediate qcd spam carrying no flavour info
                                  )
 )
process.load("MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff")
from MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff import  addMCinfo, useExistingPATMuons, useL1MatchingWindowForSinglets, changeTriggerProcessName, switchOffAmbiguityResolution, addDiMuonTriggers
    # with some customization
if MC:
        addMCinfo(process)
        # since we match inner tracks, keep the matching tight and make it one-to-one
        process.muonMatch.maxDeltaR = 0.05
        process.muonMatch.resolveByMatchQuality = True

addDiMuonTriggers(process)
useExistingPATMuons(process,'cleanPatMuons',addL1Info=False)
changeTriggerProcessName(process, 'HLT')
switchOffAmbiguityResolution(process) # Switch off ambiguity resolution: allow multiple reco muons to match to the same trigger muon
useL1MatchingWindowForSinglets(process)

process.muonL1Info.maxDeltaR     = 0.3
process.muonL1Info.fallbackToME1 = True
process.muonMatchHLTL1.maxDeltaR     = 0.3
process.muonMatchHLTL1.fallbackToME1 = True
process.muonMatchHLTL2.maxDeltaR = 0.3
process.muonMatchHLTL2.maxDPtRel = 10.0
process.muonMatchHLTL3.maxDeltaR = 0.1
process.muonMatchHLTL3.maxDPtRel = 10.0
process.muonMatchHLTCtfTrack.maxDeltaR = 0.1
process.muonMatchHLTCtfTrack.maxDPtRel = 10.0
process.muonMatchHLTTrackMu.maxDeltaR = 0.1
process.muonMatchHLTTrackMu.maxDPtRel = 10.0

from PhysicsTools.PatAlgos.tools.trackTools import *
######## adding tracks refitted with different mass
from RecoTracker.TrackProducer.TrackRefitters_cff import *
from TrackingTools.MaterialEffects.RungeKuttaTrackerPropagator_cfi import *
process.RungeKuttaTrackerPropagatorForPions = TrackingTools.MaterialEffects.RungeKuttaTrackerPropagator_cfi.RungeKuttaTrackerPropagator.clone( Mass = cms.double(0.13957), ComponentName = cms.string('RungeKuttaTrackerPropagatorForPions') )
process.refittedGeneralTracksPion = RecoTracker.TrackProducer.TrackRefitter_cfi.TrackRefitter.clone( Propagator = "RungeKuttaTrackerPropagatorForPions" )
makeTrackCandidates( process,                                # patAODTrackCands
                     label = 'TrackCands',                   # output collection will be 'allLayer0TrackCands', 'allLayer1TrackCands', 'selectedLayer1TrackCands'
                     tracks = cms.InputTag('generalTracks'), # input track collection
                     particleType = 'pi+',                   # particle type (for assigning a mass) # not working, everything is a pion
                     preselection = 'pt > 0.4',              # preselection cut on candidates. Only methods of 'reco::Candidate' are available
                     selection = 'pt > 0.4 && p > 0.5',      # Selection on PAT Layer 1 objects ('selectedLayer1TrackCands')
		     isolation = {},                         # Isolations to use ('source':deltaR; set to {} for None)
                     isoDeposits = [],
                     mcAs = None                             # Replicate MC match as the one used for Muons
             );                                   

l1cands = getattr(process, 'patTrackCands')
l1cands.addGenMatch = False

######## adding tracks refitted with Kaon mass
#process.RungeKuttaTrackerPropagator.Mass = cms.double(0.493677)
process.RungeKuttaTrackerPropagatorForKaons = TrackingTools.MaterialEffects.RungeKuttaTrackerPropagator_cfi.RungeKuttaTrackerPropagator.clone(
        Mass = cms.double(0.493677), ComponentName = cms.string('RungeKuttaTrackerPropagatorForKaons') )
process.refittedGeneralTracksKaon = RecoTracker.TrackProducer.TrackRefitter_cfi.TrackRefitter.clone( Propagator = "RungeKuttaTrackerPropagatorForKaons" )
###################################################
makeTrackCandidates( process,                                # patAODTrackCands
                     label = 'TrackKaonCands',               # output collection will be 'allLayer0TrackCands', 'allLayer1TrackCands', 'selectedLayer1TrackCands'
                     tracks = cms.InputTag('generalTracks'), # input track collection               // AP changed from generalTracks
                     particleType = 'K+',                    # particle type (for assigning a mass)  // AP changed from pi to K # not working, everything is a pion
                     preselection = 'pt > 0.4',              # preselection cut on candidates. Only methods of 'reco::Candidate' are available
                     selection = 'pt > 0.4 && p > 0.5',      # Selection on PAT Layer 1 objects ('selectedLayer1TrackCands')
		     isolation = {},                         # Isolations to use ('source':deltaR; set to {} for None)
                     isoDeposits = [],
                     mcAs = None                             # Replicate MC match as the one used for Muons
             );                                                      

l1cands = getattr(process, 'patTrackKaonCands')
l1cands.addGenMatch = False

process.load("RecoTracker.DeDx.dedxHarmonic2_cfi")
process.dedxHarmonic2Kaon = RecoTracker.DeDx.dedxHarmonic2_cfi.dedxHarmonic2.clone (
        tracks = 'refittedGeneralTracksKaon',
        trajectoryTrackAssociation = 'refittedGeneralTracksKaon'
)
# dE/dx hits
process.load("RecoVertex.BeamSpotProducer.BeamSpot_cff")
process.PATfilter = cms.EDFilter("Y4140FilterPAT")

process.mkcands = cms.EDAnalyzer("MuMuKKPAT",
                                 HLTriggerResults = cms.untracked.InputTag("TriggerResults","","HLT"),
                                 inputGEN  = cms.untracked.InputTag("genParticles"),
                                 VtxSample   = cms.untracked.string('offlinePrimaryVertices'),
                                 SameSign = cms.untracked.bool(False),
                                 DoMonteCarloTree = cms.untracked.bool( MC ),
                                 MonteCarloParticleId = cms.untracked.int32(20443),
                                 MonteCarloExclusiveDecay = cms.untracked.bool( MCExclusiveDecay ),
                                 MonteCarloMotherId = cms.untracked.int32( MCMotherId ), 
                                 MonteCarloDaughtersN = cms.untracked.int32( 3 ), # 3 for exclusive B0->psi'KPi
                                 
                                 DoMuMuMassConstraint = cms.untracked.bool(True),
                                 SkipJPsi = cms.untracked.bool(False),
                                 SkipPsi2S = cms.untracked.bool(False), 
                                 MinNumMuPixHits = cms.untracked.int32(1),
                                 MinNumMuSiHits = cms.untracked.int32(8),
                                 MaxMuNormChi2 = cms.untracked.double(7),
                                 MaxMuD0 = cms.untracked.double(10.0),
                                 sharedFraction = cms.untracked.double(0.5),

                                 MinJPsiMass = cms.untracked.double(2.8),       # SEMRA changed
                                 MaxJPsiMass = cms.untracked.double(3.4),       # SEMRA changed
				 MinPhiMass = cms.untracked.double (0.97),      # SEMRA added
 				 MaxPhiMass = cms.untracked.double (1.07),      # SEMRA added
				 MaxJPsiPhiXMass = cms.untracked.double (4.35), # SEMRA added
				 MinJPsiPhiB0Mass = cms.untracked.double (5.1), # SEMRA added
				 MaxJPsiPhiB0Mass = cms.untracked.double (5.6), # SEMRA added 

                                 MinNumTrSiHits = cms.untracked.int32(4),
                                 MinTrPt = cms.untracked.double(0.350),
                                 Chi2NDF_Track =  cms.untracked.double(7.0),
				 
				 MaxMuMuTrackDR = cms.untracked.double(1.5), 
                                 MaxB0CandTrackDR = cms.untracked.double(1.5),     
                                 UseB0Dr = cms.untracked.bool(True),             

                                 resolvePileUpAmbiguity = cms.untracked.bool(False),
                                 addMuMulessPrimaryVertex = cms.untracked.bool(True),
                                 #addMuMulessPrimaryVertex = cms.untracked.bool(False),
                                 addB0lessPrimaryVertex = cms.untracked.bool(True), 
                                 Debug_Output = cms.untracked.bool(False), # true
                                
                                 # Trigger path
                                 TriggersForMatching = cms.untracked.vstring(
                                         "HLT_DoubleMu4_Jpsi_Displaced_v9","HLT_DoubleMu4_Jpsi_Displaced_v10", "HLT_DoubleMu4_Jpsi_Displaced_v11", "HLT_DoubleMu4_Jpsi_Displaced_v12", #for Bs0
					 "HLT_Dimuon8_JPsi_v3", "HLT_Dimuon8_Jpsi_v4", "HLT_Dimuon8_Jpsi_v5", "HLT_Dimuon8_Jpsi_v6", "HLT_Dimuon8_Jpsi_v7", # for Y4140
                                 ),
				 FiltersForMatching = cms.untracked.vstring(
                                         "hltDisplacedmumuFilterDoubleMu4Jpsi", "hltDisplacedmumuFilterDoubleMu4Jpsi", "hltDisplacedmumuFilterDoubleMu4Jpsi", "hltDisplacedmumuFilterDoubleMu4Jpsi"
					 "hltVertexmumuFilterDimuon8Jpsi", "hltVertexmumuFilterDimuon8Jpsi", "hltVertexmumuFilterDimuon8Jpsi", "hltVertexmumuFilterDimuon8Jpsi", "hltVertexmumuFilterDimuon8Jpsi",
                                )
                                 
                                 
                         )


process.TFileService = cms.Service("TFileService",
    fileName = cms.string('set_below.root')
)
if (not MC) :
    process.TFileService.fileName = cms.string('MuOniaParked_Run2012B_MuMuKKPAT_ntpl.root')
elif MC :
    if MCMotherId == 511 :
            if (not official) :
                    process.TFileService.fileName = cms.string('BdToPsiKpi_18Mar_MuMuPiKPAT_ntpl.root')
            else :
                    process.TFileService.fileName = cms.string('officialBdToPsiKpi_18Mar_MuMuPiKPAT_ntpl.root')
    elif MCMotherId == 531 :
        process.TFileService.fileName = cms.string('BsToPsiMuMu_03Mar_MuMuPiKPAT_ntpl.root')


# turn off MC matching for the process
from PhysicsTools.PatAlgos.tools.coreTools import *
# old: removeMCMatching(process, ['All'], outputInProcess = False)
removeMCMatching(process,['All'],"",None,[])

process.patDefaultSequence.remove(process.patJetCorrFactors)
process.patDefaultSequence.remove(process.patJetCharge)
process.patDefaultSequence.remove(process.patJetPartonMatch)
process.patDefaultSequence.remove(process.patJetGenJetMatch)
process.patDefaultSequence.remove(process.patJetPartons)
process.patDefaultSequence.remove(process.patJetFlavourAssociation)
process.patDefaultSequence.remove(process.patJets)
process.patDefaultSequence.remove(process.patMETs)
process.patDefaultSequence.remove(process.selectedPatJets)
process.patDefaultSequence.remove(process.cleanPatJets)
process.patDefaultSequence.remove(process.countPatJets)

process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string('onia2MuMuPAT.root'),
                               outputCommands = cms.untracked.vstring('drop *',
                                             'keep patMuons_patMuonsWithTrigger_*_NTUPLE', # All PAT muons including general tracks and matches to triggers
                                                              )
                       )

process.filter = cms.Sequence(
        process.goodOfflinePrimaryVertices
        + process.primaryVertexFilter
        + process.noscraping
)

process.ntup = cms.Path(
        process.offlineBeamSpot #* process.dedxHitInfo
        * process.filter
        * process.patDefaultSequence
        * process.patMuonsWithTriggerSequence
        * process.PATfilter
        * process.mkcands
)

process.schedule = cms.Schedule(process.ntup)


