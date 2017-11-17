import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras

# Define the output directory
base            = '/afs/cern.ch/work/k/kmei/myAnalysis/CMSSW_9_2_10/src/'
outputDir       = base+'output/'

process = cms.Process("myAnalysis", eras.Run2_25ns)

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        #'file:myfile.root'
		'file:/eos/user/k/kmei/XXTo4J_M-300_CTau-30mm_reco_102_1_ne5.root'
    )
)

process.myAnalysis = cms.EDAnalyzer('TrackAnalyzer',
	tracks = cms.InputTag('generalTracks'),
	triggers = cms.InputTag('TriggerResults','','HLT'),
	caloJets = cms.InputTag('ak4CaloJets'),
	genParticles = cms.InputTag('genParticles'),
	primaryVertices = cms.InputTag('offlinePrimaryVertices'),
	secondaryVertices = cms.InputTag('inclusiveSecondaryVertices')
)

process.myAnalysis.outputFileName = cms.untracked.string("%smyOutput.root" % outputDir)

process.TFileService = cms.Service("TFileService",
	fileName = cms.string("%smyOtherOutput.root" % outputDir)
)


process.p = cms.Path(process.myAnalysis)
