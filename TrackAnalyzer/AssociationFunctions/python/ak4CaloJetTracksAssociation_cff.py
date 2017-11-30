import FWCore.ParameterSet.Config as cms

from TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAlong_cfi import *

from RecoJets.JetAssociationProducers.j2tParametersCALO_cfi import *
from RecoJets.JetAssociationProducers.j2tParametersVX_cfi import *

#Run the jet track asssociator for the calojets - only need caloJets (ignore all parts of the Particle Flow associators in the original file)
ak4JetTracksAssociatorAtVertex = cms.EDProducer("JetTracksAssociatorAtVertex",
	j2tParametersVX,
	jets = cms.InputTag("ak4CaloJets")
)

ak4JetTracksAssociatorAtCaloFAce = cms.EDProducer("JetTracksAssociatorAtCaloFace",
	j2tParametersCALO,
	jets = cms.InputTag("ak4CaloJets")
)

#Define the extender as see fit
ak4JetExtender = cms.EDProducer("JetExtender",
	jets = cms.InputTag("ak4CaloJets"),
	jet2TracksAtCALO = cms.InputTag("ak4JetTracksAssociatorAtCaloFace"),
	jet2TracksAtVX = cms.InputTag("ak4JetTracksAssociatorAtVertex"),
	coneSize = cms.double(0.4)
)

#Define the final 
ak4JTA = cms.Sequence(ak4JetTracksAssociatorAtVertex*
					  ak4JetTracksAssociatorAtCaloFace*
					  ak4JetExtender)
