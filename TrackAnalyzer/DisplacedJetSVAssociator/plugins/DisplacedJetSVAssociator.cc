//System Include Files

#include <memory>

//User Include Files

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "TrackAnalyzer/DisplacedJetSVAssociator/interface/JetVertexAssociation.h"


class DisplacedJetSVAssociator : public edm::EDProducer{

public:
	explicit DisplacedJetSVAssociator( const edm::ParameterSet& );
	~DisplacedJetSVAssociator();

	static void fillDescriptions( edm::ConfigurationDescriptions& descriptions );

	typedef std::vector<float> jetVertexScores;

private:
	virtual void beginJob() override;
	virtual void produce( edm::Event&, const edm::EventSetup& ) override;
	virtual void endJob() override;

	edm::InputTag caloJetsTag_;
	edm::InputTag secondaryVerticesTag_;
	edm::InputTag primaryVerticesTag_;

	std::string algorithmName_;
	std::string outputLabel_;
	float jetPtCut_;
	int debug_;
};


DisplacedJetSVAssociator::DisplacedJetSVAssociator( const edm::ParameterSet& iConfig ) {
	caloJetsTag_			= iConfig.getUntrackedParameter<edm::InputTag>( "caloJets" );
	secondaryVerticesTag_	= iConfig.getUntrackedParameter<edm::InputTag>( "secondaryVertices" );
	primaryVerticesTag_		= iConfig.getUntrackedParameter<edm::InputTag>( "primaryVertices" );
	algorithmName_			= iConfig.getUntrackedParameter<std::string>( "algoName" );
	jetPtCut_				= iConfig.getUntrackedParameter<double>( "jetPtCut" );
	debug_					= iConfig.getUntrackedParameter<int>( "debug" );
	outputLabel_			= iConfig.getUntrackedParameter<std::string>( "outputLabel" );

	produces<std::vector<std::vector<float> > >( outputLabel_ );
}

DisplacedJetSVAssociator::~DisplacedJetSVAssociator()
{
}

void DisplacedJetSVAssociator::produce( edm::Event& iEvent, const edm::EventSetup& iSetup ) {
	using namespace edm;

	std::auto_ptr<std::vector<jetVertexScores>> jetVertexAssociationCollection( new std::vector<jetVertexScores> );

	edm::Handle<reco::CaloJetCollection>	caloJets;
	edm::Handle<reco::VertexCollection>		secondaryVertices;
	edm::Handle<reco::VertexCollection>		primaryVertices;

	iEvent.getByLabel( caloJetsTag_, caloJets );
	iEvent.getByLabel( secondaryVerticesTag_, secondaryVertices );
	iEvent.getByLabel( primaryVerticesTag_, primaryVertices );

	const reco::CaloJetCollection&	jets	= *( caloJets.product() );
	const reco::VertexCollection& 	sv		= *( secondaryVertices.product() );
	const reco::VertexCollection&	pv		= *( primaryVertices.product() );

	const reco::Vertex primaryVertex = pv[0];
	
	for( reco::CaloJetCollection::const_iterator iJet = jets.begin(); iJet != jets.end(); ++iJet ) {
	
		math::XYZTLorentzVectorD	p4		= iJet->detectorP4();
		float						jet_pt	= p4.pt();
		float						jet_eta	= p4.eta();
		float						jet_phi	= p4.phi();

		if( jet_pt < jetPtCut_ ) continue;

		jetVertexScores vertexScores;
		float bestScore = -1;
		reco::Vertex bestVertex;

		for( reco::VertexCollection::const_iterator ss = sv.begin(); ss != sv.end(); ++ss ) {
			float score = 0;
			for( std::vector<reco::TrackBaseRef>::const_iterator tt = ss->tracks_begin(); tt != ss->tracks_end(); ++tt ) {
				float track_outerEta = (*tt)->outerEta();
				float track_outerPhi = (*tt)->outerPhi();

				float dR = reco::deltaR( jet_eta, jet_phi, track_outerEta, track_outerPhi );

				if( dR > 1.0 ) continue;
				else score += 1.0 / dR;
			}

			vertexScores.push_back( score );

			if( score > bestScore ) {
				bestScore = score;
				bestVertex = *ss;
			}
		}

		( *jetVertexAssociationCollection).push_back( vertexScores );
	}
	
	//iEvent.put( jetVertexAssociationCollection, outputLabel_ );
}

void DisplacedJetSVAssociator::beginJob() {
}

void DisplacedJetSVAssociator::endJob() {
}

void DisplacedJetSVAssociator::fillDescriptions( edm::ConfigurationDescriptions& descriptions ) {
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault( desc );
}

DEFINE_FWK_MODULE( DisplacedJetSVAssociator );
