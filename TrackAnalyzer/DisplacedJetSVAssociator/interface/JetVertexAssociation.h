// Jet include files
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/CaloJet.h"

// Tracks include files
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

// Vertex include files
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BTauReco/interface/VertexTypes.h"

// Kinematics
#include "DataFormats/Math/interface/deltaR.h"

class JetVertexAssociation {

	public:

	JetVertexAssociation( std::string nameString, const reco::Vertex& PV, const int& debug_ ) {
		primaryVertex 		= PV;
		name				= nameString;
		debug				= debug_;
	}

	// Fillers
	void addCaloJet( const reco::CaloJet& );
	void addVertex( const reco::Vertex& );
	void setPrimaryVertex( const reco::Vertex& );

	// Accessors
	int getNVertices() 			{ return vertexCollection.size(); }
	int getNJets()				{ return caloJetCollection.size(); }
	float getBestVertexScore()	{ return bestVertexScore; }

	std::string getName() 		{ return name; }

	// Score
	const std::pair< const reco::Vertex, const float > getBestVertex( const reco::CaloJet&, const std::string& );
	float getVertexJetScore( const reco::CaloJet&, const reco::Vertex&, const std::string& );

	private:

	reco::Vertex			primaryVertex;
	reco::CaloJetCollection	caloJetCollection;
	reco::VertexCollection	vertexCollection;
	reco::Vertex			bestVertex;
	float					bestVertexScore;
	std::string				name;
	int						debug;
};

void JetVertexAssociation::setPrimaryVertex( const reco::Vertex& pv ) {
	primaryVertex = pv;
}

void JetVertexAssociation::addCaloJet( const reco::CaloJet& jet ) { caloJetCollection.push_back( jet ); }

void JetVertexAssociation::addVertex( const reco::Vertex& vertex ) { vertexCollection.push_back(vertex); }

float JetVertexAssociation::getVertexJetScore( const reco::CaloJet& jet, const reco::Vertex& vertex, const std::string& algo ) {
	math::XYZTLorentzVectorD 	p4 	= jet.detectorP4();
	float						jet_eta = p4.eta();
	float						jet_phi = p4.phi();
	float						score = 0;

	for( reco::Vertex::trackRef_iterator tt = vertex.tracks_begin(); tt!= vertex.tracks_end(); ++tt) {
		float eta = (*tt)->eta();
		float phi = (*tt)->phi();

		float dR  = reco::deltaR( jet_eta, jet_phi, eta, phi );

		if( dR > 1.0 ) continue;
		else score += 1.0/dR;
	}

	return score;
}

const std::pair<const reco::Vertex, const float> JetVertexAssociation::getBestVertex( const reco::CaloJet& jet, const std::string& algo ) {
	float bestScore = -1;
	reco::Vertex bestVertex;
	
	for( reco::VertexCollection::const_iterator ss = vertexCollection.begin(); ss != vertexCollection.end(); ++ss ) {
		float score = getVertexJetScore( jet, *ss, algo );

		if( score > bestScore ) {
			bestScore = score;
			bestVertex = *ss;
		}
	}

	if( bestScore > 0 ) {
		const std::pair<const reco::Vertex, const float> vertexPair( bestVertex, bestScore );
		return vertexPair;
	}

	const std::pair<const reco::Vertex, const float> vertexPair( primaryVertex, 0 );
	return vertexPair;
}
