typedef std::vector<DisplacedJet> DisplacedJetCollection;


class DisplacedJetEvent {
    //Constructor designating the calojets, primary vertices, and kinematics cuts
	
	public: 
	
	DisplacedJetEvent( const bool &isMC,
	                   const reco::CaloJetCollection &caloJets,
					   const reco::VertexCollection &primaryVertices,
					   const float &minPT_,
					   const float &minEta_,
					   const edm::EventSetup &iSetup_,
					   const int &debug_ ) :
					   selPV(*(primaryVertices.begin())),
					   minPT(minPT_),
					   minEta(minEta_),
					   debug(debug_),
					   iSetup(iSetup_) {
	
        caloHT            = 0;
		caloMET           = 0;
		nIVFReconstructed = 0;
		nIVFGenMatched    = 0;

        //Construct the empty displaced jet objects to merge info later

        if( debug > 1 ) std::cout<<"[DEBUG 1] Construting the event from calojets"<<std::endl;
		
		reco::CaloJetCollection::const_iterator itCaloJet = caloJets.begin();
		caloLeadingJetPT    = -1;
		caloSubLeadingJetPT = -1;
		
		int jetIndex = 0;
		
		for( ; itCaloJet != caloJets.end(); ++itCaloJet, ++jetIndex ) {
            float pt  = itCaloJet->pt();
			float eta = itCaloJet->eta();

			if( pt > 40 && std::fabs(eta) < 3.0 )  caloHT += pt; //Only add jets with pT higher than 40 GeV and in either the endcap or the barrel
			if( pt < minPT || std::fabs(eta) > minEta ) continue; //Remove all jets with less PT than the minimum and if out of bounds of analysis #FIXME: Why do we do this?
			if( pt > 99999 || pt < 0 || isnan(pt) ) {
			    std::cerr<<"Badly defined jet pT - check event"<<std::endl;
				continue;
			}//Eliminate all jets with negative pT, a pT that's too large, or a nonexistent pT
			if( std::fabs(eta) > 10 || std::fabs(itCaloJet->phi()) > 3.142 ) {
			    std::cerr<<"Jet in a region that is not detectable by detector (either by eta: "<<std::fabs(eta)<<" or phi: "<<std::fabs(itCaloJet->phi())<<")"<<std::endl;
				continue;
			}

			//Now define the jets with the leading pT and the subleading pT
			
			//If no other jet has been saved
			if( pt  > caloLeadingJetPT && caloLeadingJetPT == -1 ) {
			    caloLeadingJetPT = pt;
			}
			//If another jet has been saved, compare and if current jet has higher pT
			else if( pt > caloLeadingJetPT && caloLeadingJetPT > 0 ) {
			    caloSubLeadingJetPT = caloLeadingJetPT;
				caloLeadingJetPT = pt;
			}

			else if( pt > caloSubLeadingJetPT ) {
			    caloSubLeadingJetPT = pt;
			}

			DisplacedJet djet( *itCaloJet, selPV, isMC, jetIndex, iSetup, debug );

			djets.push_back( djet );
			djetIndex++;
			djet.calcJetAlpha(djet.getVertexMatchedTracks(), primaryVertices);
        }
	}

	//Define the variables that are used in the class
	
	//Kinematic Variables:
	float caloHT;
    std::vector<float> caloDHT;
	float caloMET;

	//Ordered Quantities
	float caloLeadingJetPT;
	float caloSubLeadingJetPT;

    //IVF Related
	int nIVFReconstructed;
	int nIVFGenMatched;

	//Index variable for defining how large the displaced jet array is
	
	int djetIndex = 0;

	DisplacedJetCollection djets;

	//Track related variables used to add vertexing track informtion
	int nTracks       = 0;
	int nTracksDisp   = 0;
	int nTracksPrompt = 0;
	float sumTrackPt  = 0;

	//Variables Used in the Constructor
	const reco::Vertex     selPV;
	const float            minPT;
	const float            minEta;
	const int              debug;
	const edm::EventSetup& iSetup;

	//Variables used in the helper functions
    reco::TrackRefVector  vertexMatchedTrackRefs;
	reco::TrackCollection vertexMatchedTracks;
	reco::TrackCollection caloMatchedTracks;

	//Functions used for doing vertex matching
	reco::TrackRefVector   getVertexMatchedTrackRefs() { return vertexMatchedTrackRefs; }
	reco::TrackCollection  getVertexMatchedTracks()    { return vertexMatchedTracks; }
	reco::TrackCollection  getCaloMatchedTracks()      { return caloMatchedTracks; }
	//DisplacedJet& myDisplacedJet;

	//Defining publicly accessible functions
	DisplacedJetCollection& getDisplacedJets() { return djets; }
	void                    mergeTrackAssociations ( const reco::JetTracksAssociation::Container&,
													 const reco::JetTracksAssociation::Container& );
	DisplacedJet&           findDisplacedJetByPtEtaPhi( const float&, const float&, const float& );
};

//Need this to get add the vertex track info FIXME:Add the variables needed in this function to the class
void DisplacedJetEvent::mergeTrackAssociations( const reco::JetTracksAssociation::Container& caloMatched, const reco::JetTracksAssociation::Container& vertexMatched) {

	std::vector<reco::JetBaseRef> caloJets   = reco::JetTracksAssociation::allJets(caloMatched);
	std::vector<reco::JetBaseRef> vertexJets = reco::JetTracksAssociation::allJets(vertexMatched);

	if( caloJets.size() != vertexJets.size() ) {
		throw cms::Exception("TrackAssociationJetMatchingFailure");
	}

	for( std::vector<reco::JetBaseRef>::const_iterator itJetBaseRef = caloJets.begin(); itJetBaseRef != caloJets.end(); ++itJetBaseRef ) {
		//Writing code that finds the displaced jet by pt, eta and phi?
	}
}

DisplacedJet& DisplacedJetEvent::findDisplacedJetByPtEtaPhi( const float& pt,
															 const float& eta,
															 const float& phi ) {
	if( debug ) std::cout<<"[DEBUG] Looking for a displaced jet by pt, eta and phi"<<std::endl;

	bool displacedJetFound = false;
	std::vector<DisplacedJet>::iterator myDjet;

	for( std::vector<DisplacedJet>::iterator itDjet = djets.begin(); itDjet != djets.end(); ++itDjet ) {
		bool matchPt  = (itDjet->caloPt  == pt);
		bool matchEta = (itDjet->caloEta == eta);
		bool matchPhi = (itDjet->caloPhi == phi);

		if( matchPt && matchEta && matchPhi ) {
			displacedJetFound = true;
			myDjet = itDjet;
			break;
		}
	}
	assert(displacedJetFound);
	return *myDjet;
}
