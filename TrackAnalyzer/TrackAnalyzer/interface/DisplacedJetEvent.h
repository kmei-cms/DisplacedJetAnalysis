typedef std::vector<DisplacedJet> DisplacedJetCollection;

class DisplacedJetEvent {
    //Constructor designating the calojets, primary vertices, and kinematics cuts
	
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
		
		for( ; itCaloJet != caloJets.end(); ++jetIter, ++jetIndex ) {
            float pt  = jetIter->pt();
			float eta = jetIter->eta();

			if( pt > 40 && std::fabs(eta) < 3.0 )  caloHT += pt; //Only add jets with pT higher than 40 GeV and in either the endcap or the barrel
			if( pt < minPT || std::fabs(eta) > minEta ) continue; //Remove all jets with less PT than the minimum and if out of bounds of analysis #FIXME: Why do we do this?
			if( pt > 99999 || pt < 0 || isnan(pt) ) {
			    std::cerr<<"Badly defined jet pT - check event"<<std<<endl;
				continue;
			}//Eliminate all jets with negative pT, a pT that's too large, or a nonexistent pT
			if( std::fabs(eta) > 10 || std::fabs(itCaloJet->phi()) < 3.142 ) {
			    std::cerr<<"Jet in a region that is not detectable by detector (either by eta or phi"<<std::endl;
				continue;
			}

			//Now define the jets with the leading pT and the subleading pT
			
			//If no other jet has been saved
			if( pt  > caloLeadignJetPT && caloLeadingJetPT == -1 ) {
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

			DisplacedJet djet( *itCaloJet, selPV, isMC, jetIDCounter, iSetup, debug );

			djets.push_back( djet );
			djetIndex++;
        }
	}


	//Define the variables that are used in the class
	
	//Kinematic Variables:
	float caloHT;
    std::vector<float> caloDHT;
	float caloMET;

	//Ordered Quantities
	float caloLeadingJetPT;
	float caloSubLeadingJetPt;

    //IVF Related
	int nIVFReconstructed;
	int nIVFGenMatched;

	//Index variable for defining how large the displaced jet array is
	
	int djetIndex = 0;
	
}
