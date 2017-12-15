typedef std::vector<DisplacedTrack> DisplacedTrackCollection;
typedef std::vector<Displaced2TrackVertex> DisplacedV0Collection;

class DisplacedJet {
    
	public:
	
	DisplacedJet( const reco::CaloJet   &jet_,
	              const reco::Vertex    &primaryVertex,
				  const bool            &isMC_,
				  const int             &jetID_,
				  const edm::EventSetup &iSetup_,
				  const int             &debug_ ) :
				  debug(debug_),
				  iSetup(iSetup_),
				  jet(jet_),
				  isMC(isMC_),
				  jetID(jetID_),
				  selPV(primaryVertex) {
	    if( debug ) std::cout<<"[DEBUG] In the DisplacedJet code and creating a displaced jet with pT="<<jet_.pt()<<", eta="<<jet_.eta()<<", phi="<<jet_.phi()<<std::endl;

		//Initialize variables in the displaced jet class
		
		//General track association variables

		nTracks                         = 0;
		nTracksPrompt                   = 0;
		nTracksDisp                     = 0;
		
        //Regional Tracking Counts

		nTracksRegPrompt               = -1;
		nTracksRegDisp                 = -1;
		nTracksRegPromptUp             = -1;
		nTracksRegPromptDown           = -1;
		nTracksRegDispUp               = -1;
		nTracksRegDispDown             = -1;

		//Initialize all the caloJet variables

		caloPt    = jet.pt();
		caloEta   = jet.eta();
		caloPhi   = jet.phi();
		caloPx    = jet.px();
		caloPy    = jet.py();
		caloPz    = jet.pz();

	    //Initialize the jet vertex fraction variable
		alpha                          = 0;
		alphaMax                       = 0;

	}

    //Initialization variables defined by the object
	const int              debug;
	const edm::EventSetup  &iSetup;
	const reco::CaloJet    jet;
	const bool             isMC;
	const int              jetID;
	const reco::Vertex     selPV;

    //Track Association Variables - Reg stands for regional, Disp stands for displaced
    int nTracks;
	int nTracksPrompt;
	int nTracksDisp;
	int nTracksRegPrompt;
	int nTracksRegDisp;
	int nTracksRegPromptUp;
	int nTracksRegPromptDown;
	int nTracksRegDispUp;
	int nTracksRegDispDown;

    float sumTrackPt;

	//Calojet Related Parameters
	float caloPt;
	float caloEta;
	float caloPhi;
	float caloPx;
	float caloPy;
	float caloPz;
	
	float detPt;
	float detEta;
	float detPhi;
	float caloEMEnergyFrac;
	float calohadEnergyFrac;

	float alpha;
	float alphaMax;

	//Jet Impact Parameter Variables
	std::vector<float> ip3dVector;
	std::vector<float> ip3dsVector;
	std::vector<float> ip2dVector;
	std::vector<float> ip2dsVector;

    //Jet Info Extraction Functions
	DisplacedTrackCollection getDisplacedTracks()     { return displacedTracks; }
	reco::TrackCollection    getVertexMatchedTracks() { return vertexMatchedTracks; }
	reco::Vertex             getIVFVertexSelected()   { return selIVF; }
	reco::Vertex             getSVVertex()            { return selSV; }

    //Define functions for this class
	void calcJetAlpha( const reco::TrackCollection&,
	                   const reco::VertexCollection& );
	void addVertexTrackInfo( const reco::TrackRefVector& );
	void addIPTagInfo( const reco::TrackIPTagInfo& );

	//Track Counting Related Variables
	void addTrackAngles(const DisplacedTrackCollection&);
	
	float jetNTracksNoPixel;
	float jetNTracksPixel;
	float jetPtSumTracksNoPixel;
	float jetPtSumTracksPixel;

	float ptSumCosTheta2D;
	float ptSumCosTheta3D;
	float ptSumCosThetaDet2D;
	float ptSumCosThetaDet3D;

	float sumCosTheta2D;
	float sumCosTheta3D;
	float sumCosThetaDet2D;
	float sumCosThetaDet3D;

	float meanCosTheta2D;
	float meanCosTheta3D;
	float meanCosThetaDet2D;
	float meanCosThetaDet3D;

	float medianCosTheta2D;
	float medianCosTheta3D;
	float medianCosThetaDet2D;
	float medianCosThetaDet3D;
	
	float trackSumMomCosTheta2D;
	float trackSumMomCosTheta3D;
	float trackSumMomMag2D;
	float trackSumMomMag3D;

	float ipPosSumMag2D;
	float ipPosSumMag3D;

	std::vector<float> cosTheta2DVector;
	std::vector<float> cosTheta3DVector;
	std::vector<float> cosThetaDet2DVector;
	std::vector<float> cosThetaDet3DVector;

	std::vector<float> trackPtVector;
	std::vector<float> trackEtaVector;
	std::vector<float> trackPhiVector;

	//Hit Related Function and Variables
	void addHitInfo( const reco::TrackCollection tracks );
	
	float jetMedianInnerHitPos;
	float jetMedianOuterHitPos;
	float jetMeanInnerHitPos;
	float jetMeanOuterHitPos;

	float jetMedianInnerHitPosInPixel;
	float jetMedianOuterHitPosInPixel;
	float jetMeanInnerHitPosInPixel;
	float jetMeanOuterHitPosInPixel;

	float jetMedianInnerHitPosOutPixel;
	float jetMedianOuterHitPosOutPixel;
	float jetMeanInnerHitPosOutPixel;
	float jetMeanOuterHitPosOutPixel;

	float jetMedianTrackValidHitFrac;
	float jetMeanTrackValidHitFrac;

	// Jet Distribution Calculator
	float getJetMedian( const std::vector<float>&, bool );
	float getJetMean( const std::vector<float>&, bool );

    DisplacedTrackCollection displacedTracks;

	private:

    reco::TrackCollection vertexMatchedTracks;
	reco::TrackRefVector  vertexMatchedTrackRefs;
	reco::Vertex          selIVF;
	reco::Vertex          selSV;

	std::vector<reco::btag::TrackIPData> lifetimeIPData;

};
//Function to add the impact parameter info to each displaced jet
void DisplacedJet::addIPTagInfo( const reco::TrackIPTagInfo &ipTagInfo ) {
	if( debug ) std::cout<<"[DEBUG] Adding secondary vertex IP info into the displaced jet"<<std::endl;

	//Pull the Impact Parameter Data
	lifetimeIPData = ipTagInfo.impactParameterData();

	//Loop over the IP info on each track
	for( std::vector<reco::btag::TrackIPData>::const_iterator itIP = lifetimeIPData.begin(); itIP != lifetimeIPData.end(); ++itIP ) {
		if( debug ) std::cout<<"[DEBUG] Filling up IP Info"<<std::endl;

		float ip3d  = itIP->ip3d.value();
		float ip3ds = itIP->ip3d.significance();
		float ip2d  = itIP->ip2d.value();
		float ip2ds = itIP->ip2d.significance();
		
		float jetAxisDis    = itIP->distanceToJetAxis.value();
		float jetAxisDisSig = itIP->distanceToJetAxis.significance();

		//Fill the vectors
		ip3dVector.push_back(ip3d);
		ip2dVector.push_back(ip2d);
		ip3dsVector.push_back(ip3ds);
		ip2dsVector.push_back(ip2ds);

		std::cout<<jetAxisDis<<" "<<jetAxisDisSig<<" "<<std::endl;
	}
}


//Function to count the number fo tracks based on the association at the vertex
void DisplacedJet::addVertexTrackInfo( const reco::TrackRefVector &trackRefs ) {
    
	//Initialize some variables
	vertexMatchedTrackRefs = trackRefs;
	nTracks                = 0;
	sumTrackPt             = 0;
	nTracksPrompt          = 0;
	nTracksDisp            = 0;

    if( debug ) std::cout<<"[DEBUG] Adding vertex matched track info "<<std::endl;

	reco::TrackRefVector::const_iterator itTrack = trackRefs.begin();
	for( ; itTrack != trackRefs.end(); ++itTrack ) {

	    //Build the Displaced Track object
	    DisplacedTrack dTrack( *itTrack, selPV, iSetup, debug);
		displacedTracks.push_back(dTrack);
	

	//Apply a pT cut to the vertexMatchedTrackCollection
	
	    float pt = (*itTrack)->pt();
	    if( pt > 1.0 ) {
	        nTracks++; //Increment the number of tracks
		    sumTrackPt += pt;
		    vertexMatchedTracks.push_back(**itTrack);

		    //Increment 
	  	    if( dTrack.ip2d < 0.1 ) nTracksPrompt++;
		    if( dTrack.ip2d > 0.05 && std::fabs( dTrack.ip2dSig ) > 5.0 ) nTracksDisp++;
		}
	}
}

//Function to comput variables related to the track angles and the caloJets

void DisplacedJet::addHitInfo( const reco::TrackCollection tracks ) {
	
	//Iterate over the tracks and store in these vectors
	std::vector<float> fractionValidHits;
	std::vector<float> innerHitPos;
	std::vector<float> outerHitPos;
	std::vector<float> innerHitPosInPixel;
	std::vector<float> innerHitPosOutPixel;
	std::vector<float> outerHitPosInPixel;
	std::vector<float> outerHitPosOutPixel;

	jetNTracksNoPixel     = 0;
	jetNTracksPixel       = 0;
	jetPtSumTracksNoPixel = 0;
	jetPtSumTracksPixel   = 0;


	for( reco::TrackCollection::const_iterator itTrack = tracks.begin(); itTrack != tracks.end(); ++itTrack ) {
		static GetTrackTrajInfo getTrackTrajInfo;
		const reco::Track &const_track = *itTrack;
		std::vector<GetTrackTrajInfo::Result> trajInfo = getTrackTrajInfo.analyze( iSetup, const_track );

		//Calculate the fraction of valid hits
		
		float trackNValidHits = 0.0;
		float pt = itTrack->pt();

		for ( std::vector<GetTrackTrajInfo::Result>::const_iterator itTraj = trajInfo.begin(); itTraj != trajInfo.end(); ++itTraj ) {
			if ( (*itTraj).valid ) {
				trackNValidHits++;
			}
		}

		fractionValidHits.push_back( float(trackNValidHits) / float(trajInfo.size()) );

		// Check the inner and outer hits
		if( trajInfo[0].valid && trajInfo.back().valid ) {
			
			// Get the state on the surface
			const TrajectoryStateOnSurface& tsosInnerHit = trajInfo[0].detTSOS;
			const TrajectoryStateOnSurface& tsosOuterHit = trajInfo.back().detTSOS;

			// Get the position of the hit
			const GlobalPoint&				innerPos     = tsosInnerHit.globalPosition();
			const GlobalPoint&              outerPos     = tsosOuterHit.globalPosition();

			TVector3 outerV3( outerPos.x(), outerPos.y(), outerPos.z() );
			TVector3 outerV2( outerPos.x(), outerPos.y(), 0.0 );
			TVector3 ipV3( itTrack->vx(), itTrack->vy(), itTrack->vz() );
			TVector3 ipV2( itTrack->vx(), itTrack->vy(), 0.0 );
			TVector3 innerV3 ( innerPos.x(), innerPos.y(), innerPos.z() );
			TVector3 innerV2 ( innerPos.x(), innerPos.y(), 0.0 );

			float innerRadius = innerV2.Mag();
			float outerRadius = outerV2.Mag();

			innerHitPos.push_back( innerRadius );
			outerHitPos.push_back( outerRadius );

			bool hasPixelHits = itTrack->hitPattern().numberOfValidPixelHits() > 0;

			//If the inner hit is in the pixel layers
			if( hasPixelHits ) {
				innerHitPosInPixel.push_back( innerRadius );
				outerHitPosInPixel.push_back( outerRadius );
				jetNTracksPixel++;
				jetPtSumTracksPixel += pt;
			}
			else {
				innerHitPosOutPixel.push_back( innerRadius );
				outerHitPosOutPixel.push_back( outerRadius );
				jetNTracksNoPixel++;
				jetPtSumTracksNoPixel += pt;
			}
		}
	}
	bool check_medians = false;
	if( check_medians ) {
		//ADD CODE TO CHECK .getJetMedian if necessary
		std::cout<<"No need to check get jet median (unless something comes up that contradicts this sentiment)"<<std::endl;
	}

	jetMedianInnerHitPos 		 = getJetMedian( innerHitPos, false );
	jetMedianOuterHitPos         = getJetMedian( outerHitPos, false );
	jetMeanInnerHitPos			 = getJetMean( innerHitPos, false );
	jetMeanOuterHitPos			 = getJetMean( outerHitPos, false );

	jetMedianInnerHitPosInPixel  = getJetMedian( innerHitPosInPixel, false );
	jetMedianOuterHitPosInPixel  = getJetMedian( outerHitPosInPixel, false );
	jetMeanInnerHitPosInPixel    = getJetMean( innerHitPosInPixel, false );
	jetMeanOuterHitPosInPixel    = getJetMean( outerHitPosInPixel, false );

	jetMedianInnerHitPosOutPixel = getJetMedian( innerHitPosOutPixel, false );
	jetMedianOuterHitPosOutPixel = getJetMedian( outerHitPosOutPixel, false );
	jetMeanInnerHitPosOutPixel   = getJetMean( innerHitPosOutPixel, false );
	jetMeanOuterHitPosOutPixel   = getJetMean( outerHitPosOutPixel, false );

	jetMedianTrackValidHitFrac   = getJetMedian( fractionValidHits, false );
	jetMeanTrackValidHitFrac     = getJetMean( fractionValidHits, false );
}

//Compute Variables Related to the Track Angles and the Calo Jet
void DisplacedJet::addTrackAngles( const DisplacedTrackCollection& tracks ) {
	
	//PT Weighted Angle Distribution
	ptSumCosTheta2D 		= 0;
	ptSumCosTheta3D 		= 0;
	ptSumCosThetaDet2D 		= 0;
	ptSumCosThetaDet3D		= 0;

	//Absolute Sum Angle Distribution
	sumCosTheta2D			= 0;
	sumCosTheta3D			= 0;
	sumCosThetaDet2D		= 0;
	sumCosThetaDet3D		= 0;

	//Find the Angle between the Jet Momentum and Track Momentum
	
	float sumTrackPtValid = 0;
	for( DisplacedTrackCollection::const_iterator itTrack = tracks.begin(); itTrack != tracks.end(); ++itTrack ) {
		if( itTrack->isValid ) {
			float pt = itTrack->pt;
			if( pt < 1 ) continue; //Only look at tracks with at least 1 GeV of pT

			//Track cosine with respect to the momentum at the reference point
			//Momentum is already defined at the reference point

			float cosTheta3D 	= itTrack->angleMomentumAndPVAtInnerHit3D;
			float cosTheta2D	= itTrack->angleMomentumAndPVAtInnerHit2D;
			
			float cosThetaDet3D = itTrack->angleMomentumAndPVAtOuterHit3D;
			float cosThetaDet2D = itTrack->angleMomentumAndPVAtOuterHit2D;

			trackEtaVector.push_back( itTrack->eta );
			trackPhiVector.push_back( itTrack->phi );
			trackPtVector.push_back( itTrack->pt );

			cosTheta2DVector.push_back( cosTheta2D );
			cosTheta3DVector.push_back( cosTheta3D );
			cosThetaDet2DVector.push_back( cosThetaDet2D );
			cosThetaDet3DVector.push_back( cosThetaDet3D );

			//Increment Sum Pt
			sumTrackPtValid += pt;

			//Pt Weighted Sum
			ptSumCosTheta2D 	+= sin( cosTheta2D ) * pt;
			ptSumCosTheta3D 	+= sin( cosTheta3D ) * pt;
			ptSumCosThetaDet2D 	+= sin( cosThetaDet2D ) * pt;
			ptSumCosThetaDet3D 	+= sin( cosThetaDet3D ) * pt;

			//Not Pt Weighted Sum
			sumCosTheta2D		+= sin( cosTheta2D );
			sumCosTheta3D		+= sin( cosTheta3D );
			sumCosThetaDet2D	+= sin( cosThetaDet2D );
			sumCosThetaDet3D	+= sin( cosThetaDet3D );
		}

		//Calculate Means
		meanCosTheta2D 		= getJetMean( cosTheta2DVector, true );
		meanCosTheta3D		= getJetMean( cosTheta3DVector, true );
		meanCosThetaDet2D	= getJetMean( cosThetaDet2DVector, true );
		meanCosThetaDet3D	= getJetMean( cosThetaDet3DVector, true );

		//Calculate Medians
		medianCosTheta2D	= getJetMedian( cosTheta2DVector, true );
		medianCosTheta3D	= getJetMedian( cosTheta3DVector, true );
		medianCosThetaDet2D	= getJetMedian( cosThetaDet2DVector, true );
		medianCosThetaDet3D = getJetMedian( cosThetaDet3DVector, true );
	}

	//Initialize Track Sum Variables
	trackSumMomCosTheta2D 	= 0.0;
	trackSumMomCosTheta3D 	= 0.0;
	trackSumMomMag2D		= 0.0;
	trackSumMomMag3D		= 0.0;

	//Initialize Impact Parameter Vector Sums
	ipPosSumMag2D			= 0.0;
	ipPosSumMag3D			= 0.0;

	//Normalize all variables by the sum Pt
	ptSumCosTheta2D 	/= sumTrackPtValid ? sumTrackPtValid : 1;
	ptSumCosTheta3D 	/= sumTrackPtValid ? sumTrackPtValid : 1;
	ptSumCosThetaDet2D	/= sumTrackPtValid ? sumTrackPtValid : 1;
	ptSumCosThetaDet3D 	/= sumTrackPtValid ? sumTrackPtValid : 1;
}

//Function to calculate the jet variable alpha, defined by Josh as the ratio of (Vertex Tracks Sum pT Matching the Jet) / (General Tracks Sum Pt Matching the Jet)
void DisplacedJet::calcJetAlpha( const reco::TrackCollection  &tracks,
                                 const reco::VertexCollection &primaryVertices ) { 

    //First calculate the sum track pT based on the tracks matched
	
	//Step 1. Take the scalar sum pT of tracks from the primary vertices relative to the total pT matched to the jets
	
	//Initialize some variables
	
	float sumJetPt               = 0;
	float sumJetPtMax            = 0;
	float sumJetPt_temp          = 0;

	float leadingVertexPtSq      = 0;
	float leadingVertexPtSq_temp = 0;

	//Loop over the vertices in the event (use beam spot constraint vertices since displaced jets are least likely to come from those)
	
	reco::VertexCollection::const_iterator itVertex = primaryVertices.begin();
	for( ; itVertex != primaryVertices.end(); ++itVertex ) {
	    //Calculate the sums for the specific vertex
		sumJetPt_temp          = 0;
		leadingVertexPtSq_temp = 0;

		std::vector<reco::TrackBaseRef>::const_iterator itTrack = itVertex->tracks_begin();
		
		//Loop over the tracks in the vertex
		for( ; itTrack != itVertex->tracks_end(); ++itTrack ) {
		    float tempPt = (*itTrack)->pt();

			if( (*itTrack)->pt() < 1.0 ) continue; // Used to apply a cut on the track pT if necessary
			leadingVertexPtSq_temp += tempPt * tempPt; //Sort by highest pt squared

			//Check for matching to the jet
			float dr = reco::deltaR( (*itTrack)->eta(), (*itTrack)->phi(), caloEta, caloPhi );
			if( dr < 0.4 ) sumJetPt_temp += (*itTrack)->pt();
		}
        //Step 2:Pick the vertex with the high sum pt squared
	
	    if( leadingVertexPtSq_temp > leadingVertexPtSq ) {
	        leadingVertexPtSq = leadingVertexPtSq_temp;
			sumJetPt = sumJetPt_temp;
	    } //Loop over the tracks from vertex

		//Also keep the highest posible sum
		sumJetPtMax = std::max(sumJetPtMax, sumJetPt);
	} //loop over vertices

	alpha       = sumJetPt;
	alphaMax    = sumJetPtMax;
	std::cout<<"AlphaMax: "<<alphaMax<<std::endl;
}

//The following three functions are needed for getHitInfo

float DisplacedJet::getJetMedian( const std::vector<float>& values, bool is_signed ) {
	if( values.size() == 0 ) return -999;

	int size = values.size();
	float *sorted = new float[size];

	for( int it = 0; it < size; ++it ) {
		sorted[it] = is_signed ? values[it] : std::fabs( values[it] );
	}

	//Sort the Array
	for( int jt = size - 1; jt > 0; --jt ) {
		for( int kt = 0; kt < jt; ++kt ) {
			if( sorted[kt] > sorted[kt+1] ) {
				float dTemp = sorted[kt];
				sorted[kt] = sorted[kt+1];
				sorted[kt+1] = dTemp;
			}
		}
	}

	//Middle or Average of Middle Values in the Sorted Array
	float dMedian = 0.0;
	if( ( size%2 ) == 0 ) dMedian = ( sorted[size/2] + sorted[ (size/2) - 1] ) / 2.0;
	else dMedian = sorted[size/2];

	delete [] sorted;

	//Safety Check
	if( !is_signed ) {
		if( dMedian < 0 ) std::cout<<"ASSERT FAIL: "<<dMedian<<std::endl;
		assert( dMedian >= 0 );
	}
	
	return dMedian;
}

float DisplacedJet::getJetMean( const std::vector<float> &values, bool is_signed) {
	if( values.size() == 0 ) return 0.0;

	float sum = 0.0;
	for( std::vector<float>::const_iterator itValue = values.begin(); itValue != values.end(); ++itValue )
		sum += is_signed? *itValue : std::fabs( *itValue );

	float mean = sum / float( values.size() );
	return mean;
}

