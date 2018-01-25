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
	DisplacedTrackCollection getDisplacedTracks()     	 { return displacedTracks; }
	reco::TrackCollection    getVertexMatchedTracks() 	 { return vertexMatchedTracks; } 
	reco::TrackRefVector	 getVertexMatchedTrackRefs() { return vertexMatchedTrackRefs; }
	reco::Vertex             getIVFVertexSelected()   	 { return selIVF; }
	reco::Vertex             getSVVertex()            	 { return selSV; }

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

	//Vertex Related Function and Variables
	void addV0Info( const reco::TrackRefVector tracks );

	int jetOneTrackNuclearCount;
	int jetTwoTrackNuclearCount;
	int jetVertexNearBPIX1;
	int jetVertexNearBPIX2;
	int jetVertexNearBPIX;
	int jetTightNuclear;
	int jetLooseNuclear;
	int jetNV0NoHitBehindVertex;
	int jetNV0HitBehindVertex;
	int jetNV0KShort;
	int jetNV0Lambda;
	int jetV0HIndex;

	//Cluster Related Variables
	DisplacedCluster *v0Cluster = NULL;
	int	  jetV0ClusterSize;
	float jetV0ClusterLxy;
	float jetV0ClusterLxySig;
	float jetV0ClusterLxyz;
	float jetV0ClusterLxyzSig;
	float jetV0ClusterX;
	float jetV0ClusterY;
	float jetV0ClusterZ;
	float jetV0ClusterChi2;
	float jetV0ClusterIntercept;
	float jetV0ClusterAngle;
	float jetV0ClusterAngleMom;
	int   jetV0ClusterNTracks;

	//IVF Related Variables
	float selIVFIsPVScore;
	bool  selIVFIsPV;
	float ivfX;
	float ivfY;
	float ivfZ;
	float ivfXError;
	float ivfYError;
	float ivfZError;
	int   ivfNTracks;
	float ivfMass;
	float ivfLxySig;
	float ivfLxyzSig;
	float ivfLxy;
	float ivfLxyz;
	float ivfMatchingScore;

	
	//Clique Related Function and Variables
	typedef std::vector<std::vector<int>> vertexGraph;
	vertexGraph buildInputGraph( const std::vector<std::pair<int,int>>& graph_edges, const std::vector<int>& uniq_tracks, const bool& print);
	void findVertexCliques();
	int findRefTrack( const reco::TrackRef ref, const reco::TrackRefVector vector );
	int calcHIndex( const vertexGraph& graph );
	void calcClusterSize( const DisplacedV0Collection& vertices, const float& errorWindow );
	//void calcNJetClusterSize( const DisplacedV0Collection& vertices, std::vector<DisplacedJet>& djets, const float& errorWindow );
	void addSVTagInfo( const reco::SecondaryVertexTagInfo& );
	void addIVFCollection( const reco::VertexCollection& );
	
	//Secondary Vertex Related Variables
	int svNVertex;
	int svNTracks;
	int scNDof;
	float svX;
	float svY;
	float svZ;
	float svChi2;
	float svXError;
	float svYError;
	float svZError;

	float svAngle2D;
	float svAngle3D;
	float svMass;
	float svPt;
	float svEta;
	float svPhi;
	float svLxyzSig;
	float svLxySig;
	float svLxyz;
	float svLxy;

	//Combinatorial Vertices for the Jet
	vertexGraph v0Graph;
	DisplacedV0Collection displacedV0Vector;
	DisplacedV0Collection displacedV0VectorCleaned;
	std::vector<TransientVertex> transientV0Vector;
	std::vector<TransientVertex> transientV0VectorCleaned;


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

	float metric2D( float x, float y ) {
		return std::sqrt( x*x + y*y );
	}

	float metric3D( float x, float y, float z ) {
		return std::sqrt ( x*x + y*y + z*z );
	}

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

void DisplacedJet::addSVTagInfo( const reco::SecondaryVertexTagInfo& svTagInfo ) {
	svNVertex = svTagInfo.nVertices();

	//Track the Best Reconstructed Vertex
	int svIndex		= 0;
	int mostTracks	= 0;
	int tieBreaker	= 0;

	//Loop Over All the Reconstructed Vertices and Choose the Best One
	for( int itVertex = 0; itVertex < svNVertex; itVertex++ ) {
		reco::Vertex vertex = svTagInfo.secondaryVertex( itVertex );
		float pt = vertex.p4().pt();
		int nTracksSV = vertex.nTracks();

		//Take the vertex with the most tracks, tie breaker is the sum pt of vertex
		if( ( nTracksSV > mostTracks ) || ( nTracksSV == mostTracks && pt > tieBreaker ) ) {
			mostTracks = nTracksSV;
			tieBreaker = pt;
			svIndex = itVertex;
		}
	}

	//Pick the Selected Vertex
	const reco::Vertex & selVertex = svTagInfo.secondaryVertex( svIndex );

	//Set the gloabl vertex to the selected vertex
	selSV	= selVertex;

	//Store the Quantities of the Selected Vertex
	
	svX			= selVertex.x();
	svY			= selVertex.y();
	svZ			= selVertex.z();

	svNTracks	= selSV.nTracks();
	svChi2		= selSV.chi2();
	scNDof		= selSV.ndof();

	svXError	= selSV.xError();
	svYError	= selSV.yError();
	svZError	= selSV.zError();

	float pvxE	= selPV.xError();
	float pvyE	= selPV.yError();
	float pvzE	= selPV.zError();

	float xE	= metric2D( svXError, pvxE );
	float yE	= metric2D( svYError, pvyE );
	float zE	= metric2D( svZError, pvzE );

	float dx	= selPV.x() - svX;
	float dy 	= selPV.y() - svY;
	float dz	= selPV.z() - svZ;

	TVector3 pvVector3D( selPV.x(), selPV.y(), selPV.z() );
	TVector3 pvVector2D( selPV.x(), selPV.y(), 0.0 );
	TVector3 svVector3D( selSV.x(), selSV.y(), selSV.z() );
	TVector3 svVector2D( selSV.x(), selSV.y(), 0.0 );

	TVector3 svMom3D( selVertex.p4().x(), selVertex.p4().y(), selVertex.p4().z() );
	TVector3 svMom2D( selVertex.p4().x(), selVertex.p4().y(), 0.0 );

	float sign2D = ( svMom2D * ( svVector2D - pvVector2D ) ) > 0 ? -1 : 1;
	float sign3D = ( svMom3D * ( svVector3D - pvVector3D ) ) > 0 ? -1 : 1;

	TVector3 pvToVertex2D( sign2D * dx, sign2D * dy, 0.0 );
	TVector3 pvToVertex3D( sign3D * dx, sign3D * dy, sign3D * dz );

	svAngle2D = pvToVertex2D.Angle( svMom2D );
	svAngle3D = pvToVertex3D.Angle( svMom3D );

	svMass 		= selVertex.p4().mass();
	svPt 		= selVertex.p4().pt();
	svEta		= selVertex.p4().eta();
	svPhi		= selVertex.p4().phi();
	svLxySig	= std::sqrt( dx * dx + dy * dy ) / std::sqrt( xE * xE + yE * yE );
	svLxyzSig	= std::sqrt( dx * dx + dy * dy + dz * dz ) / std::sqrt( xE * xE + yE * yE + zE * zE );
	svLxy		= std::sqrt( dx * dx + dy * dy );
	svLxyz		= std::sqrt( dx * dx + dy * dy + dz * dz );

}

//Function to add the primary vertex information
void DisplacedJet::addV0Info( const reco::TrackRefVector tracks ) {

	//Get the transient track builder
	edm::ESHandle<TransientTrackBuilder> builder;
	iSetup.get<TransientTrackRecord>().get( "TransientTrackBuilder", builder );

	//Loop over all the pairs of tracks
	int nTracks = displacedTracks.size();
	for( int itTrack1 = 0; itTrack1 < nTracks; ++itTrack1 ) {
		for( int itTrack2 = 0; itTrack2 < nTracks; ++ itTrack2 ) {

			DisplacedTrack track1 = displacedTracks[itTrack1];
			DisplacedTrack track2 = displacedTracks[itTrack2];

			if( !track1.isValid || !track2.isValid ) continue;

			if( std::fabs( track1.ip2dSig ) < 2 || std::fabs( track2.ip2dSig ) < 2 ) continue;

			Displaced2TrackVertex vertex( track1, track2, selPV, iSetup, debug );

			if( vertex.chi2 > 20 || !vertex.isValid ) continue;

			displacedV0Vector.push_back( vertex );		
		}
	}

	int v0VectorSize = displacedV0Vector.size();

	for( int itVertex = 0; itVertex < v0VectorSize; ++itVertex ) {

		Displaced2TrackVertex vertex = displacedV0Vector[itVertex];

		static const float MIN_DISTANCE = 0.05;
		bool oneTrackNuclear = ( vertex.dInnerHitTrack1 < MIN_DISTANCE ) || ( vertex.dInnerHitTrack2 < MIN_DISTANCE );
		bool twoTrackNuclear = ( vertex.dInnerHitTrack1 < MIN_DISTANCE ) && ( vertex.dInnerHitTrack2 < MIN_DISTANCE );

		//Make sure that the inner hits are in the pixel barrel and the vertex is transversely close to a given layer
		//Josh always assumes the vertex is not near BPIX3 - maybe ask why?
		bool vertexNearBeamPipe = vertex.lxy < 2.3 	&& vertex.lxy > 2.1 && vertex.tot_xyE < 0.03;
		bool vertexNearBPIX1 	= vertex.lxy < 5 	&& vertex.lxy > 4	&& vertex.tot_xyE < 0.03;
		bool vertexNearBPIX2	= vertex.lxy < 7.8	&& vertex.lxy > 6.8	&& vertex.tot_xyE < 0.03;
		bool vertexNearBPIX		= vertexNearBeamPipe || vertexNearBPIX1 || vertexNearBPIX2;

		if( oneTrackNuclear ) jetOneTrackNuclearCount++;
		if( twoTrackNuclear ) jetTwoTrackNuclearCount++;
		if( vertexNearBPIX1 ) jetVertexNearBPIX1++;
		if( vertexNearBPIX2 ) jetVertexNearBPIX2++;
		if( vertexNearBPIX  ) jetVertexNearBPIX++;
		if( vertexNearBPIX && oneTrackNuclear ) jetLooseNuclear++;
		if( vertexNearBPIX && twoTrackNuclear ) jetTightNuclear++;

		//Restrict the inner hit behind calls to tracks that can enter the Median IP significance calculation.
		if( vertex.innerHitBehindVertex1 && vertex.innerHitBehindVertex2 ) jetNV0HitBehindVertex++;

		//Count all vertices (including the ones from tracks who do not enter the Median IP significance calculation)
		//Require the vertices aren't fake with hits behind the position
		//Must have 2D significance of at least 1 as to not be consistent with the beamspot
		//Must be outside 1/2 mm

		if( !vertex.innerHitBehindVertex1 && !vertex.innerHitBehindVertex2 && vertex.chi2 < 20 && vertex.lxy > 0.05 && std::fabs( vertex.lxySig ) > 3 ) {
			jetNV0NoHitBehindVertex++;
			displacedV0VectorCleaned.push_back( vertex );
		}

		//Particle Comparisons
		jetNV0KShort += vertex.isKShort ? 1 : 0;
		jetNV0Lambda += vertex.isLambda ? 1 : 0;
	}

	//Create the graph of vertex associatiation and find the h index (find VertexCliques();
	//Fill information for the clustering of the vertices;
	calcClusterSize( displacedV0VectorCleaned, 2 );
}

void DisplacedJet::calcClusterSize( const DisplacedV0Collection& vertices, const float& errorWindow ) {

	const int nVtx = vertices.size();

	if( nVtx == 0 ) return;

	DisplacedCluster *maxCluster = NULL;

	//Check every vertex in a jet for the largest number of vertices within an error window
	//Call the vertexing behind the checked center
	
	for( int center  = 0; center < nVtx; ++center ) {
		
		Displaced2TrackVertex centerVtx = vertices[center];
		if( !centerVtx.isValid ) continue;

		//Candidate cluster
		DisplacedCluster cluster_temp( centerVtx, selPV, iSetup, debug );

		//Loop over all the possible neighbors

		for( int neighbor = 0; neighbor < nVtx; ++neighbor ) {
		
			if( center == neighbor )  continue;
			Displaced2TrackVertex neighVtx = vertices[neighbor];

			if( !neighVtx.isValid ) continue;

			float dx		= centerVtx.x - neighVtx.x;
			float dy		= centerVtx.y - neighVtx.y;
			float distance	= metric2D( dx, dy );
			float sig		= distance / metric2D( neighVtx.tot_xyE, centerVtx.tot_xyE );

			if( sig < errorWindow ) cluster_temp.addVertex( neighVtx );
		}

		//Set the max cluster for the first iteration
		
		if( maxCluster == NULL ) {
			maxCluster = &cluster_temp;
			continue;
		}

		// Otherwise, it needs to have more vertices

		if( cluster_temp.nV0 > maxCluster->nV0 ) maxCluster = &cluster_temp;
	}

	//Build the related cluster quantities
	maxCluster->buildClusterQuantities();

	//Set quantities related to the max cluster
	
	v0Cluster				= maxCluster;
	jetV0ClusterSize		= maxCluster->nV0;
	jetV0ClusterLxy			= maxCluster->center.lxy;
	jetV0ClusterLxyz		= maxCluster->center.lxyz;
	jetV0ClusterLxySig		= maxCluster->center.lxySig;
	jetV0ClusterLxyzSig		= maxCluster->center.lxyzSig;
	jetV0ClusterX			= maxCluster->center.x;
	jetV0ClusterY			= maxCluster->center.y;
	jetV0ClusterZ			= maxCluster->center.z;
	jetV0ClusterChi2		= maxCluster->fitChi2;
	jetV0ClusterIntercept	= maxCluster->interceptPv;
	jetV0ClusterAngle		= maxCluster->cosAngleToFit;
	jetV0ClusterAngleMom	= maxCluster->cosAngleToMomentum;
	jetV0ClusterNTracks		= maxCluster->nTracks;
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

//The following function adds info for the intermediate vertex 
void DisplacedJet::addIVFCollection( const reco::VertexCollection& vertices ) {
	JetVertexAssociation JVAIVF ("IVF", selPV, debug );
	for( reco::VertexCollection::const_iterator itVertex = vertices.begin(); itVertex != vertices.end(); ++itVertex ) {
		JVAIVF.addVertex( *itVertex );
	}

	const std::pair<reco::Vertex, float> 	bestVertexPair 	= JVAIVF.getBestVertex( jet, "oneOverR" );
	const reco::Vertex						bestVertex		= bestVertexPair.first;
	const float								bestVertexScore = bestVertexPair.second;

	selIVF = bestVertex;

	float x = bestVertex.x();
	float y = bestVertex.y();
	float z = bestVertex.z();
	float dx = x - selPV.x();
	float dy = y - selPV.y();
	float dz = z - selPV.z();

	selIVFIsPVScore = std::sqrt( (dx/x) * (dx/x) + (dy/y) * (dy/y) + (dz/z) * (dz/z) );
	selIVFIsPV		= selIVFIsPVScore < .05;

	float svxE = bestVertex.xError();
	float svyE = bestVertex.yError();
	float svzE = bestVertex.zError();
	float pvxE = selPV.xError();
	float pvyE = selPV.yError();
	float pvzE = selPV.zError();
	float xE   = std::sqrt( svxE * svxE + pvxE * pvxE );
	float yE   = std::sqrt( svyE * svyE + pvyE * pvyE );
	float zE   = std::sqrt( svzE * svzE + pvzE * pvzE );

	ivfX		= selIVFIsPV ? 0 : x;
	ivfY		= selIVFIsPV ? 0 : y;
	ivfZ		= selIVFIsPV ? 0 : z;
	ivfXError	= selIVFIsPV ? 0 : svxE;
	ivfYError	= selIVFIsPV ? 0 : svyE;
	ivfZError	= selIVFIsPV ? 0 : svzE;

	ivfNTracks 	= selIVFIsPV ? 0 : bestVertex.nTracks();
	ivfMass		= selIVFIsPV ? 0 : bestVertex.p4().mass();
	ivfLxySig	= selIVFIsPV ? 0 : std::sqrt( dx * dx + dy * dy ) / std::sqrt( xE * xE + yE * yE );
	ivfLxyzSig	= selIVFIsPV ? 0 : std::sqrt( dx * dx + dy * dy + dz * dz ) / std::sqrt( xE * xE + yE * yE + zE * zE );
	ivfLxy		= selIVFIsPV ? 0 : std::sqrt( dx * dx + dy * dy );
	ivfLxyz		= selIVFIsPV ? 0 : std::sqrt( dx * dx + dy * dy + dz * dz );

	ivfMatchingScore = bestVertexScore;
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

