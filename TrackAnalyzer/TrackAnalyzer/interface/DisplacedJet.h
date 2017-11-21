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
	    if( debug ) std::cout<<"[DEBUG] In the DisplacedJet code and creating a dispalced jet with pT="<<jet_.pt()<<", eta="<<jet_.eta()<<", phi="<<jet_.phi()<<std::endl;

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


    //Jet Info Extraction Functions
	DisplacedTrackCollection getDisplacedTracks()     { return displacedTracks; }
	reco::TrackCollection    getVertexMatchedTracks() { return vertexMatchedTracks; }
	reco::Vertex             getIVFVertexSelected()   { return selIVF; }
	reco::Vertex             getSVVertex()            { return selSV; }

    //Define functions for this class
	void calcJetAlpha( const reco::TrackCollection&,
	                   const reco::VertexCollection& );
	void addVertexTrackInfo( const reco::TrackRefVector& );

    DisplacedTrackCollection displacedTracks;

	private:

    reco::TrackCollection vertexMatchedTracks;
	reco::TrackRefVector  vertexMatchedTrackRefs;
	reco::Vertex          selIVF;
	reco::Vertex          selSV;

};

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

//Function to calculate the jet variable alpha, defined by Josh as the ratio of (Vertex Tracks Sum pT Matching the Jet / (General Tracks Sum Pt Matching the Jet)
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
}
