typedef std::vector<DisplacedTrack> DisplacedTrackCollection;
typedef std::vector<Displaced2TrackVertex> DisplacedV0Collection;

class DisplacedJet {
    
	public:
	
	DisplacedJet( const reco::CaloJet   &jet_,
	              const reco::Vertex    &primaryVertex_,
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

		nTracksRegPrompt                = 0;
		nTracksRegDisp                  = 0;
		nTracksRe

	}

    //Track Association Variables
    int nTracks;
	int nTracksPrompt;
	int nTracksDisp;
	int nTracksRegPrompt;
	int nTracksRegDisp;
	int nTracksRegPromptUp;
	int nTracksRegPromptDown;
	int nTracksRegDispUp;
	int nTracksRegDispDown;


}

