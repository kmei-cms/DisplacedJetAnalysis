class Displaced2TrackVertex {
   
   public:

    Displaced2TrackVertex( DisplacedTrack        &track1_,
                           DisplacedTrack        &track2_,
						   const reco::Vertex    &pv,
						   const edm::EventSetup &iSetup,
						   const int             &debug_) :
						   track1( track1_ ),
						   track2( track2_ ),
						   selPV ( pv ),
						   debug ( debug_ ) {
	
    }

	//Initialize variables used in the constructor
	DisplacedTrack track1;
	DisplacedTrack track2;
	const reco::Vertex selPV;
	int debug;	

};
