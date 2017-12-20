//Container for Displaced2TracksVertices contained in a combinatorial cluster
typedef std::vector<Displaced2TrackVertex> DisplacedV0Collection;

class DisplacedCluster {

	public:
			
		DisplacedCluster( const Displaced2TrackVertex center_,
						  const reco::Vertex& pv,
						  const edm::EventSetup& iSetup,
						  const int& debug_ ) :
						  center( center_ ),
						  selPV( pv ),
						  debug( debug_ ) {
			
			//Recalculated lxy and lxySignificance
			//For information related to the center, access the center directly

			nV0					= 0;
			nTracks				= 0;
			fitChi2				= 0;
			slope				= 0;
			interceptPv			= -999999;
			interceptZero		= -999999;
			cosAngleToFit		= -10;
			cosAngleToMomentum	= -10;
			meanX 				= -9999;
			meanY				= -9999;
			meanZ				= -9999;
			sumPX				= -9999;
			sumPY				= -9999;
			sumPZ				= -9999;
		}

		//Declare Methods
		void addVertex( Displaced2TrackVertex& vtx );
		void buildClusterQuantities();
		void findUniqueTracks();
		void fitVertices();
		bool containsVertex( const Displaced2TrackVertex& vtx );

		//Quantities Related ot the vertex at the center of the cluster
		float fitChi2;
		float slope;
		float interceptPv;
		float interceptZero;
		float cosAngleToFit;
		float cosAngleToMomentum;

		int nV0;
		int nTracks;

		float meanX;
		float meanY;
		float meanZ;
		float sumPX;
		float sumPY;
		float sumPZ;

		//Constructor Initialized Qunatities
		Displaced2TrackVertex	 center;
		const reco::Vertex&		 selPV;
		const int				 debug;

		//Container of Vertices
		std::vector<Displaced2TrackVertex> v0s;
		std::vector<DisplacedTrack> uniq_tracks;
	
};

//Add a single vertex
void DisplacedCluster::addVertex( Displaced2TrackVertex& vtx ) {
	v0s.push_back( vtx );
	nV0++;
}

//Fit the vertices curently contaiend in the cluster
void DisplacedCluster::fitVertices() {
	int nVertices = v0s.size();

	//Define vector pointing to the center of the cluster
	TVector3 pvToCenter( -1 * ( selPV.x() - meanX ), -1 * ( selPV.y() - meanY ), 0.0 );
	TVector3 vertexMomentum( center.px, center.py, 0.0 );

	if( nVertices == 0 ) return;

	cosAngleToMomentum = pvToCenter.Angle( vertexMomentum );

	if( nVertices < 2 ) return;

	float x[nVertices];
	float y[nVertices];
	//float xE[nVertices];
	//float yE[nVertices];

	//Loop over all the vertices and fill the arrays for the TGraph
	//Relative to the primary vertex system
	
	for( int itVertex = 0; itVertex < nVertices; ++itVertex ) {
		x[itVertex] 	= v0s[itVertex].x - meanX + .000001 * itVertex;
		y[itVertex] 	= v0s[itVertex].y - meanY + .000001 * itVertex;
		//xE[itVertex] 	= v0s[itVertex].xE;
		//yE[itVertex]	= v0s[itVertex].yE;
	}

	//Build the Graph
	TGraphErrors graph( nVertices, x, y, 0, 0);
	graph.Fit( "pol1", "Q" );

	//Extract the parameters
	TF1 *fit = graph.GetFunction("pol1");
	if( fit == NULL ) return;
	float p0 = fit->GetParameter(0);
	float p1 = fit->GetParameter(1);

	fitChi2 = fit->GetChisquare();

	//Vector pointing along the fit line to the center of the cluster
	TVector3 fitExtrapPlus( 1, p0+p1, 0 );
	TVector3 fitExtrapMinus( -1, p0 - p1 , 0 );

	TVector3 vertexSumMomentum( sumPX, sumPY, 0.0 );

	//Pick the direction in the direction of the vertex momentum
	if( vertexMomentum * fitExtrapPlus > 0 ) {
		cosAngleToFit 	= pvToCenter.Angle( fitExtrapPlus );
		interceptPv 	= pvToCenter.Mag() * sin( pvToCenter.Angle( fitExtrapPlus ) );
	}

	else if ( vertexMomentum * fitExtrapMinus > 0 ) {
		cosAngleToFit	= pvToCenter.Angle( fitExtrapMinus );
		interceptPv		= pvToCenter.Mag() * sin( pvToCenter.Angle( fitExtrapMinus ) );
	}
	else {
		cosAngleToFit 	= pvToCenter.Angle( fitExtrapMinus ) < pvToCenter.Angle( fitExtrapPlus ) ? pvToCenter.Angle( fitExtrapMinus ) : pvToCenter.Angle( fitExtrapPlus );
		interceptPv 	= pvToCenter.Angle( fitExtrapMinus ) < pvToCenter.Angle( fitExtrapPlus ) ? pvToCenter.Mag() * sin( pvToCenter.Angle( fitExtrapMinus ) ) : pvToCenter.Mag() * sin( pvToCenter.Angle( fitExtrapPlus ) ); 
	}
}

//Since a cluster can contain combinatorial vertices, there may be redundant tracks
//This function finds the unique tracks by reference and stores them in the uniq_tracks vector
//, eliminating redundancy
void DisplacedCluster::findUniqueTracks() {
	if( nV0 == 0 || v0s.size() == 0 ) return;
	
	for( DisplacedV0Collection::const_iterator itVertex = v0s.begin(); itVertex != v0s.end(); ++itVertex ) {
		bool found_match1 = false;
		bool found_match2 = false;

		if( !itVertex->isValid ) continue;

		for( std::vector<DisplacedTrack>::const_iterator itTrack = uniq_tracks.begin(); itTrack != uniq_tracks.end(); ++itTrack ) {
			if( uniq_tracks.size() == 0 ) break;

			//Compare the track references stored in the DisplacedTrack object
			DisplacedTrack tr1 = itVertex->track1;
			DisplacedTrack tr2 = itVertex->track2;

			//Compare by reference
			bool refMatch1 	= ( itTrack->trackRef == tr1.trackRef );
			bool refMatch2	= ( itTrack->trackRef == tr2.trackRef );

			//Update the matching
			found_match1 	= ( found_match1 || refMatch1 );
			found_match2	= ( found_match2 || refMatch2 );

			//If both are found, we are done
			if( found_match1 && found_match2 ) continue;
		}

		DisplacedTrack tr1 = itVertex->getTrack1();
		DisplacedTrack tr2 = itVertex->getTrack2();

		if( !found_match1 ) uniq_tracks.push_back( tr1 );
		if( !found_match2 ) uniq_tracks.push_back( tr2 );
	}

	nTracks = uniq_tracks.size();
}

//Check if a vertex is already contained in a cluster
bool DisplacedCluster::containsVertex( const Displaced2TrackVertex& vtx ) {
	if( nV0 == 0 ) return false;
	bool found_match = false;

	for( DisplacedV0Collection::const_iterator itVertex = v0s.begin(); itVertex != v0s.end(); ++itVertex ) {
		bool refMatch1 = vtx.track1.trackRef == itVertex->track1.trackRef;
		bool refMatch2 = vtx.track2.trackRef == itVertex->track2.trackRef;
		found_match = ( found_match || (refMatch1 && refMatch2 ) );
	}

	return found_match;
}


//Call to calculate all the relevant quantities
//Calls the necessary methods in the appropriate order

void DisplacedCluster::buildClusterQuantities() {
	if( nV0 == 0 || v0s.size() == 0 ) return;

	meanX = 0; 
	meanY = 0;
	meanZ = 0;

	for( DisplacedV0Collection::const_iterator itVertex = v0s.begin(); itVertex != v0s.end(); ++itVertex ) {
		meanX += itVertex->x / float(nV0);
		meanY += itVertex->y / float(nV0);
		meanZ += itVertex->z / float(nV0);

		sumPX += itVertex->px;
		sumPY += itVertex->py;
		sumPZ += itVertex->pz;
	}

	fitVertices();
}
