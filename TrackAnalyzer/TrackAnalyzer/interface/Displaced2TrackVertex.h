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

		//Build the object with two tracks (maybe use pair instead of vector?)
		
		std::vector<reco::TransientTrack> trackPair;
		trackPair.push_back(track1.transientTrack);
		trackPair.push_back(track2.transientTrack);

		//Fit the pair and check validity
		KalmanVertexFitter fitter(true); //Refits the vertex
		if( track1.isValid && track2.isValid ) {
			tVertex = fitter.vertex( trackPair );
			isValid = tVertex.isValid();
		}
		else {
			isValid = false;
		}

		//Make sure the fit itself is actually valid
		if( isValid ) {
			vertex 	= tVertex; //Convert transient vertex to reco::Vertex
			chi2 	= vertex.chi2();

			//Add Kinematics to Object
			mass	= vertex.p4().mass();
			charge	= track1.q + track2.q;
			pt		= vertex.p4().pt();
			eta 	= vertex.p4().eta();
			phi		= vertex.p4().phi();
			px		= vertex.p4().px();
			py 		= vertex.p4().py();
			pz		= vertex.p4().pz();
			dr		= reco::deltaR( track1.eta, track1.phi, track2.eta, track2.phi );

			//Calculate the lambda hypothesis mass?
			ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> vec1;
			ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> vec2;

			if( track1.pt > track2.pt && charge == 0 ) {
				vec1.SetM( 0.938272 ); //Proton Mass
				vec2.SetM( 0.139570 ); //Pion Mass
			}

			else if ( track1.pt < track2.pt && charge == 0 ) {
				vec1.SetM( 0.139570 ); //Pion Mass
				vec2.SetM( 0.938272 ); //Proton Mass
			}
			else {
				vec1.SetM( 0.0 );
				vec2.SetM( 0.0 );
			}

			//Set the momentum from the tracks
			vec1.SetPx( track1.px );
			vec1.SetPy( track1.py );
			vec1.SetPz( track1.pz );

			vec2.SetPx( track2.px );
			vec2.SetPy( track2.py );
			vec2.SetPz( track2.pz );

			//Build the Mass Candidate
			math::XYZTLorentzVectorD sumLambda;
			sumLambda += vec1;
			sumLambda += vec2;
			massLambda = sumLambda.mass();

			//Position Calculations
			vtxPos		= vertex.position();
			x			= vtxPos.x();
			y			= vtxPos.y();
			z			= vtxPos.z();
			xE			= vertex.xError();
			yE			= vertex.yError();
			zE			= vertex.zError();

			//Distance from the selected primary vertex
			dx			= vtxPos.x() - selPV.x();
			dy			= vtxPos.y() - selPV.y();
			dz			= vtxPos.z() - selPV.z();

			//2D and 3D distances from the primary vertex
			lxy			= metric2D( dx, dy );
			lxyz		= metric3D( dx, dy, dz );

			//Total Error
			tot_xE		= metric2D( xE, selPV.xError() );
			tot_yE		= metric2D( yE, selPV.yError() );
			tot_zE		= metric2D( zE, selPV.zError() );

			//2D and 3D Error
			tot_xyE		= metric2D( tot_xE, tot_yE );
			tot_xyzE	= metric3D( tot_xE, tot_yE, tot_zE );

			//Significances calculated using the errors
			lxySig		= lxy / tot_xyE;
			lxyzSig		= lxyz / tot_xyzE;

			//Make the vectors pointing from the vertex to the inner and next hit
			TVector3	vertex_to_innerHit1( track1.innerX - x, track1.innerY - y, track1.innerZ - z );
			TVector3	vertex_to_innerHit2( track2.innerX - x, track2.innerY - y, track2.innerZ - z );
			TVector3	vertex_to_nextHit1( track1.nextX - x, track1.nextY - y, track1.nextZ - z );
			TVector3	vertex_to_nextHit2( track2.nextX - x, track2.nextY - y, track2.nextZ - z );

			//Inner product should be positive if vertex is real
			innerHitBehindVertex1	= ( vertex_to_nextHit1 * vertex_to_innerHit1 ) < 0;
			innerHitBehindVertex2	= ( vertex_to_nextHit2 * vertex_to_innerHit2 ) < 0;

			//Track inner hit should be near the vertex (reduces chances that it is a nuclear interaction)
			dInnerHitTrack1 = metric3D( x - track1.innerX, y - track1.innerY, z - track1.innerZ);
			dInnerHitTrack2 = metric3D( x - track2.innerX, y - track2.innerY, z - track2.innerZ);
			dInnerHitTrack1Track2 = metric3D ( track1.innerX - track2.innerX, track1.innerY - track2.innerY, track1.innerZ - track2.innerZ );

			//Total Missing Hits
			sumLostHits		= track1.nLostHits + track2.nLostHits;
			sumValidHits	= track1.nValidHits + track2.nValidHits;

			//Check KShort and Lambda Window
			isKShort = ( std::fabs( mass - .493 ) < 0.02 ) && charge == 0 && track1.pt > 1.0 && track2.pt > 1.0 && lxySig > 50.0;
			isLambda = ( std::fabs( massLambda - 1.1156 ) < 0.02 ) && charge == 0 && track1.pt > 1.0 && track2.pt > 1.0 && lxySig > 50.0;
		}
    }

	//Initialize variables used in the constructor
	DisplacedTrack track1;
	DisplacedTrack track2;
	const reco::Vertex selPV;
	int debug;	

	//Associated Objects
	reco::Vertex	vertex;
	TransientVertex tVertex;
	math::XYZPoint	vtxPos;

	//Quality
	bool 	isValid;
	float 	chi2;
	int		sumLostHits;
	int		sumValidHits;

	//Kinematics
	int 	charge;
	float	pt;
	float	eta;
	float	phi;
	float	px;
	float	py;	
	float	pz;
	float	mass;
	float	massLambda;
	float	dr;

	//Position
	float	x;
	float	y;
	float	z;
	float	xE;
	float	yE;
	float	zE;
	float 	dx;
	float	dy;
	float	dz;
	float	lxy;
	float	lxyz;
	float	lxySig;
	float	lxyzSig;

	//Errors combiend with the primary vertex error
	float	tot_xE;
	float 	tot_yE;
	float 	tot_zE;
	float	tot_xyE;
	float	tot_xyzE;

	//Relative Position to Tracks
	float	dInnerHitTrack1;
	float	dInnerHitTrack2;
	float	dInnerHitTrack1Track2;

	//Particle consistency
	bool	isKShort;
	bool	isLambda;

	//Nuclear Interaction Events
	
	bool	innerHitBehindVertex1;
	bool	innerHitBehindVertex2;

	//Simplify Quadrature Calculations
	float 	metric2D( float x, float y ) {
		return std::sqrt(x*x + y*y);
	}
	
	float 	metric3D( float x, float y, float z ) {
		return std::sqrt(x*x + y*y + z*z);
	}
	
	DisplacedTrack getTrack1() const { return track1; }
	DisplacedTrack getTrack2() const { return track2; }

};
