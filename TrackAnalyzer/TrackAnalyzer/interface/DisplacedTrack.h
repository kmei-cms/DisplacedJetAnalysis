class DisplacedTrack {

    public:

	DisplacedTrack( const reco::TrackRef  &ref,
	                const reco::Vertex    &pv,
					const edm::EventSetup &iSetup,
					const int             &debug_) :
					trackRef(ref),
					selPV(pv),
					track(*ref),
					debug(debug_) {

	    //Set Positiions for Primary Vertex

		pvPos            = selPV.position();

		//Set Kinematics of the Track
		
		pt               = track.pt();
		px               = track.px();
		py               = track.py();
		pz               = track.pz();
		eta              = track.eta();
		phi              = track.phi();

        //Quality Information of the Track 
		
		q                = track.charge();
		qoverp           = track.qoverp();
		qoverpSig        = track.qoverp() / track.qoverpError();
		chi2             = track.chi2();
		nLostHits        = track.numberOfLostHits();
		nValidHits       = track.numberOfValidHits();
		algoInt          = track.algo();

		//Impact Parameters Relative to the Primary Vertex
		
		dxy              = track.dxy(pvPos);
		dz               = track.dz(pvPos);
		dxySig           = track.dxy(pvPos) / track.dxyError();
		dzSig            = track.dz(pvPos)  / track.dzError();

		//Initialize the Transient Track Builder Required for ertexing
		edm::ESHandle<TransientTrackBuilder> builder;
		iSetup.get<TransientTrackRecord>().get( "TransientTrackBuilder",builder );
		transientTrack = builder->build( *ref );

		//Trajector Information for Accessing Hits
		static GetTrackTrajInfo getTrackTrajInfo;
		std::vector<GetTrackTrajInfo::Result> trajInfo = getTrackTrajInfo.analyze( iSetup, track );

		//Now look at the hit information if the hits are valid
		
		if( trajInfo[0].valid && trajInfo[1].valid && trajInfo.back().valid ) {
		    isValid       = true;
			//Trajectory State on Surface
			tsosInnerHit  = trajInfo[0].detTSOS;
			tsosNextHit   = trajInfo[1].detTSOS;
			tsosLastHit   = trajInfo.back().detTSOS;

			//Get the detector layer from the trajectory info
			const DetLayer &detLayerInner = *( trajInfo[0].detLayer );
			
			//Detector Layers
			GeomDetEnumerators::SubDetector subDetLayerInner = detLayerInner.subDetector();

			//Check if the track has pixel hits
			hasPixelHits  = ( subDetLayerInner == GeomDetEnumerators::PixelBarrel || subDetLayerInner == GeomDetEnumerators::PixelEndcap );

			//Positions
			innerPos      = tsosInnerHit.globalPosition();
			nextPos       = tsosNextHit.globalPosition();
			lastPos       = tsosLastHit.globalPosition();

			//Momenta
			innerPosMom   = tsosInnerHit.globalMomentum();
			lastPosMom    = tsosLastHit.globalMomentum();
		
		    TVector3 pvVector3D( selPV.x(), selPV.y(), selPV.z() );
	 	    TVector3 pvVector2D( selPV.x(), selPV.y(), 0.0 );

		    TVector3 innerPos3D( innerPos.x(), innerPos.y(), innerPos.z() );
		    TVector3 innerPos2D( innerPos.x(), innerPos.y(), 0.0 );
		
		    TVector3 lastPos3D ( lastPos.x(), lastPos.y(), lastPos.z() );
		    TVector3 lastPos2D ( lastPos.x(), lastPos.y(), 0.0 );

		    TVector3 innerMom3D  ( innerPosMom.x(), innerPosMom.y(), innerPosMom.z() );
		    TVector3 innerMom2D  ( innerPosMom.x(), innerPosMom.y(), 0.0 );

		    TVector3 outerMom3D( lastPosMom.x(), lastPosMom.y(), lastPosMom.z() );
		    TVector3 outerMom2D( lastPosMom.x(), lastPosMom.y(), 0.0 );
		    
			TVector3 refMom3D( px, py, pz );
		    TVector3 refMom2D( px, py, 0.0 );

			//Angular Variables Related to the Primary Vertex
			angleMomentumAndPVAtInnerHit2D = ( -1 * ( pvVector2D - innerPos2D ) ).Angle( refMom2D );
			angleMomentumAndPVAtInnerHit3D = ( -1 * ( pvVector3D - innerPos3D ) ).Angle( refMom3D );
			angleMomentumAndPVAtOuterHit2D = ( -1 * ( pvVector2D - innerPos2D ) ).Angle( innerMom2D );
			angleMomentumAndPVAtOuterHit3D = ( -1 * ( pvVector3D - innerPos3D ) ).Angle( innerMom3D );
		}

		else {
            isValid = false;
		}

        //Calculate the Impact Parameters
		
		if( transientTrack.isValid() ) {
		
		//Separate the calculations for the transverse and 3D impact parameter measurements
		    std::pair<bool,Measurement1D> ip2dMeasurement = IPTools::absoluteTransverseImpactParameter( transientTrack, selPV );
			std::pair<bool,Measurement1D> ip3dMeasurement = IPTools::absoluteImpactParameter3D( transientTrack, selPV );

			ip2d             = ip2dMeasurement.second.value();
			ip2dSig          = ip2dMeasurement.second.significance();
			ip3d             = ip3dMeasurement.second.value();
			ip3dSig          = ip3dMeasurement.second.significance();

		}

		else {
		    isValid = false;
		}
	}

    //Define variables used in the constructor
	
	const reco::TrackRef     &trackRef;
	const reco::Vertex       &selPV;
	const reco::Track        &track;
	const int                debug;

	//Define all variables used in this class
	
	//Info
	bool isValid;

	//Kinematics
	float pt;
	float px;
	float py;
	float pz;
	float eta;
	float phi;

	float q;
	float qoverp;
	float qoverpSig;
	float chi2;
	int algoInt;

	//Quality Info
	float nLostHits;
	float nValidHits;
	bool  hasPixelHits;

	//Impact Parameter Related
	float dxy;
	float dz;
	float dxySig;
	float dzSig;
	float ip2d;
	float ip3d;
	float ip2dSig;
	float ip3dSig;

    //Track Associated Objects
	reco::Vertex::Point                   pvPos;
	reco::TransientTrack                  transientTrack;
	std::vector<GetTrackTrajInfo::Result> trajInfo;

	//Position of the hits in the track
	GlobalPoint   innerPos;
	GlobalPoint   nextPos;
	GlobalPoint   lastPos;

	GlobalVector  innerPosMom;
	GlobalVector  nextPosMom;
	GlobalVector  lastPosMom;

	//Trajector State of the Inner and Next Hit
	TrajectoryStateOnSurface tsosInnerHit;
	TrajectoryStateOnSurface tsosNextHit;
	TrajectoryStateOnSurface tsosLastHit;

	//Angle Related Variables for the Track
	float angleMomentumAndPVAtInnerHit2D;
	float angleMomentumAndPVAtInnerHit3D;
	float angleMomentumAndPVAtOuterHit2D;
	float angleMomentumAndPVAtOuterHit3D;
	

	
};
