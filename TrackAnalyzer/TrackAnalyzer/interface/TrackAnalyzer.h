class TrackAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit TrackAnalyzer(const edm::ParameterSet&);
      ~TrackAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
	  virtual std::vector<std::pair<std::string,int>> createSkimmedTriggerResults(const edm::Event&, edm::Handle<edm::TriggerResults> triggerResultsObject);

      // ----------member data ---------------------------
	  
	  //Variables needed for output file manipulation 
	  
	  Int_t        debugger_ = 0;
	  std::string  outputFileName_;
	  std::string  jetTreeName_;
	  std::string  trackTreeName_;
	  std::string  vertexTreeName_;
	  std::string  genTreeName_;
	  
	  TFile *outputFile_;
	  
      //Define all user trees here 
	  TTree *runTree_; //This is the tree with information like run, lumi, event
	  TTree *jetTree_;
	  TTree *trackTree_;
	  TTree *vertexTree_;
	  TTree *genTree_;

	  //Define all variables needed for branch definition

      //For Run Statistics
      Int_t run = -1;
	  Int_t lumi = -1;
	  Int_t event = -1;
      std::vector<std::pair<std::string,Int_t>> skimmedTriggerResults;

      bool isMC_;


	  //Tokens necessary for the analysis
	  edm::EDGetTokenT<reco::TrackCollection> trackCollectionTag_; //For the generalTracks collection
	  edm::EDGetTokenT<edm::TriggerResults> triggerResultsTag_; //For the HLT trigger results
	  edm::EDGetTokenT<reco::CaloJetCollection> caloJetCollectionTag_; //For the caloJets collection
      edm::EDGetTokenT<edm::View<reco::GenParticle>> genParticleCollectionTag_; //For the generation information
	  edm::EDGetTokenT<reco::VertexCollection> vertexCollectionTag_; //For all of the primary vertices
	  edm::EDGetTokenT<edm::View<reco::Vertex>> secondaryVertexCollectionTag_; //For all of the secondary vertices

	  std::vector<std::string> triggerPaths_;

	  TH1D* histo_tracks_pT;
	  TH1D* histo_caloJets_pT;
	  TH1D* histo_gen_pT;

	  TH1D* histo_pv_z;
	  TH1D* histo_sv_z;
	  TH1F* histo_alphaMax;

      int ipdf = 1;
};

namespace LHAPDF {
    void     initPDFSet(int nset, const std::string& filename, int member=0);
	int      numberPDF(int nset);
	void     usePDFMember(int nset, int member);
	double   xfx(int nset, double x, double Q, int fl);
	double   getXmin(int nset, int member);
	double   getXmax(int nset, int member);
	double   getQ2min(int nset, int member);
	double   getQ2max(int nset, int member);
	void     extrapolate(bool extrapolate=true);
}
