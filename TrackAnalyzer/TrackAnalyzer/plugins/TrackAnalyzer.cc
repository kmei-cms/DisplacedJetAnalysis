// -*- C++ -*-
//
// Package:    TrackAnalyzer/TrackAnalyzer
// Class:      TrackAnalyzer
// 
/**\class TrackAnalyzer TrackAnalyzer.cc TrackAnalyzer/TrackAnalyzer/plugins/TrackAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Kelvin Mei
//         Created:  Mon, 06 Nov 2017 22:06:13 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

// included in the tutorial for track analysis
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// included for the ability to read out triggers
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"

// included to get the ak4CaloJet objects
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/JetTracksAssociation.h"

// included to get the generator information
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"

// included to get the generator information with respects to pileup 
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

// included to get primary vertex information
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/BTauReco/interface/VertexTypes.h"

// included to get secondary vertex information
#include "DataFormats/BTauReco/interface/SecondaryVertexTagInfo.h"
#include "DataFormats/BTauReco/interface/TrackIPTagInfo.h"
#include "DataFormats/BTauReco/interface/TrackIPData.h"
#include "DataFormats/BTauReco/interface/CandSecondaryVertexTagInfo.h"
#include "DataFormats/BTauReco/interface/TaggingVariable.h"

// including all user defined objects
//#include "TrackAnalyzer/RecoTools/plugins/DisplacedJetOverloader.h"
//#include "TrackAnalyzer/TrackAnalyzer/interface/TrackAnalyzerJet.h"
//#include "TrackAnalyzer/TrackAnalyzer/interface/TrackAnalyzerEvent.h"


//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class TrackAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit TrackAnalyzer(const edm::ParameterSet&);
      ~TrackAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------


	  //Tokens necessary for the analysis
	  edm::EDGetTokenT<edm::View<reco::Track>> trackCollectionTag_; //For the generalTracks collection
	  edm::EDGetTokenT<edm::TriggerResults> triggerResultsTag_; //For the HLT trigger results
	  edm::EDGetTokenT<edm::View<reco::CaloJet>> caloJetCollectionTag_; //For the caloJets collection
      edm::EDGetTokenT<edm::View<reco::GenParticle>> genParticleCollectionTag_; //For the generation information
	  edm::EDGetTokenT<edm::View<reco::Vertex>> vertexCollectionTag_; //For all of the primary vertices
	  edm::EDGetTokenT<edm::View<reco::Vertex>> secondaryVertexCollectionTag_; //For all of the secondary vertices


	  std::vector<std::string> triggerPaths_;

	  TH1D* histo_tracks_pT;
	  TH1D* histo_caloJets_pT;
	  TH1D* histo_gen_pT;

	  TH1D* histo_pv_z;
	  TH1D* histo_sv_z;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
TrackAnalyzer::TrackAnalyzer(const edm::ParameterSet& iConfig) :
   trackCollectionTag_(consumes<edm::View<reco::Track>> (iConfig.getParameter<edm::InputTag>("tracks"))),
   triggerResultsTag_(consumes<edm::TriggerResults> (iConfig.getParameter<edm::InputTag>("triggers"))),
   caloJetCollectionTag_(consumes<edm::View<reco::CaloJet>> (iConfig.getParameter<edm::InputTag>("caloJets"))),
   genParticleCollectionTag_(consumes<edm::View<reco::GenParticle>> (iConfig.getParameter<edm::InputTag>("genParticles"))),
   vertexCollectionTag_(consumes<edm::View<reco::Vertex>> (iConfig.getParameter<edm::InputTag>("primaryVertices"))),
   secondaryVertexCollectionTag_(consumes<edm::View<reco::Vertex>> (iConfig.getParameter<edm::InputTag>("secondaryVertices")))

{
   //now do what ever initialization is needed
   usesResource("TFileService");
   edm::Service<TFileService> fs;
   histo_tracks_pT = fs->make<TH1D>("pt", "PT", 50, 0, 50);
   histo_caloJets_pT = fs->make<TH1D>("calo_pt","Calo PT", 50, 0, 200); 
   histo_gen_pT = fs->make<TH1D>("gen_pt","Gen PT", 50, 0, 200);

   histo_pv_z = fs->make<TH1D>("pv_z","PV z", 50, -10, 10);
   histo_sv_z = fs->make<TH1D>("sv_z","SV z", 50, -10, 10);

}


TrackAnalyzer::~TrackAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
TrackAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   //Initialization procedure
   int run   = iEvent.id().run();
   int lumi  = iEvent.id().luminosityBlock();
   int event = iEvent.id().event();

   std::cout<<"Run: "<<run<<"; LumiBlock: "<<lumi<<"; Event: "<<event<<std::endl;
   
   //Analysis loop to iterate over the different tracks
   edm::Handle <View<reco::Track>> tracks;
   iEvent.getByToken(trackCollectionTag_,tracks);
   
   for( View<reco::Track>::const_iterator itTrack = tracks->begin(); itTrack != tracks->end(); ++itTrack )
   {
       histo_tracks_pT->Fill(itTrack->pt());
   }


   //Analysis loop to iterate over the different trigger paths (decide whether the event passes the trigger or not)

   edm::Handle<edm::TriggerResults> triggers;
   iEvent.getByToken(triggerResultsTag_,triggers);

   const edm::TriggerNames &names = iEvent.triggerNames(*triggers);

   for( unsigned int itTrigger = 0; itTrigger < triggers->size(); ++itTrigger )
   {
       std::string name = names.triggerName(itTrigger);
	   std::string triggerTag = "Displaced";
	   if (name.find(triggerTag) != std::string::npos)
	   {
	       bool accept = (triggers->accept(itTrigger));
	       std::cout<<name<<" "<<accept<<std::endl;
	   }
   }

   //Analysis loop to iterate over the calojets
   
   edm::Handle<View<reco::CaloJet>> caloJets;
   iEvent.getByToken(caloJetCollectionTag_,caloJets);

   for( View<reco::CaloJet>::const_iterator itCaloJet = caloJets->begin(); itCaloJet != caloJets->end(); ++itCaloJet )
   {
       histo_caloJets_pT->Fill(itCaloJet->pt());
   }

   //Analysis loop to iterate over the gen particles
   edm::Handle<View<reco::GenParticle>> genParticles;
   iEvent.getByToken(genParticleCollectionTag_,genParticles);

   for( View<reco::GenParticle>::const_iterator itGenParticle = genParticles->begin(); itGenParticle != genParticles->end(); ++itGenParticle )
   {
       histo_gen_pT->Fill(itGenParticle->pt());
   }

   //TrackAnalyzerEvent myTrackAnalyzerEvent(caloJets);
   //fillGenInfo(myTrackAnalyzerEvent,genParticles);

   //Analysis loop to iterate over all the vertices
   edm::Handle<View<reco::Vertex>> primaryVertices;
   iEvent.getByToken(vertexCollectionTag_, primaryVertices);

   edm::Handle<View<reco::Vertex>> secondaryVertices;
   iEvent.getByToken(secondaryVertexCollectionTag_, secondaryVertices);


   for( View<reco::Vertex>::const_iterator itPrimaryVertex = primaryVertices->begin(); itPrimaryVertex != primaryVertices->end(); ++itPrimaryVertex )
   {
       histo_pv_z->Fill(itPrimaryVertex->z());
   }
   
   for( View<reco::Vertex>::const_iterator itSecondaryVertex = secondaryVertices->begin(); itSecondaryVertex != secondaryVertices->end(); ++itSecondaryVertex )
   {
       histo_sv_z->Fill(itSecondaryVertex->z());
   }

#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void 
TrackAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TrackAnalyzer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void TrackAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

// ----------- method takes the generator information into a tree ------------------------
//void TrackAnalyzer::fillGenInfo(TrackAnalyzerEvent &djEvent, const reco::GenParticleCollection &gen) 
//{
//    std::cout<<"Filling generator tree with information."<<std::endl;
//}	







//define this as a plug-in
DEFINE_FWK_MODULE(TrackAnalyzer);
