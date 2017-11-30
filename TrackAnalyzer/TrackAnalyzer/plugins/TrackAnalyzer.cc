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
#include <vector>
#include <string>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

// included in the tutorial for track analysis
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

// included for more complicated track analysis (vertex matching)
#include "TrackingTools/DetLayers/interface/DetLayer.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/PatternTools/interface/TransverseImpactPointExtrapolator.h"
#include "RecoTracker/DebugTools/interface/GetTrackTrajInfo.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

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

//including all files needed for proper root implementation
#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"

// including all user defined objects
#include "TrackAnalyzer/TrackAnalyzer/interface/DisplacedTrack.h"
#include "TrackAnalyzer/TrackAnalyzer/interface/Displaced2TrackVertex.h"
#include "TrackAnalyzer/TrackAnalyzer/interface/DisplacedJet.h"
#include "TrackAnalyzer/TrackAnalyzer/interface/DisplacedJetEvent.h"
//#include "TrackAnalyzer/RecoTools/plugins/DisplacedJetOverloader.h"
//#include "TrackAnalyzer/TrackAnalyzer/interface/TrackAnalyzerJet.h"
//#include "TrackAnalyzer/TrackAnalyzer/interface/TrackAnalyzerEvent.h"

// including LHAPDF (in TrackAnalyzer.h")
#include "TrackAnalyzer/TrackAnalyzer/interface/TrackAnalyzer.h"


TrackAnalyzer::TrackAnalyzer(const edm::ParameterSet& iConfig) :
   trackCollectionTag_(consumes<edm::View<reco::Track>> (iConfig.getParameter<edm::InputTag>("tracks"))),
   triggerResultsTag_(consumes<edm::TriggerResults> (iConfig.getParameter<edm::InputTag>("triggers"))),
   caloJetCollectionTag_(consumes<edm::View<reco::CaloJet>> (iConfig.getParameter<edm::InputTag>("caloJets"))),
   genParticleCollectionTag_(consumes<edm::View<reco::GenParticle>> (iConfig.getParameter<edm::InputTag>("genParticles"))),
   vertexCollectionTag_(consumes<edm::View<reco::Vertex>> (iConfig.getParameter<edm::InputTag>("primaryVertices"))),
   secondaryVertexCollectionTag_(consumes<edm::View<reco::Vertex>> (iConfig.getParameter<edm::InputTag>("secondaryVertices")))

{
   // initialize the PDFs
   
   // LHAPDF::initPDFSet( ipdf, "NNPDF23_lo_as_0130.qed.LHgrid"); //FIXME : Why does this not work CMSSW_9_2_10?
   
   debugger_       = iConfig.getUntrackedParameter<int>("debugger");
   outputFileName_ = iConfig.getUntrackedParameter<std::string>("outputFileName");
   jetTreeName_    = iConfig.getUntrackedParameter<std::string>("jetTreeName");
   trackTreeName_  = iConfig.getUntrackedParameter<std::string>("trackTreeName");
   vertexTreeName_ = iConfig.getUntrackedParameter<std::string>("vertexTreeName");
   genTreeName_    = iConfig.getUntrackedParameter<std::string>("genTreeName");
   
   //now do what ever initialization is needed
   usesResource("TFileService");
   edm::Service<TFileService> fs;
   histo_tracks_pT = fs->make<TH1D>("pt", "PT", 50, 0, 50);
   histo_caloJets_pT = fs->make<TH1D>("calo_pt","Calo PT", 50, 0, 200); 
   histo_gen_pT = fs->make<TH1D>("gen_pt","Gen PT", 50, 0, 200);

   histo_pv_z = fs->make<TH1D>("pv_z","PV z", 50, -10, 10);
   histo_sv_z = fs->make<TH1D>("sv_z","SV z", 50, -10, 10);
   
   histo_alphaMax = fs->make<TH1F>("alphaMax","alpha", 20, 0, 10);

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
   run   = iEvent.id().run();
   lumi  = iEvent.id().luminosityBlock();
   event = iEvent.id().event();

   std::cout<<"Run: "<<run<<"; LumiBlock: "<<lumi<<"; Event: "<<event<<std::endl;

   //Analysis loop to iterate over the different trigger paths (decide whether the event passes the trigger or not)

   edm::Handle<edm::TriggerResults> triggers;
   iEvent.getByToken(triggerResultsTag_,triggers);

   skimmedTriggerResults.clear();
   skimmedTriggerResults = createSkimmedTriggerResults(iEvent, triggers);
   
   runTree_->Fill();
   
   //Analysis loop to iterate over the different tracks
   edm::Handle<edm::View<reco::Track>> tracks;
   iEvent.getByToken(trackCollectionTag_,tracks);
   
   for( View<reco::Track>::const_iterator itTrack = tracks->begin(); itTrack != tracks->end(); ++itTrack )
   {
       histo_tracks_pT->Fill(itTrack->pt());
   }

   //Analysis loop to iterate over the calojets
  
   edm::Handle<edm::View<reco::CaloJet>> caloJets;
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
   

   float cut_jetPt  = 10.0;
   float cut_jetEta = 4.5;

   DisplacedJetEvent djEvent( isMC_, *(caloJets.product()), *(primaryVertices.product()), cut_jetPt, cut_jetEta, iSetup, debugger_);

#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}

// ------------ method called to fill the trigger results in a vector -------------------
std::vector<std::pair<std::string,int>>
TrackAnalyzer::createSkimmedTriggerResults(const edm::Event& iEvent, edm::Handle<edm::TriggerResults> triggerResultsObject )
{
   const edm::TriggerNames &names = iEvent.triggerNames(*triggerResultsObject);
   std::vector<std::pair<std::string,int>> tempSkimmedTriggerResults;
   tempSkimmedTriggerResults.clear();

   for ( unsigned int itTrigger = 0; itTrigger < triggerResultsObject->size(); ++itTrigger )
   {
       std::string name = names.triggerName(itTrigger);
	   std::string triggerTag = "Displaced";
	   if ( name.find(triggerTag) != std::string::npos )
	   {
	       int accept = int(triggerResultsObject->accept(itTrigger));
		   tempSkimmedTriggerResults.push_back(std::make_pair(name,accept));
	   }
   }

   return tempSkimmedTriggerResults;
}



// ------------ method called once each job just before starting event loop  ------------
void 
TrackAnalyzer::beginJob()
{
    outputFile_ = new TFile(outputFileName_.c_str(), "RECREATE");
	runTree_    = new TTree("RunStats", "General run statistics");

	//Define all the branches for the trees you will be using #FIXME - move each tree initialization to a separate function.
	
    runTree_->Branch("run", &run, "run/I");
	runTree_->Branch("lumi", &lumi, "lumi/I");
	runTree_->Branch("event", &event, "event/I");
	runTree_->Branch("skimTrigResults", &skimmedTriggerResults);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TrackAnalyzer::endJob() 
{
    
	outputFile_->cd();
    //histo_tracks_pT->Write();
	runTree_->Write();
	outputFile_->Close();
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
