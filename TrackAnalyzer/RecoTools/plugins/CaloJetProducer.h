#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "RecoTools/interface/CaloJetProducer.h"

template <typename T,typename G>
class CaloJetProducer : public edm::EDProducer {
   public:
  explicit CaloJetProducer(const edm::ParameterSet& iConfig):
    src_(consumes<std::vector<T> >(iConfig.getParameter<edm::InputTag>("src"))) 
    {
      produces<std::vector<T> >();

    }
  ~CaloJetProducer() 
    {

    }
    
   protected:
    
     virtual void produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
       using namespace edm;
       
       std::auto_ptr<std::vector<T> > out(new std::vector<T> );
       Handle<std::vector<T> > srcH;
       
       if(iEvent.getByToken(src_,srcH)) 
	 for(unsigned int i=0;i<srcH->size();++i) {
	   T object = srcH->at(i);
	   out->push_back(object);
	 }
       iEvent.put(out);
     } 


      edm::EDGetTokenT<std::vector<T> > src_;           //input Collection
};


#endif
