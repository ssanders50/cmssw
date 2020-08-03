#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "Math/Vector3D.h"

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/ESWatcher.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"

#include "DataFormats/HeavyIonEvent/interface/Centrality.h"

#include "DataFormats/HeavyIonEvent/interface/EvtPlane.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Track/interface/CoreSimTrack.h"
#include "SimDataFormats/EncodedEventId/interface/EncodedEventId.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/TrackingHit/interface/UpdatablePSimHit.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CondFormats/DataRecord/interface/HeavyIonRPRcd.h"
#include "CondFormats/DataRecord/interface/HeavyIonRcd.h"
#include "CondFormats/HIObjects/interface/CentralityTable.h"
#include "CondCore/DBOutputService/interface/PoolDBOutputService.h"
#include "CondFormats/HIObjects/interface/RPFlatParams.h"

#include "RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneFlatten.h"
#include "TList.h"
#include "TString.h"
#include <ctime>
#include <cstdlib>

#include "RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneList.h"
#include "RecoHI/HiEvtPlaneAlgos/interface/LoadEPDB.h"

using namespace std;
using namespace hi;

#include <vector>
using std::vector;


//
// class declaration
//

class HiEvtPlaneFlatProducer : public edm::stream::EDProducer<> {
public:
  explicit HiEvtPlaneFlatProducer(const edm::ParameterSet&);
  ~HiEvtPlaneFlatProducer() override;

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  // ----------member data ---------------------------

  std::string centralityVariable_;
  std::string centralityLabel_;
  std::string centralityMC_;

  edm::InputTag centralityBinTag_;
  edm::EDGetTokenT<int> centralityBinToken;

  edm::InputTag centralityTag_;
  edm::EDGetTokenT<reco::Centrality> centralityToken;

  edm::InputTag vertexTag_;
  edm::EDGetTokenT<std::vector<reco::Vertex>> vertexToken;

  edm::InputTag inputPlanesTag_;
  edm::EDGetTokenT<reco::EvtPlaneCollection> inputPlanesToken;

  edm::InputTag trackTag_;
  edm::EDGetTokenT<reco::TrackCollection> trackToken;
  edm::Handle<reco::TrackCollection> trackCollection_;

  edm::ESWatcher<HeavyIonRcd> hiWatcher;
  edm::ESWatcher<HeavyIonRPRcd> hirpWatcher;

  const int FlatOrder_;
  int NumFlatBins_;
  double minvtx_;
  double maxvtx_;
  int flatnvtxbins_;
  double flatminvtx_;
  double flatdelvtx_;
  double caloCentRef_;
  double caloCentRefWidth_;
  int CentBinCompression_;
  HiEvtPlaneFlatten * flat[NumEPNames];
  bool useOffsetPsi_;
  double nCentBins_;
  int pcnt;
};
//
// constants, enums and typedefs
//

typedef std::vector<TrackingParticle>                   TrackingParticleCollection;
typedef TrackingParticleRefVector::iterator               tp_iterator;


//
// static data member definitions
//

//
// constructors and destructor
//
HiEvtPlaneFlatProducer::HiEvtPlaneFlatProducer(const edm::ParameterSet& iConfig):
  centralityVariable_ ( iConfig.getParameter<std::string>("centralityVariable")),
  centralityBinTag_ ( iConfig.getParameter<edm::InputTag>("centralityBinTag")),
  centralityTag_ ( iConfig.getParameter<edm::InputTag>("centralityTag")),
  vertexTag_  ( iConfig.getParameter<edm::InputTag>("vertexTag")),
  inputPlanesTag_ ( iConfig.getParameter<edm::InputTag>("inputPlanesTag")),
  trackTag_ ( iConfig.getParameter<edm::InputTag>("trackTag")),
  FlatOrder_ ( iConfig.getParameter<int>("FlatOrder")),
  NumFlatBins_ ( iConfig.getParameter<int>("NumFlatBins")),
  minvtx_ ( iConfig.getParameter<double>("minvtx") ),
  maxvtx_ ( iConfig.getParameter<double>("maxvtx") ),
  flatnvtxbins_ ( iConfig.getParameter<int>("flatnvtxbins") ),
  flatminvtx_ ( iConfig.getParameter<double>("flatminvtx") ),
  flatdelvtx_ ( iConfig.getParameter<double>("flatdelvtx") ),
  caloCentRef_ ( iConfig.getParameter<double>("caloCentRef")),
  caloCentRefWidth_ ( iConfig.getParameter<double>("caloCentRefWidth")),
  CentBinCompression_ ( iConfig.getParameter<int>("CentBinCompression")),
  useOffsetPsi_ ( iConfig.getParameter<bool>("useOffsetPsi"))
{
//  UseEtHF = kFALSE;
  nCentBins_ = 200.;

  if(iConfig.exists("nonDefaultGlauberModel")){
    centralityMC_ = iConfig.getParameter<std::string>("nonDefaultGlauberModel");
  }
  centralityLabel_ = centralityVariable_+centralityMC_;

  centralityBinToken = consumes<int>(centralityBinTag_);

  centralityToken = consumes<reco::Centrality>(centralityTag_);

  vertexToken = consumes<std::vector<reco::Vertex>>(vertexTag_);

  trackToken = consumes<reco::TrackCollection>(trackTag_);

  inputPlanesToken = consumes<reco::EvtPlaneCollection>(inputPlanesTag_);

  //register your products
  produces<reco::EvtPlaneCollection>();
   //now do what ever other initialization is needed
  for(int i = 0; i<NumEPNames; i++) {
    flat[i] = new HiEvtPlaneFlatten();
    flat[i]->init(FlatOrder_,NumFlatBins_,flatnvtxbins_,flatminvtx_,flatdelvtx_,EPNames[i],EPOrder[i]);
  }
  pcnt=0;

}


HiEvtPlaneFlatProducer::~HiEvtPlaneFlatProducer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  for(int i = 0; i<NumEPNames; i++) {
    delete flat[i];
  }

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
HiEvtPlaneFlatProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  using namespace reco;



  if( hiWatcher.check(iSetup) ) {
    //
    //Get Size of Centrality Table
    //
    edm::ESHandle<CentralityTable> centDB_;
    iSetup.get<HeavyIonRcd>().get(centralityLabel_,centDB_);
    nCentBins_ = centDB_->m_table.size();
    for(int i = 0; i<NumEPNames; i++) {
      if(caloCentRef_>0) {
	int minbin = (caloCentRef_-caloCentRefWidth_/2.)*nCentBins_/100.;
	int maxbin = (caloCentRef_+caloCentRefWidth_/2.)*nCentBins_/100.;
	minbin/=CentBinCompression_;
	maxbin/=CentBinCompression_;
	if(minbin>0 && maxbin>=minbin) {
	  if(EPDet[i]==HF || EPDet[i]==Castor) flat[i]->setCaloCentRefBins(minbin,maxbin);
	}
      }
    }
  }
  //
  //Get flattening parameter file.  
  //
  if( hirpWatcher.check(iSetup) ) {
    edm::ESHandle<RPFlatParams> flatparmsDB_;
    iSetup.get<HeavyIonRPRcd>().get(flatparmsDB_);
    LoadEPDB db(flatparmsDB_,flat);
  } //rp record change
  
  //
  //Get Centrality
  //
  int bin = 0;
  int cbin = 0;
  edm::Handle<int> cbin_;
  iEvent.getByToken(centralityBinToken, cbin_);
  cbin = *cbin_;
  bin = cbin/CentBinCompression_;
  //double centval = cbin*100./(double) nCentBins_;
  //
  //Get Vertex
  //
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vertexToken,vertices);
  
  //best vertex
  double bestvz=-999.9, bestvx=-999.9, bestvy=-999.9;
  const reco::Vertex & vtx = (*vertices)[0];
  bestvz = vtx.z(); bestvx = vtx.x(); bestvy = vtx.y();
  math::XYZPoint bestvtx(bestvx,bestvy,bestvz);
  if(bestvz<minvtx_ || bestvz>maxvtx_) return;

  //
  //Get Event Planes
  //
  
  edm::Handle<reco::EvtPlaneCollection> evtPlanes_;
  iEvent.getByToken(inputPlanesToken,evtPlanes_);
  
  if(!evtPlanes_.isValid()){
    return ;
  }

  auto evtplaneOutput = std::make_unique<EvtPlaneCollection>();
  EvtPlane * ep[NumEPNames];
  for(int i = 0; i<NumEPNames; i++) {
    ep[i]=nullptr;
  }
  int indx = 0;
  for (EvtPlaneCollection::const_iterator rp = evtPlanes_->begin();rp !=evtPlanes_->end(); rp++) {
	double s = rp->sumSin(0);
	double c = rp->sumCos(0);
	uint m = rp->mult();
	double soff = s;
	double coff = c;
	double psiOffset = -10;
	double psiFlat = -10;
	if(rp->angle(0)>-5) {
	  if(useOffsetPsi_) {
	    soff = flat[indx]->getSoffset(s, bestvz, bin);
	    coff = flat[indx]->getCoffset(c, bestvz, bin);
	    psiOffset = flat[indx]->getOffsetPsi(soff, coff);
	  }
	  psiFlat = flat[indx]->getFlatPsi(psiOffset,bestvz,bin);
	}
	ep[indx]= new EvtPlane(indx, 2, psiFlat, soff, coff,rp->sumw(), rp->sumw2(), rp->sumPtOrEt(), rp->sumPtOrEt2(), m);
	ep[indx]->addLevel(0, rp->angle(0), s, c);
	ep[indx]->addLevel(3,0., rp->sumSin(3), rp->sumCos(3));
	if(useOffsetPsi_) ep[indx]->addLevel(1, psiOffset, soff, coff);
	++indx;
    }
  
  for(int i = 0; i< NumEPNames; i++) {
    if(ep[i]!=nullptr) evtplaneOutput->push_back(*ep[i]);
    
  }
  ++pcnt;
  iEvent.put(std::move(evtplaneOutput));
  for(int i = 0; i<indx; i++) delete ep[i];
}


//define this as a plug-in
DEFINE_FWK_MODULE(HiEvtPlaneFlatProducer);
