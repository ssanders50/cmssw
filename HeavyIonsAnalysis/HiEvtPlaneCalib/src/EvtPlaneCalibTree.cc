#include <memory>


// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "Math/Vector3D.h"

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "DataFormats/CastorReco/interface/CastorTower.h"

#include "DataFormats/HeavyIonEvent/interface/Centrality.h"

#include "DataFormats/HeavyIonEvent/interface/EvtPlane.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CondFormats/DataRecord/interface/HeavyIonRPRcd.h"
#include "CondFormats/DataRecord/interface/HeavyIonRcd.h"
#include "CondFormats/HIObjects/interface/CentralityTable.h"
#include "CondCore/DBOutputService/interface/PoolDBOutputService.h"
#include "CondFormats/HIObjects/interface/RPFlatParams.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/Common/interface/RefProd.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefVector.h"

#include "RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneFlatten.h"
#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TTree.h"
#include "TH1I.h"
#include "TF1.h"
#include "TList.h"
#include "TString.h"
#include "TMath.h"
#include <time.h>
#include <cstdlib>
#include <vector>

#include "RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneList.h"
#include "RecoHI/HiEvtPlaneAlgos/interface/LoadEPDB.h"
#include "CondFormats/HIObjects/interface/CentralityTable.h"

using namespace std;
using namespace hi;

static const int nptbins = 28;
static const int MaxNumFlatBins = 200;

static const float ptbins[] = {
     0.0,   0.1,   0.2,   0.3,   0.4,   0.5,   0.6,   0.7,   0.8,   0.9, 
     1.0,   1.2,   1.4,   1.6,   2.0,   2.5,   3.0,   4.0,   6.0,  10.0,  
    15.0,  20.0,  30.0,  50.0, 100.0, 200.0, 400.0, 600.0,1000.0};

class EvtPlaneCalibTree : public edm::EDAnalyzer {
public:
  explicit EvtPlaneCalibTree(const edm::ParameterSet&);
  ~EvtPlaneCalibTree();
  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  edm::Service<TFileService> fs;
  
  std::string centralityVariable_;
  std::string centralityLabel_;
  std::string centralityMC_;
  
  edm::InputTag centralityBinTag_;
  edm::EDGetTokenT<int> centralityBinToken;
  edm::Handle<int> cbin_;

  edm::InputTag centralityTag_;
  edm::EDGetTokenT<reco::Centrality> centralityToken;
  edm::Handle<reco::Centrality> centrality_;
  
  edm::InputTag vertexTag_;
  edm::EDGetTokenT<std::vector<reco::Vertex>> vertexToken;
  edm::Handle<std::vector<reco::Vertex>> vertex_;
  
  edm::InputTag caloTag_;
  edm::EDGetTokenT<CaloTowerCollection> caloToken;
  edm::Handle<CaloTowerCollection> caloCollection_;
  
  edm::InputTag castorTag_;
  edm::EDGetTokenT<std::vector<reco::CastorTower>> castorToken;
  edm::Handle<std::vector<reco::CastorTower>> castorCollection_;
  
  edm::InputTag trackTag_;
  edm::EDGetTokenT<reco::TrackCollection> trackToken;
  edm::Handle<reco::TrackCollection> trackCollection_;
  
  edm::InputTag inputPlanesTag_;
  edm::EDGetTokenT<reco::EvtPlaneCollection> inputPlanesToken;
  edm::Handle<reco::EvtPlaneCollection> inputPlanes_;
  
  bool foundCentTag;
  bool useECAL_;
  bool useHCAL_;
  bool useTrack_;
  bool loadDB_;
  double minet_;
  double maxet_;
  double effm_;
  double minpt_;
  double maxpt_;
  double minvtx_;
  double maxvtx_;
  double delvtx_;
  double dzerr_;
  double chi2_;
  bool FirstEvent;
  int FlatOrder_;
  uint runno_; 
  float centval;
  int bin;
  int cbin;
  uint ntrkval;
  float vtx;
  double caloCentRef_;
  double caloCentRefWidth_;
  int caloCentRefMinBin_;
  int caloCentRefMaxBin_;
  int NumFlatBins_;
  int CentBinCompression_;
  
  int HFEtScale_;
  bool useNtrkBins_; 
  bool bypassCentrality_;
  int obins;
  TTree * tree;
  TH1D * fparams;
  TH1I * iparams;
  TH1D * hNtrkoff;
  TH1D * hcent;
  TH1D * hcentbins;
  TH2F * hpt;
  TH2F * hpt2;
  TH2F * hptcnt;
  
  //TH2D * hcastor;
  float ws[NumEPNames];
  float wc[NumEPNames];
  float ws_no_w[NumEPNames];
  float wc_no_w[NumEPNames];
  float s[NumEPNames];
  float c[NumEPNames];
  float pt[NumEPNames];
  float pt2[NumEPNames];
  int cnt[NumEPNames];
  double nCentBins_;
  bool genFlatPsi_;
  bool useOffsetPsi_;
  bool storeNames_;
  
private:
  
};
EvtPlaneCalibTree::EvtPlaneCalibTree(const edm::ParameterSet& iConfig) {
  foundCentTag = kTRUE;
  nCentBins_ = 200.;
  
  centralityVariable_ = iConfig.getParameter<std::string>("centralityVariable");
  if(iConfig.exists("nonDefaultGlauberModel")){
    centralityMC_ = iConfig.getParameter<std::string>("nonDefaultGlauberModel");
  }
  centralityLabel_ = centralityVariable_+centralityMC_;
  
  centralityBinTag_ = iConfig.getParameter<edm::InputTag>("centralityBinTag_");
  centralityBinToken = consumes<int>(centralityBinTag_);
  
  centralityTag_ = iConfig.getParameter<edm::InputTag>("centralityTag_");
  centralityToken = consumes<reco::Centrality>(centralityTag_);
  
  vertexTag_  = iConfig.getParameter<edm::InputTag>("vertexTag_");
  vertexToken = consumes<std::vector<reco::Vertex>>(vertexTag_);

  caloTag_ = iConfig.getParameter<edm::InputTag>("caloTag_");
  caloToken = consumes<CaloTowerCollection>(caloTag_);

  castorTag_ = iConfig.getParameter<edm::InputTag>("castorTag_");
  castorToken = consumes<std::vector<reco::CastorTower>>(castorTag_);

  trackTag_ = iConfig.getParameter<edm::InputTag>("trackTag_");
  trackToken = consumes<reco::TrackCollection>(trackTag_);

  inputPlanesTag_ = iConfig.getParameter<edm::InputTag>("inputPlanesTag_");
  inputPlanesToken = consumes<reco::EvtPlaneCollection>(inputPlanesTag_);

  genFlatPsi_ = iConfig.getUntrackedParameter<bool>("genFlatPsi_",true);
  FlatOrder_ = iConfig.getUntrackedParameter<int>("FlatOrder_", 9);
  NumFlatBins_ = iConfig.getUntrackedParameter<int>("NumFlatBins_",40);
  CentBinCompression_ = iConfig.getUntrackedParameter<int>("CentBinCompression_",5);
  caloCentRef_ = iConfig.getUntrackedParameter<double>("caloCentRef_",80.);
  caloCentRefWidth_ = iConfig.getUntrackedParameter<double>("caloCentRefWidth_",5.);
  if(NumFlatBins_ > MaxNumFlatBins) {
    cout<<"NumFlatBins set to max of "<<MaxNumFlatBins<<endl;
  }
  useOffsetPsi_ = iConfig.getUntrackedParameter<bool>("useOffsetPsi_",true);
  minet_ = iConfig.getUntrackedParameter<double>("minet_",-1.);
  maxet_ = iConfig.getUntrackedParameter<double>("maxet_",-1.);
  effm_ = iConfig.getUntrackedParameter<double>("effm_",0.0);
  minpt_ = iConfig.getUntrackedParameter<double>("minpt_",0.5);
  maxpt_ = iConfig.getUntrackedParameter<double>("maxpt_",3.0);
  minvtx_ = iConfig.getUntrackedParameter<double>("minvtx_",-25.);
  maxvtx_ = iConfig.getUntrackedParameter<double>("maxvtx_",25.);
  delvtx_ = iConfig.getUntrackedParameter<double>("delvtx_",5.);
  dzerr_ = iConfig.getUntrackedParameter<double>("dzerr_",10.);
  chi2_  = iConfig.getUntrackedParameter<double>("chi2_",40.);
  storeNames_ = 1;
  
  FirstEvent = kTRUE;
  hcent = fs->make<TH1D>("cent","cent",101,0,100);
  hcentbins = fs->make<TH1D>("centbins","centbins",201,0,200);
  
  //produces<reco::EvtPlaneCollection>("recoLevel");
  TString epnamesF = EPNames[0].data();
  epnamesF = epnamesF+"/F";
  for(int i = 0; i<NumEPNames; i++) if(i>0) epnamesF = epnamesF + ":" + EPNames[i].data() + "/F";  
  
  TString epnamesI = EPNames[0].data();
  epnamesI = epnamesI+"/I";
  for(int i = 0; i<NumEPNames; i++) if(i>0) epnamesI = epnamesI + ":" + EPNames[i].data() + "/I";  
  
  
  cout<<"=========================="<<endl;
  cout<<"EvtPlaneProducer:         "<<endl;
  cout<<"  NumFlatBins_:           "<<NumFlatBins_<<endl;
  cout<<"  CentBinCompression_:    "<<CentBinCompression_<<endl;
  cout<<"  minet_:                 "<<minet_<<endl;
  cout<<"  maxet_:                 "<<maxet_<<endl;
  cout<<"  minpt_:                 "<<minpt_<<endl;
  cout<<"  maxpt_:                 "<<maxpt_<<endl; 
  cout<<"  minvtx_:                "<<minvtx_<<endl;
  cout<<"  maxvtx_:                "<<maxvtx_<<endl;
  cout<<"  delvtx_:                "<<delvtx_<<endl;
  cout<<"  dzerr_                  "<<dzerr_<<endl;
  cout<<"  chi2_                   "<<chi2_<<endl;
  cout<<"=========================="<<endl;
  
  
  fparams = fs->make<TH1D>("fparams","fparams",10,0,10);
  iparams = fs->make<TH1I>("iparams","iparams",10,0,10);
  hNtrkoff = fs->make<TH1D>("Ntrkoff","Ntrkoff",1001,0,3000);
  int nvtxbins = (maxvtx_ - minvtx_)/delvtx_ + 0.1;
  
  
  hpt    = fs->make<TH2F>("hpt",   "hpt",   NumFlatBins_,0,NumFlatBins_,nvtxbins,minvtx_,maxvtx_);
  hpt2   = fs->make<TH2F>("hpt2",  "hpt2",  NumFlatBins_,0,NumFlatBins_,nvtxbins,minvtx_,maxvtx_);
  hptcnt = fs->make<TH2F>("hptcnt","hptcnt",NumFlatBins_,0,NumFlatBins_,nvtxbins,minvtx_,maxvtx_);
  
  
  fparams->SetBinContent(1,minet_);
  fparams->SetBinContent(2,maxet_);
  fparams->SetBinContent(3,minpt_);
  fparams->SetBinContent(4,maxpt_);
  fparams->SetBinContent(5,minvtx_);
  fparams->SetBinContent(6,maxvtx_);
  fparams->SetBinContent(7,caloCentRef_);
  fparams->SetBinContent(8,caloCentRefWidth_);
  iparams->SetBinContent(1,FlatOrder_);
  iparams->SetBinContent(2,NumFlatBins_);
  iparams->SetBinContent(3,CentBinCompression_);
  tree = fs->make<TTree>("tree","EP tree");
  tree->Branch("Cent",    &centval,    "cent/F");
  tree->Branch("Vtx",     &vtx,        "vtx/F");
  tree->Branch("Run",     &runno_,     "run/i");
  tree->Branch("bin",     &bin,        "bin/I");
  tree->Branch("cbin",     &cbin,        "cbin/I");
  tree->Branch("NtrkOff",&ntrkval,     "ntrkoff/i");
  tree->Branch("ws", &ws, epnamesF.Data());
  tree->Branch("wc", &wc, epnamesF.Data());
  tree->Branch("ws_no_w", &ws_no_w, epnamesF.Data());
  tree->Branch("wc_no_w", &wc_no_w, epnamesF.Data());
  tree->Branch("pt", &pt, epnamesF.Data());
  tree->Branch("pt2", &pt2, epnamesF.Data());
  tree->Branch("cnt",&cnt,epnamesI.Data());
  runno_ = 0;
}


EvtPlaneCalibTree::~EvtPlaneCalibTree()
{
  
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  
}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
EvtPlaneCalibTree::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  centval = -1;
  vtx = -100;
  ntrkval = -1;

  using namespace edm;
  using namespace std;
  using namespace reco;
  
  Bool_t newrun = kFALSE;
  if((int) runno_ != (int) iEvent.id().run()) newrun = kTRUE;
  runno_ = iEvent.id().run();
  //////////////////
  if(FirstEvent || newrun) {
    FirstEvent = kFALSE;
    newrun = kFALSE;
    //
    //Get Size of Centrality Table
    //
    edm::ESHandle<CentralityTable> centDB_;
    iSetup.get<HeavyIonRcd>().get(centralityLabel_,centDB_);
    nCentBins_ = (int) centDB_->m_table.size();
  } //First event
  
  //
  //Get Centrality
  //
  
  iEvent.getByToken(centralityBinToken, cbin_);
  cbin = *cbin_;
  bin = cbin/CentBinCompression_; 
  
  double cscale = 100./nCentBins_;
  centval = cscale*cbin;
  
  hcent->Fill(centval);
  hcentbins->Fill(cbin);
  bin = cbin/CentBinCompression_; 
  
  //
  //Get Vertex
  //
  int vs_sell;
  float vzr_sell;
  iEvent.getByLabel(vertexTag_,vertex_);
  const reco::VertexCollection * vertices3 = vertex_.product();
  vs_sell = vertices3->size();
  if(vs_sell>0) {
    vzr_sell = vertices3->begin()->z();
  } else
    vzr_sell = -999.9;
  //
  vtx = vzr_sell;
  if(vzr_sell>minvtx_ && vzr_sell<maxvtx_) {    
    //
    //Get Event Planes
    //
    iEvent.getByToken(inputPlanesToken,inputPlanes_);
    
    if(!inputPlanes_.isValid()){
      cout << "Error! Can't get hiEvtPlane product!" << endl;
      return ;
    }
    
    int i = 0;
    for (EvtPlaneCollection::const_iterator ep = inputPlanes_->begin();ep !=inputPlanes_->end(); ep++) {
      if(i>=NumEPNames || i!= ep->indx()) {
	cout<<"DATA INCONSISTENCY.  Too many or wrong EPs found in collection. Expect a crash!"<<endl;
      }

      ws[i]=0;
      wc[i]=0;
      ws_no_w[i] = 0;
      wc_no_w[i] = 0;
      cnt[i]=0;
      pt[i]=0;
      pt2[i]=0;


      double c = ep->sumCos(0);
      double s = ep->sumSin(0);
      double c_no_w = ep->sumCos(3);
      double s_no_w = ep->sumSin(3);
      uint m = ep->mult();
      float fpt = ep->sumPtOrEt();
      float fpt2 = ep->sumPtOrEt2();

      if(ep->angle(0)>-5) {      
	ws[i]=s;
	wc[i]=c;
	ws_no_w[i] = s_no_w;
	wc_no_w[i] = c_no_w;
	cnt[i]=m;
	pt[i]=fpt;
	pt2[i]=fpt2;
      }
      ++i;
    }    
    tree->Fill();   
  }
}

// ------------ method called once each job just before starting event loop  ------------
void 
EvtPlaneCalibTree::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
EvtPlaneCalibTree::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(EvtPlaneCalibTree);
