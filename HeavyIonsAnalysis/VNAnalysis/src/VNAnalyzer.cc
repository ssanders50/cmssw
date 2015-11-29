// -*- C++ -*-
//
// Package:    VNAnalyzer
// Class:      VNAnalyzer
// 


// system include files
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

#include "DataFormats/HeavyIonEvent/interface/Centrality.h"
#include "DataFormats/HeavyIonEvent/interface/EvtPlane.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CondFormats/DataRecord/interface/HeavyIonRPRcd.h"
#include "CondFormats/DataRecord/interface/HeavyIonRcd.h"
#include "CondFormats/HIObjects/interface/CentralityTable.h"
#include "CondFormats/HIObjects/interface/RPFlatParams.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"
#include "HeavyIonsAnalysis/VNAnalysis/interface/TrackEfficiency.h"

#include "TROOT.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TTree.h"
#include "TH1I.h"
#include "TF1.h"
#include "TMath.h"
#include "TRandom.h"
#include <time.h>
#include <cstdlib>
	
#include <vector>
using std::vector;
using std::rand;
using namespace std;
#include "RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneFlatten.h"
#include "RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneList.h"
#include "RecoHI/HiEvtPlaneAlgos/interface/LoadEPDB.h"
using namespace hi;

static const int nptbinsDefault = 16;
static const double ptbinsDefault[]={
  0.2,  0.3,  0.4,  0.5,  0.6,  0.8,  1.0,  1.2,  1.6,  2.0,
  2.5,  3.0,  3.5,  4.0,  5.0,  6.0,  8.0};
static const int netabinsDefault = 12;
static const double etabinsDefault[]={-2.4, -2.0, -1.6, -1.2,
				      -0.8, -0.4, 0.0,  0.4,  0.8,
				      1.2,  1.6,  2.0,  2.4};


//
// class declaration
//

class VNAnalyzer : public edm::EDAnalyzer {
public:
  explicit VNAnalyzer(const edm::ParameterSet&);
  ~VNAnalyzer();
      
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  // ----------member data ---------------------------
  int eporder_;


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

  edm::InputTag trackTag_;
  edm::EDGetTokenT<reco::TrackCollection> trackToken;
  edm::Handle<reco::TrackCollection> trackCollection_;

  edm::InputTag inputPlanesTag_;
  edm::EDGetTokenT<reco::EvtPlaneCollection> inputPlanesToken;
  edm::Handle<reco::EvtPlaneCollection> inputPlanes_;

  edm::Service<TFileService> fs;

  double caloCentRef_;
  double caloCentRefWidth_;
  int caloCentRefMinBin_;
  int caloCentRefMaxBin_;

  double nCentBins_;

  int vs_sell;   // vertex collection size
  float vzr_sell;
  float vzErr_sell;
  TH1D * hcent;
  TH1D * hcentbins;
  TH2D * hEff[10];
  TH2D * hw[10];
  double centval;
  double ntrkval;
  double vtx;
  int Noff;

  Double_t epang[NumEPNames];
  Double_t epsin[NumEPNames];
  Double_t epcos[NumEPNames];

  Double_t qx[NumEPNames];
  Double_t qy[NumEPNames];
  Double_t q[NumEPNames];
  Double_t epmult[NumEPNames];
  Double_t vn[NumEPNames];

  Double_t rescor[NumEPNames];
  Double_t rescorErr[NumEPNames];


  unsigned int runno_;

  TH1D * hNtrkoff;
  int nEtaBins;
  TH1I * hrun;
  string rpnames[NumEPNames];
  string effTable_;
  TTree * tree;
  TrackEfficiency *teff;
  int FlatOrder_;
  int NumFlatBins_;
  int CentBinCompression_;
  int Noffmin_;
  int Noffmax_;
  int EPOrder_;
  TH2D * qxtrk;
  TH2D * qytrk;
  TH2D * qcnt;
  TH2D * wqxtrk;
  TH2D * wqytrk;
  TH2D * wqcnt;
  TH2D * weff;
  TH2D * avpt;
  HiEvtPlaneFlatten * flat[NumEPNames];
  bool loadDB_;
  bool FirstEvent_;

  int getNoff(const edm::Event& iEvent, const edm::EventSetup& iSetup, double cent)
  {
    int Noff = 0;
    using namespace edm;
    using namespace reco;
    qxtrk->Reset();
    qytrk->Reset();
    qcnt->Reset();
    if(teff) {
      wqxtrk->Reset();
      wqytrk->Reset();
      wqcnt->Reset();
      weff->Reset();
    }
    avpt->Reset();
    iEvent.getByToken(vertexToken,vertex_);
    const VertexCollection * recoVertices = vertex_.product();
    int primaryvtx = 0;
    math::XYZPoint v1( (*recoVertices)[primaryvtx].position().x(), (*recoVertices)[primaryvtx].position().y(), (*recoVertices)[primaryvtx].position().z() );
    double vxError = (*recoVertices)[primaryvtx].xError();
    double vyError = (*recoVertices)[primaryvtx].yError();
    double vzError = (*recoVertices)[primaryvtx].zError();
    iEvent.getByLabel(trackTag_,trackCollection_);
    for(TrackCollection::const_iterator itTrack = trackCollection_->begin();
	itTrack != trackCollection_->end();                      
	++itTrack) {    
      if ( !itTrack->quality(reco::TrackBase::highPurity) ) continue;
      if ( itTrack->charge() == 0 ) continue;
      double d0 = -1.* itTrack->dxy(v1);
      double derror=sqrt(itTrack->dxyError()*itTrack->dxyError()+vxError*vyError);
      double dz=itTrack->dz(v1);
      double dzerror=sqrt(itTrack->dzError()*itTrack->dzError()+vzError*vzError);
      if ( fabs(itTrack->eta()) > 2.4 ) continue;
      if ( fabs( dz/dzerror ) > 3. ) continue;
      if ( fabs( d0/derror ) > 3. ) continue;
      if ( itTrack->ptError()/itTrack->pt() > 0.1 ) continue;

      qxtrk->Fill(itTrack->pt(), itTrack->eta(), TMath::Cos(EPOrder_*itTrack->phi()));
      qytrk->Fill(itTrack->pt(), itTrack->eta(), TMath::Sin(EPOrder_*itTrack->phi()));
      qcnt->Fill(itTrack->pt(), itTrack->eta());
      avpt->Fill(itTrack->pt(), itTrack->eta(), itTrack->pt());

      if(teff) {
	double w = teff->getEfficiencies(itTrack->pt(),cent,itTrack->phi(),itTrack->eta());
	if(w>0.0 ) {
	  wqxtrk->Fill(itTrack->pt(), itTrack->eta(), TMath::Cos(EPOrder_*itTrack->phi())/w);
	  wqytrk->Fill(itTrack->pt(), itTrack->eta(), TMath::Sin(EPOrder_*itTrack->phi())/w);
	  wqcnt->Fill(itTrack->pt(), itTrack->eta());
	  weff->Fill(itTrack->pt(), itTrack->eta(),1/w);
          hw[(int)(cent/10.)]->Fill(itTrack->phi(),itTrack->eta(), 1/w);
          hEff[(int)(cent/10.)]->Fill(itTrack->phi(),itTrack->eta(), 1.);

	}
      }
      if ( itTrack->pt() < 0.4 ) continue;
      Noff++;
      if( itTrack->pt() < 0.5 ) continue;
      if(!teff) hEff[(int)(cent/10.)]->Fill(itTrack->phi(),itTrack->eta());
    }

    if(Noff < Noffmin_ || Noff > Noffmax_) return -2;
    return Noff;
  }


};


//
// constructors and destructor
//
VNAnalyzer::VNAnalyzer(const edm::ParameterSet& iConfig):runno_(0)
  
{
  runno_ = 0;
  loadDB_ = kTRUE;
  FirstEvent_ = kTRUE;

  for(int i = 0; i<NumEPNames; i++) {
    epang[i] = -10;
    epsin[i] = 0;
    epcos[i] = 0;
    qx[i] = 0;
    qy[i] = 0;
    q[i] = 0;
    epmult[i] = 0;
    rescor[i] = 0;
    rescorErr[i] = 0;
  }

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

  trackTag_ = iConfig.getParameter<edm::InputTag>("trackTag_");
  trackToken = consumes<reco::TrackCollection>(trackTag_);

  inputPlanesTag_ = iConfig.getParameter<edm::InputTag>("inputPlanesTag_");
  inputPlanesToken = consumes<reco::EvtPlaneCollection>(inputPlanesTag_);

  EPOrder_ = iConfig.getUntrackedParameter<int>("EPOrder_",2);
  FlatOrder_ = iConfig.getUntrackedParameter<int>("FlatOrder_", 9);
  NumFlatBins_ = iConfig.getUntrackedParameter<int>("NumFlatBins_",20);
  caloCentRef_ = iConfig.getUntrackedParameter<double>("caloCentRef_",80.);
  caloCentRefWidth_ = iConfig.getUntrackedParameter<double>("caloCentRefWidth_",5.);
  CentBinCompression_ = iConfig.getUntrackedParameter<int>("CentBinCompression_",5);
  Noffmin_ = iConfig.getUntrackedParameter<int>("Noffmin_", 0);
  Noffmax_ = iConfig.getUntrackedParameter<int>("Noffmax_", 50000);	
  effTable_ = iConfig.getParameter<std::string>("effTable_");
  teff = 0;
  if(!effTable_.empty()) teff = new TrackEfficiency(effTable_.data());

  hNtrkoff = fs->make<TH1D>("Ntrkoff","Ntrkoff",1001,0,3000);
  qxtrk = fs->make<TH2D>(Form("qxtrk_v%d",EPOrder_),Form("qxtrk_v%d",EPOrder_),nptbinsDefault,ptbinsDefault, netabinsDefault, etabinsDefault);
  qytrk = fs->make<TH2D>(Form("qytrk_v%d",EPOrder_),Form("qytrk_v%d",EPOrder_),nptbinsDefault,ptbinsDefault, netabinsDefault, etabinsDefault);
  qcnt =  fs->make<TH2D>(Form("qcnt_v%d",EPOrder_), Form("qcnt_v%d",EPOrder_),nptbinsDefault,ptbinsDefault, netabinsDefault, etabinsDefault);
  if(teff) {
    wqxtrk = fs->make<TH2D>(Form("wqxtrk_v%d",EPOrder_),Form("wqxtrk_v%d",EPOrder_),nptbinsDefault,ptbinsDefault, netabinsDefault, etabinsDefault);
    wqytrk = fs->make<TH2D>(Form("wqytrk_v%d",EPOrder_),Form("wqytrk_v%d",EPOrder_),nptbinsDefault,ptbinsDefault, netabinsDefault, etabinsDefault);
    wqcnt =  fs->make<TH2D>(Form("wqcnt_v%d",EPOrder_), Form("wqcnt_v%d",EPOrder_),nptbinsDefault,ptbinsDefault, netabinsDefault, etabinsDefault);
    weff =  fs->make<TH2D>("weff","weff",nptbinsDefault,ptbinsDefault, netabinsDefault, etabinsDefault);
  }
  avpt =  fs->make<TH2D>("avpt","avpt",nptbinsDefault,ptbinsDefault, netabinsDefault, etabinsDefault);

  hcent = fs->make<TH1D>("cent","cent",220,-10,110);
  hcentbins = fs->make<TH1D>("centbins","centbins",201,0,200);
  hrun = fs->make<TH1I>("runs","runs",100000,150001,250000);

  for(int i = 0; i<10; i++) {
    TString hn = Form("Eff_%d_%d",10*i,10*(i+1));
    hEff[i] = fs->make<TH2D>(hn.Data(),hn.Data(),50,-TMath::Pi(),TMath::Pi(),50,-2.4,2.4);
    hEff[i]->Sumw2();
    hEff[i]->SetXTitle("#phi (radians)");
    hEff[i]->SetYTitle("#eta");
    hEff[i]->SetOption("colz");

    TString hnw = Form("Effw_%d_%d",10*i,10*(i+1));
    hw[i] = fs->make<TH2D>(hnw.Data(),hnw.Data(),50,-TMath::Pi(),TMath::Pi(),50,-2.4,2.4);
    hw[i]->Sumw2();
    hw[i]->SetXTitle("#phi (radians)");
    hw[i]->SetYTitle("#eta");
    hw[i]->SetOption("colz");
  }
	      
  TString epnames = EPNames[0].data();
  epnames = epnames+"/D";
  for(int i = 0; i<NumEPNames; i++) {
    if(i>0) epnames = epnames + ":" + EPNames[i].data() + "/D";
    TFileDirectory subdir = fs->mkdir(Form("%s",EPNames[i].data()));
    flat[i] = new HiEvtPlaneFlatten();
    flat[i]->init(FlatOrder_,NumFlatBins_,EPNames[i],EPOrder[i]);

  }
  tree = fs->make<TTree>("tree","EP tree");

  tree->Branch("Cent",&centval,"cent/D");
  tree->Branch("NtrkOff",&Noff,"Noff/I");
  tree->Branch("Vtx",&vtx,"vtx/D");
  tree->Branch("epang",&epang, epnames.Data());
  tree->Branch("qx",      &qx,       epnames.Data());
  tree->Branch("qy",      &qy,       epnames.Data());
  tree->Branch("q",       &q,       epnames.Data());
  tree->Branch("vn", &vn, epnames.Data());
  tree->Branch("mult",    &epmult,  epnames.Data());
  tree->Branch("Run",     &runno_,   "run/i");
  tree->Branch("Rescor",  &rescor,   epnames.Data());
  tree->Branch("RescorErr",  &rescorErr,   epnames.Data());
  tree->Branch(Form("qxtrk_v%d",EPOrder_),   "TH2D",  &qxtrk, 128000, 0);
  tree->Branch(Form("qytrk_v%d",EPOrder_),   "TH2D",  &qytrk, 128000, 0);
  tree->Branch(Form("qcnt_v%d",EPOrder_),    "TH2D",  &qcnt, 128000, 0);
  if(teff) {
    tree->Branch(Form("wqxtrk_v%d",EPOrder_),   "TH2D",  &wqxtrk, 128000, 0);
    tree->Branch(Form("wqytrk_v%d",EPOrder_),   "TH2D",  &wqytrk, 128000, 0);
    tree->Branch(Form("wqcnt_v%d",EPOrder_),    "TH2D",  &wqcnt, 128000, 0);
    tree->Branch(Form("weff_v%d",EPOrder_),    "TH2D",  &weff, 128000, 0);
  }
  tree->Branch("avpt",    "TH2D",  &avpt, 128000, 0);

}



VNAnalyzer::~VNAnalyzer()
{
  
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
VNAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  using namespace reco;
  Bool_t newrun = kFALSE;
  if(runno_ != iEvent.id().run()) newrun = kTRUE;
  runno_ = iEvent.id().run();
  hrun->Fill(runno_);

  if(FirstEvent_ || newrun) {
    FirstEvent_ = kFALSE;
    newrun = kFALSE;
    //
    //Get Size of Centrality Table
    //
    edm::ESHandle<CentralityTable> centDB_;
    iSetup.get<HeavyIonRcd>().get(centralityLabel_,centDB_);
    nCentBins_ = (int) centDB_->m_table.size();
    for(int i = 0; i<NumEPNames; i++) {
      flat[i]->setCaloCentRefBins(-1,-1);
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
    cout<<"Load flattening parameters"<<endl;
    //
    //Get flattening parameter file.  
    //
    edm::ESHandle<RPFlatParams> flatparmsDB_;
    iSetup.get<HeavyIonRPRcd>().get(flatparmsDB_);
    LoadEPDB * db = new LoadEPDB(flatparmsDB_,flat);
    if(!db->IsSuccess()) {
      cout<<"Failed to load DB!"<<endl;
      loadDB_ = kFALSE;
    }        
    if(teff) teff->setRunNumber(runno_);
  } //First event


  //
  //Get Centrality
  //
  //int bin = 0;
  int Noff=0;

  iEvent.getByToken(centralityBinToken, cbin_);
  int cbin = *cbin_;
  //bin = cbin/CentBinCompression_; 
  double cscale = 100./nCentBins_;
  centval = cscale*cbin;

  Noff = getNoff( iEvent, iSetup,centval);
  hNtrkoff->Fill(Noff);

  hcent->Fill(centval);
  hcentbins->Fill(cbin);
  //
  //Get Vertex
  //
  int vs_sell;   // vertex collection size
  float vzr_sell;
  iEvent.getByToken(vertexToken,vertex_);
  const reco::VertexCollection * vertices3 = vertex_.product();
  vs_sell = vertices3->size();
  if(vs_sell>0) {
    vzr_sell = vertices3->begin()->z();
  } else
    vzr_sell = -999.9;
  
  vtx = vzr_sell;
  //
  //Get Event Planes
  //
  iEvent.getByToken(inputPlanesToken,inputPlanes_);
  
  if(!inputPlanes_.isValid()){
    cout << "Error! Can't get hiEvtPlaneFlat product!" << endl;
    return ;
  }
  
  Int_t indx = 0;
  for (EvtPlaneCollection::const_iterator rp = inputPlanes_->begin();rp !=inputPlanes_->end(); rp++) {
    if(indx != rp->indx() ) {
      cout<<"EP inconsistency found. Abort."<<endl;
      return;
    }
    if(rp->sumSin()!=0 || rp->sumCos()!=0) {
      epang[indx]=rp->angle();
      epsin[indx] = rp->sumSin();
      epcos[indx] = rp->sumCos();
      
      qx[indx]  = rp->qx(); 
      qy[indx]  = rp->qy();
      q[indx]   = rp->q();
      vn[indx] = rp->vn(0);
      epmult[indx] = (double) rp->mult();
      
      rescor[indx] = flat[indx]->getCentRes1((int) centval);
      rescorErr[indx] = flat[indx]->getCentResErr1((int) centval);
    }
    ++indx; 
  }


  ntrkval = Noff;
  if ( Noff == -2 ) {
    return;
  }
  hNtrkoff->Fill(Noff);

  tree->Fill(); 
}



// ------------ method called once each job just before starting event loop  ------------
void 
VNAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
VNAnalyzer::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(VNAnalyzer);

