#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH1I.h"
#include "TDirectory.h"
#include <iostream>
#include "string.h"
#include "stdio.h"

#include "EPCalib/HiEvtPlaneList.h"
#include "EPCalib/HiEvtPlaneFlatten.h"


using namespace hi;
using namespace std;

static const int nCentBins_ = 200;
static const int nptbins = 28;

static const float ptbins[] = {
     0.0,   0.1,   0.2,   0.3,   0.4,   0.5,   0.6,   0.7,   0.8,   0.9, 
     1.0,   1.2,   1.4,   1.6,   2.0,   2.5,   3.0,   4.0,   6.0,  10.0,  
    15.0,  20.0,  30.0,  50.0, 100.0, 200.0, 400.0, 600.0,1000.0};

struct sout_struct{
  double ws[NumEPNames];
  double wc[NumEPNames];
  double wsum[NumEPNames];
  double msum[NumEPNames];
  double vtx;
  double cent;
  int bin;
  unsigned int ntrk;
  unsigned int run;
}; 
struct sout_struct sout; 
Bool_t test = kFALSE;

float centval;
float vtx;
unsigned int runno_;
unsigned int ntrkval;
int bin;
// TH1F * s[NumEPNames];
// TH1F * c[NumEPNames];
// TH1I * cnt[NumEPNames];
// TH1F * pt[NumEPNames];
// TH1F * pt2[NumEPNames];
TTree * tree;
TH1D * hcent;
TH1D * hcentbin;
TH1D * hvtx;
TH1D * hPsi[NumEPNames];

TH1D * hPsiOffset[NumEPNames];
TH1D * hPsiOffset2[NumEPNames];
TH1D * hPsiFlat[NumEPNames];
TDirectory * dirs[NumEPNames];
HiEvtPlaneFlatten * flat[NumEPNames];
HiEvtPlaneFlatten * flatOffset[NumEPNames];

int FlatOrder_ ;
int NumFlatBins_ ;
int CentBinCompression_;
double minet_;
double maxet_;
double minpt_;
double maxpt_;
double minvtx_;
double maxvtx_;
double caloCentRef_;
double caloCentRefWidth_;
TH1D * flatXhist[NumEPNames];
TH1D * flatYhist[NumEPNames];
TH1D * flatCnthist[NumEPNames];
TH1D * pthist[NumEPNames];
TH1D * pt2hist[NumEPNames];
TH1D * ptcnthist[NumEPNames];

TH1D * xoffhist[NumEPNames];
TH1D * yoffhist[NumEPNames];
TH1D * xyoffcnthist[NumEPNames];
TH1D * xyoffmulthist[NumEPNames];

TH1D * hpt[NumEPNames][10];
TH1D * hscale[NumEPNames][10];


float EPcent;
float EPvtx;
unsigned int EPntrk;
unsigned int EPrun;
float wsv[NumEPNames];
float wcv[NumEPNames];
float wsv_no_w[NumEPNames];
float wcv_no_w[NumEPNames];
float pt[NumEPNames];
float pt2[NumEPNames];
float ptav[NumEPNames];
float pt2av[NumEPNames];
float ptcntav[NumEPNames];

double pt_[NumEPNames][1000];
double pt2_[NumEPNames][1000];
double ptcnt_[NumEPNames][1000];

int ptcnt[NumEPNames];
TTree * EPtree;
float EPAngs[NumEPNames];
int nentries;
int totentries;
int ncnt;
FILE * save;
#include "EPCalib/src/Loop0.h"
#include "EPCalib/src/Loop1.h"
#include "EPCalib/src/Loop2.h"
#include "EPCalib/src/Loop3.h"
#include "EPCalib/src/rescor.h"

void EPCalib(){
  save = fopen("save","wb");
  int Hbins = 0;
  int Obins = 0;
  char buf[200];
  FILE * list;
  list = fopen("tmp.lis","r");
  TString fnames[1000];
  memset(buf,0,sizeof(buf));
  int lcnt = 0;
  while(fgets(buf,200,list)!=NULL) {
    buf[strlen(buf)-1]=0;
    fnames[lcnt] = buf;
    memset(buf,0,sizeof(buf));
    ++lcnt;
  }
  fclose(list);

  TFile * tfout = new TFile("EP.root","recreate");
  TDirectory * outdir = tfout->mkdir("hiEvtPlaneFlatCalib");
  outdir->cd();
  TString epnames = EPNames[0].data();
  epnames = epnames+"/F";
  for(int i = 0; i<NumEPNames; i++) if(i>0) epnames = epnames + ":" + EPNames[i].data() + "/F";  

  EPtree = new TTree("EPtree","EP tree");
  EPtree->Branch("Cent", &EPcent,"cent/F");
  EPtree->Branch("Vtx",&EPvtx,"vtx/F");
  EPtree->Branch("NtrkOff", &EPntrk, "noff/i");
  EPtree->Branch("Run",&EPrun,"run/i");
  EPtree->Branch("EPAngs",&EPAngs,epnames.Data());
  hcent = new TH1D("hcent","centrality",200,0,100);
  hcentbin = new TH1D("hcentbin","centrality bin",400,0,400);
  hvtx = new TH1D("hvtx","vertex",200,-30,30);
  TFile * tf = new TFile(fnames[0].Data(),"read");
  TH1D * fparams = (TH1D *)tf->Get("evtPlaneCalibTree/fparams");
  TH1I * iparams = (TH1I *)tf->Get("evtPlaneCalibTree/iparams");
  minet_ = fparams->GetBinContent(1);
  maxet_ = fparams->GetBinContent(2);
  minpt_ = fparams->GetBinContent(3);
  maxpt_ = fparams->GetBinContent(4);
  minvtx_ = fparams->GetBinContent(5);
  maxvtx_ = fparams->GetBinContent(6);
  caloCentRef_ = fparams->GetBinContent(7);
  caloCentRefWidth_ = fparams->GetBinContent(8);
  FlatOrder_ = iparams->GetBinContent(1);
  NumFlatBins_ = iparams->GetBinContent(2);
  CentBinCompression_ = iparams->GetBinContent(3);

  cout<<"==============================================="<<endl;
  cout<<"minet_               "<<minet_<<endl;
  cout<<"maxet_               "<<maxet_<<endl;
  cout<<"minpt_               "<<minpt_<<endl;
  cout<<"maxpt_               "<<maxpt_<<endl;
  cout<<"minvtx_              "<<minvtx_<<endl;
  cout<<"maxvtx_              "<<maxvtx_<<endl;
  cout<<"caloCentRef_         "<<caloCentRef_<<endl;
  cout<<"caloCentRefWidth_    "<<caloCentRefWidth_<<endl;
  cout<<"FlatOrder_           "<<FlatOrder_<<endl;
  cout<<"NumFlatBins_         "<<NumFlatBins_<<endl;
  cout<<"CentBinCompression_  "<<CentBinCompression_<<endl;
  cout<<"==============================================="<<endl;
  tf->Close();
  for(int i = 0; i<NumEPNames; i++) {
    dirs[i] = outdir->mkdir(Form("%s",EPNames[i].data()));
    dirs[i]->cd();
    Double_t psirange = 4;
    if(EPOrder[i]==1 ) psirange = 3.5;
    if(EPOrder[i]==2 ) psirange = 2;
    if(EPOrder[i]==3 ) psirange = 1.5;
    if(EPOrder[i]==4 ) psirange = 1;
    if(EPOrder[i]==5) psirange = 0.8;
    if(EPOrder[i]==6) psirange = 0.6;
    for(int j = 0; j<10; j++) {
      hpt[i][j] = new TH1D(Form("pt_%d",j),Form("pt_%d",j),500,0,5);
      hpt[i][j]->SetXTitle("pt (GeV");
      hpt[i][j]->SetYTitle("Counts");
      
      hscale[i][j] = new TH1D(Form("scale_%d",j),Form("scale_%d",j),500,0,2);
      hscale[i][j]->SetXTitle("pt2 (GeV");
      hscale[i][j]->SetYTitle("Counts");
    }

    hPsi[i] = new TH1D("psi","psi",800,-psirange,psirange);
    hPsi[i]->SetXTitle("#Psi");
    hPsi[i]->SetYTitle(Form("Counts (cent<80%c)",'%'));
    
    hPsiOffset[i] = new TH1D("psiOffset","psiOffset",800,-psirange,psirange);
    hPsiOffset[i]->SetXTitle("#Psi");
    hPsiOffset[i]->SetYTitle(Form("Counts (cent<80%c)",'%'));

    hPsiOffset2[i] = new TH1D("psiOffset2","psiOffset2",800,-psirange,psirange);
    hPsiOffset2[i]->SetXTitle("#Psi");
    hPsiOffset2[i]->SetYTitle(Form("Counts (cent<80%c)",'%'));

    flat[i] = new HiEvtPlaneFlatten();
    flat[i]->init(FlatOrder_,NumFlatBins_,EPNames[i].data(),EPOrder[i]);
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
    
    flatOffset[i] = new HiEvtPlaneFlatten();
    flatOffset[i]->init(FlatOrder_,NumFlatBins_,EPNames[i].data(),EPOrder[i]);
    if(caloCentRef_>0) {
      int minbin = (caloCentRef_-caloCentRefWidth_/2.)*nCentBins_/100.;
      int maxbin = (caloCentRef_+caloCentRefWidth_/2.)*nCentBins_/100.;
      minbin/=CentBinCompression_;
      maxbin/=CentBinCompression_;
      if(minbin>0 && maxbin>=minbin) {
	if(EPDet[i]==HF || EPDet[i]==Castor) flatOffset[i]->setCaloCentRefBins(minbin,maxbin);
      }
    }
    
    Hbins = flat[i]->getHBins();
    Obins = flat[i]->getOBins();
    int nbins = flat[i]->getHBins() + 2*flat[i]->getOBins();
    flatXhist[i] = new TH1D(Form("x_%s",EPNames[i].data()),Form("x_%s",EPNames[i].data()),nbins,-0.5,nbins-0.5);
    flatYhist[i] = new TH1D(Form("y_%s",EPNames[i].data()),Form("y_%s",EPNames[i].data()),nbins,-0.5,nbins-0.5);   
    flatCnthist[i] = new TH1D(Form("cnt_%s",EPNames[i].data()),Form("cnt_%s",EPNames[i].data()),nbins,-0.5,nbins-0.5);
    xoffhist[i] = new TH1D(Form("xoff_%s",EPNames[i].data()),Form("xoff_%s",EPNames[i].data()),Obins,-0.5,Obins-0.5);
    yoffhist[i] = new TH1D(Form("yoff_%s",EPNames[i].data()),Form("yoff_%s",EPNames[i].data()),Obins,-0.5,Obins-0.5);
    xyoffcnthist[i] = new TH1D(Form("xyoffcnt_%s",EPNames[i].data()),Form("xyoffcnt_%s",EPNames[i].data()),Obins,-0.5,Obins-0.5);
    xyoffmulthist[i] = new TH1D(Form("xyoffmult_%s",EPNames[i].data()),Form("xyoffmult_%s",EPNames[i].data()),Obins,-0.5,Obins-0.5);
    pthist[i] = new TH1D(Form("pt_%s",EPNames[i].data()),Form("pt_%s",EPNames[i].data()),Obins,-0.5,Obins-0.5);
    pt2hist[i] = new TH1D(Form("pt2_%s",EPNames[i].data()),Form("pt2_%s",EPNames[i].data()),Obins,-0.5,Obins-0.5);
    ptcnthist[i] = new TH1D(Form("ptcnt_%s",EPNames[i].data()),Form("ptcnt_%s",EPNames[i].data()),Obins,-0.5,Obins-0.5);
    
    hPsiFlat[i] = new TH1D("psiFlat","psiFlat",800,-psirange,psirange);
    hPsiFlat[i]->SetXTitle("#Psi");
    hPsiFlat[i]->SetYTitle(Form("Counts (cent<80%c)",'%'));
    ptav[i] = 0;
    pt2av[i] = 0;
    ptcntav[i] = 0;
  }
  //
  // Loop 0
  //
  totentries = 0;
  ncnt = 0;
  for(int i = 0; i<NumEPNames; ++i){
    for(int j = 0; j<1000; j++) {
      pt_[i][j]=0;
      pt2_[i][j]=0;
      ptcnt_[i][j]=0;
    }
  }
  for(int findx = 0; findx< lcnt; findx++) {
    tf = new TFile(fnames[findx].Data(),"read");
    tree = (TTree *) tf->Get("evtPlaneCalibTree/tree");
    tree->SetBranchAddress("Cent",    &centval);
    tree->SetBranchAddress("Vtx",     &vtx);
    tree->SetBranchAddress("Run",     &runno_);
    tree->SetBranchAddress("NtrkOff",&ntrkval);
    tree->SetBranchAddress("bin",    &bin);
    tree->SetBranchAddress("ws", &wsv);
    tree->SetBranchAddress("wc", &wcv);
    tree->SetBranchAddress("ws_no_w", &wsv_no_w);
    tree->SetBranchAddress("wc_no_w", &wcv_no_w);
    tree->SetBranchAddress("pt", &pt);
    tree->SetBranchAddress("pt2",&pt2);
    tree->SetBranchAddress("cnt",&ptcnt);
    
    nentries = tree->GetEntries();
    cout<<"Loop 0: "<<tf->GetName()<<"\t"<<nentries<<"\t"<<totentries<<endl;
    if(test) nentries = 1000;
    Loop0();
    for(int i = 0; i<NumEPNames; i++) {
      for(int j = 0; j<1000; j++) {
	if(ptcnt_[i][j]>0) {
	  double cval = ptcnt_[i][j];
	  flat[i]->setPtDB(j,pt_[i][j]/cval);
	  flat[i]->setPt2DB(j,pt2_[i][j]/cval);
	}
      }
      if(ptcntav[i]>0) {
	ptav[i]/=ptcntav[i];
	pt2av[i]/=ptcntav[i];
      }
    }
  }
  //
  // Loop 1
  //
  totentries = 0;
  ncnt = 0;
  for(int findx = 0; findx< lcnt; findx++) {
    tf = new TFile(fnames[findx].Data(),"read");
    tree = (TTree *) tf->Get("evtPlaneCalibTree/tree");
    tree->SetBranchAddress("Cent",    &centval);
    tree->SetBranchAddress("Vtx",     &vtx);
    tree->SetBranchAddress("Run",     &runno_);
    tree->SetBranchAddress("NtrkOff",&ntrkval);
    tree->SetBranchAddress("bin",    &bin);
    tree->SetBranchAddress("ws", &wsv);
    tree->SetBranchAddress("wc", &wcv);
    tree->SetBranchAddress("ws_no_w", &wsv_no_w);
    tree->SetBranchAddress("wc_no_w", &wcv_no_w);
    tree->SetBranchAddress("pt", &pt);
    tree->SetBranchAddress("pt2",&pt2);
    tree->SetBranchAddress("cnt",&ptcnt);
    
    nentries = tree->GetEntries();
    cout<<"Loop 1: "<<tf->GetName()<<"\t"<<nentries<<"\t"<<totentries<<endl;
    if(test) nentries = 1000;
    Loop1();
  }
  
  for(int i = 0; i<NumEPNames; i++) {
    for(int j = 0; j<flat[i]->getHBins(); j++) {
      flat[i]->setXDB(j, (flat[i]->getCnt(j)>0)? flat[i]->getX(j)/flat[i]->getCnt(j): 0);
      flat[i]->setYDB(j, (flat[i]->getCnt(j)>0)? flat[i]->getY(j)/flat[i]->getCnt(j): 0);
      flatOffset[i]->setXDB(j, (flat[i]->getCnt(j)>0)? flat[i]->getX(j)/flat[i]->getCnt(j): 0);
      flatOffset[i]->setYDB(j, (flat[i]->getCnt(j)>0)? flat[i]->getY(j)/flat[i]->getCnt(j): 0);
    }
    for(int j = 0; j<flat[i]->getOBins(); j++) {
      flat[i]->setXoffDB(j, (flat[i]->getXYoffcnt(j)>0)?flat[i]->getXoff(j)/flat[i]->getXYoffcnt(j):0);
      flat[i]->setYoffDB(j, (flat[i]->getXYoffcnt(j)>0)?flat[i]->getYoff(j)/flat[i]->getXYoffcnt(j):0);
      flatOffset[i]->setXoffDB(j, (flat[i]->getXYoffcnt(j)>0)?flat[i]->getXoff(j)/flat[i]->getXYoffcnt(j):0);
      flatOffset[i]->setYoffDB(j, (flat[i]->getXYoffcnt(j)>0)?flat[i]->getYoff(j)/flat[i]->getXYoffcnt(j):0);
    }    
  }
  fclose(save);
  //
  // Loop 2
  //
  Loop2();
  for(int i = 0; i<NumEPNames; i++) {
    for(int j = 0; j<flat[i]->getHBins(); j++) {
      flatOffset[i]->setXDB(j, (flatOffset[i]->getCnt(j)>0)? flatOffset[i]->getX(j)/flatOffset[i]->getCnt(j): 0);
      flatOffset[i]->setYDB(j, (flatOffset[i]->getCnt(j)>0)? flatOffset[i]->getY(j)/flatOffset[i]->getCnt(j): 0);  
    }
    for(int j = 0; j<flat[i]->getOBins(); j++) {
      flatOffset[i]->setXoffDB(j, (flat[i]->getXYoffcnt(j)>0)?flat[i]->getXoff(j)/flat[i]->getXYoffcnt(j):0);
      flatOffset[i]->setYoffDB(j, (flat[i]->getXYoffcnt(j)>0)?flat[i]->getYoff(j)/flat[i]->getXYoffcnt(j):0);
    }
  }
  fclose(save);
  //
  // Loop 3
  //
  Loop3();  
  cout<<"finish loops"<<endl;
  for(int i = 0; i<NumEPNames; i++) {
   
    for(int j = 0; j<flatOffset[i]->getHBins();j++) {
      flatXhist[i]->SetBinContent(j+1,flatOffset[i]->getX(j));
      flatYhist[i]->SetBinContent(j+1,flatOffset[i]->getY(j));
      flatCnthist[i]->SetBinContent(j+1,flatOffset[i]->getCnt(j));
    }
    for(int j = 0; j<flatOffset[i]->getOBins();j++) {
      flatXhist[i]->SetBinContent(j+1+Hbins,flat[i]->getXoff(j));
      flatYhist[i]->SetBinContent(j+1+Hbins,flat[i]->getYoff(j));
      flatCnthist[i]->SetBinContent(j+1+Hbins,flat[i]->getXYoffcnt(j));
      xoffhist[i]->SetBinContent(j+1,flat[i]->getXoff(j));
      yoffhist[i]->SetBinContent(j+1,flat[i]->getYoff(j));
      xyoffcnthist[i]->SetBinContent(j+1,flat[i]->getXYoffcnt(j));
      xyoffmulthist[i]->SetBinContent(j+1,flat[i]->getXYoffmult(j));
      
    }
    for(int j = 0; j<flatOffset[i]->getOBins();j++) {
      flatXhist[i]->SetBinContent(j+1+Hbins+Obins,pt_[i][j]);
      flatYhist[i]->SetBinContent(j+1+Hbins+Obins,pt2_[i][j]);
      flatCnthist[i]->SetBinContent(j+1+Hbins+Obins,ptcnt_[i][j]);
      pthist[i]->SetBinContent(j+1,pt_[i][j]);
      pt2hist[i]->SetBinContent(j+1,pt2_[i][j]);
      ptcnthist[i]->SetBinContent(j+1,ptcnt_[i][j]);
    }
    dirs[i]->cd();
    for(int j = 0; j<10; j++) {
      hpt[i][j]->Write();
      hscale[i][j]->Write();
    }
    hPsi[i]->Write();
    hPsiOffset[i]->Write();
    hPsiOffset2[i]->Write();
    hPsiFlat[i]->Write();
    flatXhist[i]->Write();
    flatYhist[i]->Write();
    flatCnthist[i]->Write();
    xoffhist[i]->Write();
    yoffhist[i]->Write();
    xyoffcnthist[i]->Write();
    pthist[i]->Write();
    pt2hist[i]->Write();
    ptcnthist[i]->Write();
  }
  tfout->cd();
  EPtree->Write();
  hcent->Write();
  hcentbin->Write();
  hvtx->Write();

  tfout->Close(); 
  rescor();  
}
