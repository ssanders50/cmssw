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
void Loop1(){
    for(int i = 0; i<nentries; i++) {
      tree->GetEntry(i);
      hcent->Fill(centval);
      hvtx->Fill(vtx);

      sout.vtx = vtx;
      sout.bin = bin;
      sout.cent = centval;
      sout.run = runno_;
      sout.ntrk = ntrkval;
      ++totentries;
      for(int j = 0; j<NumEPNames; j++) {
	int indx = flat[j]->getOffsetIndx(bin,vtx);
  	sout.ws[j] = 0;
  	sout.wc[j] = 0;
  	sout.msum[j] = 0;
  	sout.wsum[j] = 0;
  	double order = EPOrder[j];
  	if(wsv[j]==0&&wcv[j]==0) continue;
	double scale = flat[j]->getEtScale(vtx,bin);
	double s = wsv[j]*scale;
	double c = wcv[j]*scale;
	double snow = wsv_no_w[j]*scale;
	double cnow = wcv_no_w[j]*scale;
	if(MomConsWeight[j][0]=='y' && ptav[j]>0) {
	  double sgn = 1.;
	  if(EPEtaMax1[j]<0) sgn = -1.;
	  s=s-sgn*(pt2_[j][indx]/pt_[j][indx])*snow;
	  c=c-sgn*(pt2_[j][indx]/pt_[j][indx])*cnow;
	}
  	double psi = atan2(s,c)/order;
  	flat[j]->fill(psi,vtx,bin);
  	flat[j]->fillOffset(s,c,ptcnt[j],vtx,bin);
  	sout.ws[j] = s;
  	sout.wc[j] = c;
  	sout.wsum[j] = ptcnt[j];
  	sout.msum[j] = ptcnt[j];
  	if(centval<80 && fabs(vtx)<15&&bin>=0) hPsi[j]->Fill(psi);
	int k = centval/10;
	if(ptcnt[j]>0) {
	  hpt[j][k]->Fill(pt[j]*scale/ptcnt[j]);
	  hscale[j][k]->Fill(scale);
	}
      }
      ++ncnt;
      fwrite(&sout,sizeof(struct sout_struct),1,save);
    }

}
