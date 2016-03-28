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
void Loop2(){
  cout<<"Loop 2:"<<endl;
  save = fopen("save","rb");
  ncnt = 0;
  while(fread(&sout,sizeof(struct sout_struct),1,save)>0){
    vtx = sout.vtx;
    centval = sout.cent;
    bin = sout.bin;
    for(int j = 0; j<NumEPNames; j++) {
      double psiOffset = -10;
      double psi = -10;
      if(sout.msum[j]>0) {
	double order = EPOrder[j];
	psi = atan2(sout.ws[j],sout.wc[j])/order;
	psiOffset = flat[j]->getOffsetPsi(flat[j]->getSoffset(sout.ws[j],vtx,bin),flat[j]->getCoffset(sout.wc[j],vtx,bin));
	if(centval<80 && fabs(vtx)<15&&bin>=0) {
	  hPsiOffset[j]->Fill(psiOffset);	
	}
	flatOffset[j]->fill(psiOffset,vtx,bin);
	flatOffset[j]->fillOffset(sout.ws[j],sout.wc[j],sout.msum[j],vtx,bin);
      }
      //if(ncnt<10&&j==5) cout<<psi<<"\t"<<psiOffset<<endl; 
    } 
    ++ncnt;
  }

}
