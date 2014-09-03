/* Simple xAna analysis example. */

#include <vector>
#include <iostream>
#include <TH1D.h>
#include <TRandom.h>
#include <TLorentzVector.h>
#include <TFile.h>
#include <cmath>
#include "puweicalc.h"
#include "untuplizer.h"
#include "ElectronSelections.h"
//#include "MuonSelections.h"
//#include "MuonSelections2.h"
#include "MuonSelectionspt.h"
#include "MuonSelections2pt.h"
#include "PhotonSelections.h"
#include "fit.h"
#include <TSystem.h>
void xAnb() {

  // get TTree from file ...
  TreeReader data("../job_summer12_DYJetsToLL_muon.root"); // v5.3.12

  // ... or make a chain of root files with TTrees
  // const char* paths[] = {"ggtree_mc_1.root", "ggtree_mc_2.root"};
  // TreeReader data(paths, 2);
  
  // useful to determine which type of variable to use for which branches
  // data.Print();
  
  // do whathever preparations are necessary for, if MC information is present
  // pileup reweighting for MC                                                                     
  /*
  unused                     
  Float_t puweiEle, puweiMu;
  */
                                
  PUWeightCalculator puCalcEle;
  PUWeightCalculator puCalcMu;
  if (data.HasMC()) {
    puCalcEle.Init("mcweights/mcwei_PU_RD1_ele.root");
    puCalcMu.Init("mcweights/mcwei_PU_RD1_muo.root");
  }
  
  //enter the interval,ptm
  int interval;
  float bin;
  float ptm;
  // char pary;
  cout<<"Please enter the number of intervals"<<endl;
  cin >> interval;
  cout<<"Please enter the number of pt maximum(approximately 100)"<<endl;
  cin >> ptm;
  // cout<<"Whether using parameter limit(please cin y/Y)"<<endl;
  //cin >>pary;
  bin=(ptm-20)/interval;
  TH1D*  MyHistoMutt[interval];
  TH1D*  MyHistoMutp[interval];
  
  //template
  TH1D* hM = new TH1D("hM", "Two-lepton invariant mass", 90, 50, 140);
  //TH1D* hM_mu = (TH1D*)hM->Clone("hM_mu");
  // TH1D* hM_mu2 = (TH1D*)hM->Clone("hM_mu2");
  // TH1D* hM_el = (TH1D*)hM->Clone("hM_el");
  for(int k=0;k<interval;k++){

    MyHistoMutt[k] = (TH1D*)hM->Clone(Form("MyHistoMutt%d",k));
    MyHistoMutp[k] = (TH1D*)hM->Clone(Form("MyHistoMutp%d",k));
    MyHistoMutt[k]-> SetName(Form("MyHistoMutt%d",k+1));
    MyHistoMutp[k]-> SetName(Form("MyHistoMutp%d",k+1));
   
 }

  //fitting 
 double  ef[interval];
 double  err[interval];
 double  sumtt[interval];
 double  errtt[interval];
 double  sumtp[interval];
 double  errtp[interval];
  //pt loop
 for(int k=0;k<interval;k++){
   float x,y;
  
   x=k*bin+20;
   y=x+bin;



  // event loop
  // for (Long64_t ev = 0; ev < 10000; ev++) {
  for (Long64_t ev = 0; ev < data.GetEntriesFast(); ev++) {
    // print progress
    if (ev % 50000 == 0)
      fprintf(stderr, "Processing event %lli of %lli\n", ev + 1, data.GetEntriesFast());
    
    data.GetEntry(ev);

    Int_t run = data.GetInt("run");
    Long64_t event=data.GetLong64("event");
    // IMPORTANT: branches with counters must be loaded BEFORE dependent branches
    // accessing HLT information
    //    Int_t nHLT = data.GetInt("nHLT");
    //    Int_t* HLTIndex = data.GetPtrInt("HLTIndex");

    //for (Int_t i = 0; i < nHLT; ++i) printf(" %i", HLTIndex[i]);
    //printf("\n");
    
    // IMPORTANT: branches with counters must be loaded BEFORE dependent branches
    Float_t* elePt    = data.GetPtrFloat("elePt");
    Float_t* eleEta   = data.GetPtrFloat("eleEta");
    Float_t* elePhi   = data.GetPtrFloat("elePhi");


    // electron selection
    /*
    vector<int> acc_ele;
    eID2012(data, acc_ele, 1);       // cut based eID (1: loose WP)

    // two-lepton loop
    for (size_t i = 0; i < acc_ele.size(); i++) {
   
      for (size_t j = i + 1; j < acc_ele.size(); j++) {
	TLorentzVector ele1, ele2;
	ele1.SetPtEtaPhiM(elePt[acc_ele[i]], eleEta[acc_ele[i]], elePhi[acc_ele[i]], 0.000511);
	ele2.SetPtEtaPhiM(elePt[acc_ele[j]], eleEta[acc_ele[j]], elePhi[acc_ele[j]], 0.000511);
	
	// fill histo
	TLorentzVector Z = ele1 + ele2;
	hM_el->Fill(Z.M());
      }
    }
*/
    // muon selection

    Float_t* muPt    = data.GetPtrFloat("muPt");
    Float_t* muEta   = data.GetPtrFloat("muEta");
    Float_t* muPhi   = data.GetPtrFloat("muPhi");

    vector<int> acc_mu;    
    
    select_muons(data, acc_mu,x,y); // tight muon ID

   

    for (size_t i = 0; i < acc_mu.size(); i++) {
      for (size_t j = i + 1; j < acc_mu.size(); j++) {
	TLorentzVector mu1, mu2;
	mu1.SetPtEtaPhiM(muPt[acc_mu[i]], muEta[acc_mu[i]], muPhi[acc_mu[i]], 0.1057);
	mu2.SetPtEtaPhiM(muPt[acc_mu[j]], muEta[acc_mu[j]], muPhi[acc_mu[j]], 0.1057);

	TLorentzVector  Z = mu1 + mu2;
        MyHistoMutt[k]->Fill(Z.M());
        if((k==5)&&(Z.M()<130)&&(Z.M()>120)){
        cout<<"run="<<run<<endl<<"event="<<event<<endl;
        //cin.get();
        gSystem->Sleep(2000);
	}
       }
    }
 
    //Tag and Probe


    vector<int> acc_muProbe;    
    select_muonsProbe(data,acc_muProbe,x,y); // loose  muon ID

    for (size_t i = 0; i < acc_mu.size(); i++) {
      for (size_t j = 0; j < acc_muProbe.size(); j++) {
	if(acc_mu[i]>=acc_muProbe[j]) continue;
        TLorentzVector mu3, mu4;
	mu3.SetPtEtaPhiM(muPt[acc_mu[i]], muEta[acc_mu[i]], muPhi[acc_mu[i]], 0.1057);
	mu4.SetPtEtaPhiM(muPt[acc_muProbe[j]], muEta[acc_muProbe[j]], muPhi[acc_muProbe[j]], 0.1057);

	TLorentzVector  Z2 = mu3 + mu4;
        MyHistoMutp[k]->Fill(Z2.M());
      }
    }
  

  }
  
  //double par0 =data.GetEntriesFast();
  //if ((pary=='y')||(pary=='Y')) par0 =data.GetEntriesFast()/interval;
 
  MyHistoMutp[k]-> Add( MyHistoMutt[k],-1); 
  FitZMass(k,MyHistoMutt[k],sumtt[k],errtt[k]);
  FitZMass(k,MyHistoMutp[k],sumtp[k],errtp[k]);

 // event loop
 // cout<<"x="<<x<<endl;
 // cout<<"y="<<y<<endl;
}
 //pt loop 
  fprintf(stderr, "Processed all events\n");
  

 for(int k=0;k<interval;k++){
   double i1=sumtt[k],i2=sumtp[k],j1=errtt[k],j2=errtp[k];


  ef[k]=2*i1/(2*i1+i2);
  err[k]=sqrt((2*i2*j1/(2*i1+i2)/(2*i1+i2))*(2*i2*j1/(2*i1+i2)/(2*i1+i2))
     +(2*i1*j2/(2*i1+i2)/(2*i1+i2))*(2*i1*j2/(2*i1+i2)/(2*i1+i2)));
 }


  TFile* outFile = new TFile("zmass_emu.root","recreate");
  //hM_el->Write();
  //hM_mu->Write();
  //hM_mu2->Write();

  for(int k=0;k<interval;k++){
  MyHistoMutt[k]->Write();

  MyHistoMutp[k]->Write();
   
 }

  outFile->Write();
  //outFile->Close();

  //fit
  TH1D* hef = new TH1D("ef", "The Tag and Probe Efficiency of Different Pt", interval, 20,ptm );
  //TH1D* heff = new TH1D("efferr<3", "(Fixed) The Tag and Probe Efficiency of Different Pt", interval, 20,ptm );
  //TH1D* hefff = new TH1D("effHaveErrExclude", "(Fixed twice) The Tag and Probe Efficiency of Different Pt", interval, 20,ptm );
  // for(int k=0;k<interval;k++){
 
  //cout<<"ef="<<ef<<endl;
  // } 
 
  for(int k=0;k<interval;k++) cout<<Form("sumtt[%d] = ",k)<<sumtt[k]<<endl;  
  for(int k=0;k<interval;k++)  cout<<Form("sumtp[%d] = ",k)<<sumtp[k]<<endl;
  for(int k=0;k<interval;k++) cout<<Form("errtt[%d] = ",k)<<errtt[k]<<endl;  
  for(int k=0;k<interval;k++)  cout<<Form("errtp[%d] = ",k)<<errtp[k]<<endl;
  for(int k=0;k<interval;k++) cout<<Form("ef[%d] = ",k)<<ef[k]<<endl;
  for(int k=0;k<interval;k++) cout<<Form("err[%d] = ",k)<<err[k]<<endl;
   
  
  
  //fix
  /*  for(int k=1;k+1<interval;k++){
    int n=0;
    if(ef[k]=0){
      float m=0;
      
      while(m=0){
	m=ef[k+1];
	n++;
      }
       cout<<Form("ef[%d] = ",k-1)<<ef[k-1]<<endl;
       cout<<Form("ef[%d] = ",k+n)<<ef[k+n]<<endl;
       cout<<Form("pow(2,%d) ",n)<<pow(2,n)<<endl;



      //ef[k]=ef[k-1]/2+ef[k+n]/pow(2,n);

     }
  }
  
  */
  // for(int k=0;k<interval;k++){
  //   cout<<Form("ef[%d] = ",k)<<ef[k]<<endl;
  // }
 

  for(int k=0;k<interval;k++){
  hef->SetBinContent(k+1,ef[k]);
  hef->SetBinError(k+1,err[k]);
  }

  /*
  for(int k=0;k<interval;k++){
    heff->SetBinContent(k+1,ef[k]);
    heff->SetBinError(k+1,err[k]);
    if((ef[k]<0)) heff->SetBinContent(k+1,0);
    if((err[k]>3))heff->SetBinError(k+1,3);
  }


   for(int k=0;k<interval;k++){
    hefff->SetBinContent(k+1,ef[k]);
    hefff->SetBinError(k+1,err[k]);
    if((ef[k]<0)) hefff->SetBinContent(k+1,0);
    if((err[k]>3)){
      hefff->SetBinError(k+1,0);
      hefff->SetBinContent(k+1,0); 
   }
  }


  */
 // TFile* outFileEf = new TFile("The Tag and Probe Efficiency of Different Pt","recreate");
  hef->Write();
  // heff->Write();
  //hefff->Write();
  //outFileEf->Write();

  cout<<"interval ="<<interval<<endl;
  cout<<"bin ="<<bin<<endl;
}


