#include "TFile.h"
#include "TH2F.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TF1.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TMath.h"
#include "TVector2.h"

#include "THnSparse.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <stdio.h>
using namespace std;


THnSparseD* h;
THnSparseD* h_trigger;

TH2F* h_Correlation_Raw;
TH2F* h_Correlation_Corrected;
void SubtractPedestal(TH1F &h1){

  double min = h1.GetMinimum();
  std::cout << " minimum " << min << std::endl;
  for(int n = 0; n< h1.GetNbinsX()+1; n++){
    h1.SetBinContent(n, h1.GetBinContent(n)-min);
  }

  return ;
}

void FillHistograms(float dphi, float deta, float phi_trigger, Bool_t signal, Double_t pT_assoc, Double_t pT_trigger){
 //this function will fill all histograms
  //Calculating acceptance correction.
  double weight = 1.0;
  weight = 2-abs(deta);
  h_Correlation_Raw->Fill(dphi,deta);
  h_Correlation_Corrected->Fill(dphi,deta,1.0/weight);
  Double_t temp = 0.0; //background
  if(signal) temp=1.0; //signal

  phi_trigger = abs(phi_trigger)/TMath::Pi();
  if(phi_trigger>0.5) phi_trigger = 1-phi_trigger;

  Double_t entries[6] = {dphi, abs(deta), phi_trigger, temp, pT_assoc, pT_trigger};
  h->Fill( entries, 1.0/weight);
  return;
}


void Run(TString namefile){

  std::cout<<"Starting Analysis" << std::endl;
  auto myFile =  TFile::Open("MC_DATA/"+namefile,"READ");

  if (!myFile || myFile->IsZombie()) {
    return;
  }
  myFile->Print();

  auto tout =(TTree*)myFile->Get("td");
  if(!tout) {
    std::cout <<"Problem getting TTree" << std::endl;
    return;
  }
  //output histograms
  auto h_pT = new TH1F("h_pT","",100,0,20);
  h_pT->SetTitle("; pT; Entries");
  auto h_Eta = new TH1F("h_Eta", "" , 50, -1.0,1.0);
  h_Eta->SetTitle("; #eta ; Entries");
  auto h_Npart = new TH1F("h_Npart", "", 100, 0, 500);
  h_Npart->SetTitle("; Npart; Events");
  auto h_Nbcol = new TH1F("h_Nbcol", "", 100, 0, 2000);
  h_Nbcol->SetTitle("; Nbcol; Events");

  auto h_Phi = new TH1F("h_Phi", "", 100, -TMath::Pi(), TMath::Pi());
  h_Phi->SetTitle("; #phi ; Entries");

  auto h_type = new TH1F("h_type", "", 2, -0.5,1.5);
  h_type->SetTitle("; type ; Entries");
    
  auto h_pythiaStatus = new TH1F("h_pythiaStatus", "", 2, -0.5,1.5);
  h_pythiaStatus->SetTitle("; pythiaStatus ; Entries");

  auto h_final = new TH1F("h_final", "", 2, -0.5,1.5);
  h_final->SetTitle("; final ; Entries");

  Float_t eta_max = 2.0; //maximum dEta
  Float_t phi_min = -0.5; // rads
  Float_t phi_max = 1.5; // rads
  Int_t n_etabins = 20;
  Int_t n_phibins = 18;

  h_Correlation_Raw = new TH2F("h_Correlation_Raw","",
				n_phibins, phi_min, phi_max,
				n_etabins,-1.9,1.9);

  h_Correlation_Corrected = new TH2F("h_Correlation_Corrected","",
			       n_phibins, phi_min, phi_max,
			       n_etabins,-1.9,1.9);

  h_Correlation_Raw->SetTitle("; #Delta#phi ; #Delta#eta");
  h_Correlation_Corrected->SetTitle("; #Delta#phi ; #Delta#eta");


  //axes are Dphi, Deta, trigger_phi, signal/background, pt_assoc

  Float_t min_pTassoc =0;
  Float_t max_pTassoc =10;
  Int_t   n_bins_pTassoc = 20;
 
  Float_t min_pTtrigger =0;
  Float_t max_pTtrigger =20;
  Int_t   n_bins_pTtrigger = 20;

  Float_t min_trigphi = 0;
  Float_t max_trigphi = 0.5; // phi/2
  Int_t n_bins_trigphi = 3;

  Int_t bins[6]    = {n_phibins, n_etabins, n_bins_trigphi, 2,     n_bins_pTassoc,  n_bins_pTtrigger };
  Double_t xmin[6] = {phi_min,   0,         min_trigphi   ,  -0.5, min_pTassoc,     min_pTtrigger };
  Double_t xmax[6] = {phi_max,   eta_max,   max_trigphi   , 1.5,   max_pTassoc,     max_pTtrigger};
  h = new THnSparseD("h"," h ; #Delta#varphi; #Delta#eta; a; pT_assoc; pT_trigger",6,bins,xmin,xmax);
 
  //axes are trig_pT, trig_phi; signal or not 
  Int_t bins_trig[2]    = {n_bins_pTtrigger, n_bins_trigphi};
  Double_t xmin_trig[2] = {min_pTtrigger,    min_trigphi};
  Double_t xmax_trig[2] = {max_pTtrigger,    max_trigphi};

  h_trigger = new THnSparseD("h_trigger", "#pT; #phi",2,bins_trig,xmin_trig,xmax_trig); 
  
  //Reader 
  TTreeReader myReader(tout);
  TTreeReaderValue<std::vector<float>> Pt(myReader, "Pt");
  TTreeReaderValue<std::vector<float>> Eta(myReader, "Eta");
  TTreeReaderValue<std::vector<float>> Phi(myReader, "Phi");
  TTreeReaderValue<Int_t> Ntot(myReader, "Ntot");
  TTreeReaderValue<Int_t> Npart(myReader, "Npart");
  TTreeReaderValue<Int_t> Nbcol(myReader, "Nbcol");
  TTreeReaderValue<std::vector<int>> pythiaStatus(myReader, "pythiaStatus");
  TTreeReaderValue<std::vector<int>> type(myReader, "type");
  TTreeReaderValue<std::vector<int>> final(myReader, "final");
  TTreeReaderValue<std::vector<int>> pdg(myReader, "pdg");
  TTreeReaderValue<std::vector<int>> Mpdg(myReader, "Mpdg");
  //looping over events

  int nevent = 0;
  
  while (myReader.Next()) {
    nevent +=1;
    if(nevent%1000==0) std::cout << " Event # " << nevent << std::endl;
    h_Npart->Fill(*Npart);
    h_Nbcol->Fill(*Nbcol);

    for(int i = 0; i<Pt->size(); i++){
      // std::cout << Eta->at(i) << std::endl;
      if(abs(Eta->at(i))>1.0) continue;
      if(Pt->at(i)<0.2) continue;
            
      h_pT->Fill(Pt->at(i));
      h_Eta->Fill(Eta->at(i));
      h_Phi->Fill(Phi->at(i));
      h_final->Fill(final->at(i));
      h_type->Fill(type->at(i));
      h_pythiaStatus->Fill(pythiaStatus->at(i));
    }
  } //end loop of events

  myReader.Restart();
  nevent =0;

  while (myReader.Next()) {
    nevent +=1;
    if(nevent%100==0) std::cout << " Event # " << nevent << std::endl;
    //if(nevent>1000) break;
    for(int i = 0; i<Pt->size(); i++){
      //if( abs(pdg->at(i))!=211 and abs(pdg->at(i))!=321 and abs(pdg->at(i))!=2212) continue; //select only pions, kaons and protons. 
      //if (pdg->at(i) != 22) continue; //if not photon, not use
      if(abs(pdg->at(i))!=111 and abs(pdg->at(i))!=211) continue;
      //if (Mpdg->at(i)!=-1) continue; //if decay, not use
      //if(final->at(i)==0) continue; //it not final state, not use
      if(Pt->at(i)>3.0 and abs(Eta->at(i))<1.0){

        double trigger_phi =Phi->at(i);
        double trigger_eta =Eta->at(i);

	double phi_trigger = abs(trigger_phi)/TMath::Pi();
	if(phi_trigger>0.5) phi_trigger = 1-phi_trigger;

	double pt_trigger = Pt->at(i);
	//looping over associated particle
	double entries[2] = {pt_trigger, phi_trigger};
	h_trigger->Fill(entries);

        for(int j=0; j<Pt->size(); j++){
	  if( abs(pdg->at(j))!=211 and abs(pdg->at(j))!=321 and abs(pdg->at(j))!=2212) continue; //select only pions kaons and protons
	  if(final->at(j)==0) continue; //skip if not final state
	  if(abs(Eta->at(j))<1.0){
	    double assoc_phi = Phi->at(j);
            double assoc_eta = Eta->at(j);
            double dphi = TVector2::Phi_mpi_pi(trigger_phi-assoc_phi);
            double deta = trigger_eta-assoc_eta;
	    if(abs(deta)>=1.9) continue;
	    dphi = dphi/TMath::Pi();
            if(dphi<-0.5) dphi +=2;

            //Filling histograms
            Bool_t signal = kFALSE;
	    if(type->at(i)>0 and type->at(j)>0) signal = kTRUE;
	    
	    double pt_assoc = Pt->at(j);

            FillHistograms(dphi, deta, trigger_phi, signal, pt_assoc, pt_trigger);
	      
	    //else if eta is large
	  } 
        } //finish loop over associated
      }
    }// finish loop over triggers
  } // finish loop over events
  
  auto fout = new TFile("ROOTFILES/fout_"+namefile,"RECREATE");
  h_Correlation_Raw->Write("Correlation_Raw");
  h_Correlation_Corrected->Write("Correlation_Corrected");
  h->Write("Sparse");
  h_trigger->Write("Sparse_trigger");
  
  fout->Close();
  return;
}


void RunOverMC(){
  // Run("Cen3040_40k_Baseline.root");
  Run("Cen3040_20k_Baseline.root");
  //Run("Cen3040_10k_baseline.root");
  // Run("Cen3040_5k_quenched.root");
  //Run("Cen020_2k_Baseline.root");
  //Run("Cen020_2k_Quenched.root");
  return;
}
