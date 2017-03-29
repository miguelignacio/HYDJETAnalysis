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

#include <vector>
#include <iostream>
#include <fstream>
#include <stdio.h>
using namespace std;


void SubtractPedestal(TH1F &h1){

  double min = h1.GetMinimum();
  std::cout << " minimum " << min << std::endl;
  for(int n = 0; n< h1.GetNbinsX()+1; n++){
    h1.SetBinContent(n, h1.GetBinContent(n)-min);
  }

  return ;
}

void Analysis(){

  std::cout<<"Starting Analysis" << std::endl;
  gStyle->SetOptStat("");

  //auto myFile = TFile::Open("MC_DATA/Cen010_1kevents.root ","READ");
  auto myFile = TFile::Open("MC_DATA/Cen3040_30kevents.root","READ");

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

  double eta_max = 1.9; //maximum dEta
  double phi_min = -0.5; // rads
  double phi_max = 1.5; // rads
  double n_etabins = 40;
  double n_phibins = 40;

  auto h_Correlation_Raw = new TH2F("h_Correlation_Raw","",
				n_phibins, phi_min, phi_max,
				n_etabins,-eta_max,eta_max);
  //These are corrected.
  auto h_Correlation_All = new TH2F("h_Correlation_All","", 
				     n_phibins, phi_min, phi_max,
				     n_etabins,-eta_max,eta_max);
  auto h_Correlation_Signal = new TH2F("h_Correlation_Signal","",
				       n_phibins, phi_min, phi_max,
				       n_etabins,-eta_max,eta_max);

  auto h_Correlation_Background = new TH2F("h_Correlation_Background","",
                                       n_phibins, phi_min, phi_max,
                                       n_etabins,-eta_max,eta_max);
  //projections:
  auto h_1D_phi_All = new TH1F("h_1D_phi_All", "", 25, phi_min, phi_max); //units of pi
  auto h_1D_phi_Signal = new TH1F("h_1D_phi_Signal", "", 25, phi_min, phi_max); //units of pi
  auto h_1D_phi_Background = new TH1F("h_1D_phi_Background", "", 25, phi_min, phi_max); //units of pi

  auto h_1D_phi_LargeEta_All = new TH1F("h_1D_phi_LargeEta_All", "", 25, phi_min, phi_max); //units of pi
  auto h_1D_phi_LargeEta_Signal = new TH1F("h_1D_phi_LargeEta_Signal", "", 25, phi_min, phi_max); //units of pi
  auto h_1D_phi_LargeEta_Background = new TH1F("h_1D_phi_LargeEta_Background", "", 25, phi_min, phi_max); //units of pi
  
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
  //looping over events

  int nevent = 0;
  
  while (myReader.Next()) {
    //    std::cout<< " NTot" << *Ntot << " Npart " << *Npart << " Nbcol " << *Nbcol << std::endl;
    nevent +=1;
    if(nevent%1000==0) std::cout << " Event # " << nevent << std::endl;
    // if (nevent>100) break;
    h_Npart->Fill(*Npart);
    h_Nbcol->Fill(*Nbcol);

    // std::cout << Pt->size() << " " << Eta->size() << " " << Phi->size () << std::endl;
    //       std::cout << final->size() << " " << tipo->size() << " " << pythiaStatus->size() << std::endl;
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
    if(nevent%500==0) std::cout << " Event # " << nevent << std::endl;
    for(int i = 0; i<Pt->size(); i++){
      if( abs(pdg->at(i))!=211 and abs(pdg->at(i))!=321 and abs(pdg->at(i))!=2212) continue; //select only pions, kaons and protons. 
      if(final->at(i)==0) continue; //skip if not final state
      if(Pt->at(i)>5.0 and Pt->at(i)<10.0 and abs(Eta->at(i))<1.0 ) {
        double trigger_phi =Phi->at(i);
        double trigger_eta =Eta->at(i);
	//looping over associated particle
        for(int j=0; j<Pt->size(); j++){
	  if( abs(pdg->at(j))!=211 and abs(pdg->at(j))!=321 and abs(pdg->at(j))!=2212) continue; //select only pions kaons and protons
	  if(final->at(j)==0) continue; //skip if not final state
          if(Pt->at(j)>1.0 and Pt->at(j)<2.0 and abs(Eta->at(j))<1.0 ){
            double assoc_phi = Phi->at(j);
            double assoc_eta = Eta->at(j);
            double dphi = TVector2::Phi_mpi_pi(trigger_phi-assoc_phi);
            double deta = trigger_eta-assoc_eta;
	    double weight =1;
	    if(abs(deta)>=1.9) continue;
	    if(deta<0) weight = 2+deta;
	    else weight = 2-deta;
	    dphi = dphi/TMath::Pi();
            if(dphi<-0.5) dphi +=2;
	    //fprintf (fData, "%f \n", dphi);	    

            //Filling histograms
	    h_Correlation_Raw->Fill(dphi,deta);
            h_Correlation_All->Fill(dphi,deta, 1.0/weight);
	    
	    if(type->at(i)>0 and type->at(j)>0) h_Correlation_Signal->Fill(dphi,deta, 1.0/weight);
	    else if(type->at(i)==0 or type->at(j)==0) h_Correlation_Background->Fill(dphi,deta,1.0/weight);
	    //projecting into eta<0.5 (maybe this can be done in other place).
	    if(abs(deta)<0.5){
                h_1D_phi_All->Fill(dphi, 1.0/weight);
		if(type->at(i)>0 and type->at(j)>0) h_1D_phi_Signal->Fill(dphi, 1.0/weight);
		else if(type->at(i)==0 or type->at(j)==0) h_1D_phi_Background->Fill(dphi,1.0/weight);
	    }
	    else if(abs(deta)>1.0){
	      h_1D_phi_LargeEta_All->Fill(dphi, 1.0/weight);
	      if(type->at(i)>0 and type->at(j)>0) h_1D_phi_LargeEta_Signal->Fill(dphi, 1.0/weight);
	      else if(type->at(i)==0 or type->at(j)==0) h_1D_phi_LargeEta_Background->Fill(dphi,1.0/weight);
	    }
	    //else if eta is large
         }
        } //finish loop over associated
      }
    }// finish loop over triggers
  } // finish loop over events
  
  auto fout = new TFile("fout_histos.root","RECREATE");
  h_Correlation_All->Write("Correlation_All");
  h_Correlation_Signal->Write("Correlation_Signal");
  h_Correlation_Background->Write("Correlation_Background");

  h_Correlation_Raw->Write("Correlation_Raw");
  h_1D_phi_All->Write("Dphi_All");
  h_1D_phi_Signal->Write("Dphi_Signal");
  h_1D_phi_Background->Write("Dphi_Background");

  h_1D_phi_LargeEta_All->Write("Dphi_LargeEta_All");
  h_1D_phi_LargeEta_Signal->Write("Dphi_LargeEta_Signal");
  h_1D_phi_LargeEta_Background->Write("Dphi_LargeEta_Background");
  
  fout->Close();
  return;
}
