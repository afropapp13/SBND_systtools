#define truepreselection_cxx
#include "truepreselection.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TString.h>

#include <iostream>   // std::cout
#include <string>     // std::string, std::to_string

#include "Constants.h"

using namespace std;
using namespace Constants;

//----------------------------------------//

void truepreselection::Loop() {

  //----------------------------------------//

   if (fChain == 0) return;
   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;

   //----------------------------------------//

   // Output file

   TString file_name = "true_events_" + fdial + "_" + TString(to_string(funiverse)) + ".root"; 
   TFile* fout = new TFile(file_name,"recreate");

   //----------------------------------------//

   // Plot declaration

   TH1D* EnuPlot[ccnc_label.size()][nu_label.size()][inte_label.size()];
   TH1D* LepPPlot[ccnc_label.size()][nu_label.size()][inte_label.size()];
   TH1D* LepCosPlot[ccnc_label.size()][nu_label.size()][inte_label.size()];
   TH1D* ProtPPlot[ccnc_label.size()][nu_label.size()][inte_label.size()];
   TH1D* ProtCosPlot[ccnc_label.size()][nu_label.size()][inte_label.size()];

   //----------------------------------------//

   // Plot loop over cc/nc, neutrino flavors, interactions

   // Loop over All, CC, NC

   for (int icc = 0; icc < (int)ccnc_label.size(); icc++) {

     // Loop over All, numu, nue, numubar, nuebar

     for (int inu = 0; inu < (int)nu_label.size(); inu++) {

       // Loop over All, QE, MEC, RES, DIS

       for (int iinte = 0; iinte < (int)inte_label.size(); iinte++) {

	 TString enu_plot_name = "Enu_" + ccnc_label[icc] + "_" + nu_label[inu] + "_" + inte_label[iinte];
	 EnuPlot[icc][inu][iinte] = new TH1D(enu_plot_name,";E_{#nu} [GeV];",20,0.,2.);

	 TString lepp_plot_name = "LepP_" + ccnc_label[icc] + "_" + nu_label[inu] + "_" + inte_label[iinte];
	 LepPPlot[icc][inu][iinte] = new TH1D(lepp_plot_name,";p_{lep} [GeV/c];",20,0.,2.);

	 TString lepcos_plot_name = "LepCos_" + ccnc_label[icc] + "_" + nu_label[inu] + "_" + inte_label[iinte];
	 LepCosPlot[icc][inu][iinte] = new TH1D(lepcos_plot_name,";cos#theta_{lep};",20,-1.,1.);

	 TString protp_plot_name = "ProtP_" + ccnc_label[icc] + "_" + nu_label[inu] + "_" + inte_label[iinte];
	 ProtPPlot[icc][inu][iinte] = new TH1D(protp_plot_name,";p_{prot} [GeV/c];",20,0.,2.);

	 TString protcos_plot_name = "ProtCos_" + ccnc_label[icc] + "_" + nu_label[inu] + "_" + inte_label[iinte];
	 ProtCosPlot[icc][inu][iinte] = new TH1D(protcos_plot_name,";E_{#nu} [GeV];",20,0.,2.);

       } // end of the loop over All, QE, MEC, RES, DIS


     } // end of the loop over All, numu, nue, numubar, nuebar


   } // end of the loop over All, CC, NC

   //----------------------------------------//

   // Loop over the events

   for (Long64_t jentry=0; jentry<nentries;jentry++) {

     //----------------------------------------//

      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      if (jentry%1000 == 0) std::cout << jentry/1000 << " k " << std::setprecision(2) << double(jentry)/nentries*100. << " %"<< std::endl;

      //----------------------------------------//

      // Loop over all the neutrinos in the event
      //cout << "mc_nnu = " << mc_nnu << endl;

      for (int inu = 0; inu < mc_nnu; inu++) {

	//----------------------------------------//

	// weights

	double weight = 1.;

	//----------------------------------------//

	// Tag based on current (CC vs NC)
	// and interaction (QE, MEC, RES, DIS, COH, other)

	int cc = true_cc->at(inu);
	//cout << "cc = " << cc << endl;

	int pdg = -1;
	int nu_pdg = true_nu_pdg->at(inu);
	//cout << nu_pdg << endl;

	if (nu_pdg == 14) { pdg = 0; }
	else if (nu_pdg == 12) { pdg = 1; }
	else if (nu_pdg == -14) { pdg = 2; }
	else if (nu_pdg == -12) { pdg = 3; }
	else { cout << "unknown nu pdg, aborting" << endl; continue; }

	int mode = -1;
	int nu_mode = true_mode->at(inu);
	if (nu_mode == 0) { mode = 1; } // QE
	else if (nu_mode == 10) { mode = 2; } // MEC
	else if (nu_mode == 1) { mode = 3; } // RES
	else if (nu_mode == 2) { mode = 4; } // DIS
	else if (nu_mode == 3) { mode = 5; } // COH
	else { mode = 6; } // other
	//cout << "nu_mode = " << nu_mode << endl;

	//----------------------------------------//

	//	int n_mcpart = true_particle_pdg[inu].size();
	int n_mcpart = true_particle_pdg->at(inu).size();
	//cout << "n_mcpart = " << n_mcpart << endl;

	int n_leptons = 0;
	int n_protons = 0;

	int lepton_index = -1;
	double lepton_p = 0.;

	int proton_index = -1;
	double proton_p = 0.;

	for (int imcpart = 0; imcpart < n_mcpart; imcpart++) {
	  
	  int pdg_mcpart = true_particle_pdg->at(inu).at(imcpart);
	  double p_mcpart = true_particle_p->at(inu).at(imcpart);

	  //cout << "pdg_mcpart = " << pdg_mcpart << " p_mcpart = " << p_mcpart << endl;

	  if ( TMath::Abs(pdg_mcpart) == 11 || TMath::Abs(pdg_mcpart) == 12 || TMath::Abs(pdg_mcpart) == 13 || TMath::Abs(pdg_mcpart) == 14) { 
	    
	    n_leptons++; 

	    if (p_mcpart > lepton_p) { 

	      lepton_p = p_mcpart;
	      lepton_index = imcpart;

	    }

	  }

	  if (pdg_mcpart == 2212) { 

	    n_protons++; 

	    if (p_mcpart > proton_p) { 

	      proton_p = p_mcpart;
	      proton_index = imcpart;

	    }

	  }
	  
	}

	//cout << "n_leptons = " << n_leptons << endl;
	//cout << "n_protons = " << n_protons << endl;
	//cout << endl;

	//cout << "lepton_index = " << lepton_index << endl;
	//cout << "proton_index = " << proton_index << endl;

	//----------------------------------------//

	// Each plot needs to be filled 8 eights
	// to capture all the combinations

	EnuPlot[0][0][0]->Fill(true_Ev->at(inu),weight);
	EnuPlot[0][pdg][mode]->Fill(true_Ev->at(inu),weight);
	EnuPlot[cc][0][mode]->Fill(true_Ev->at(inu),weight);
	EnuPlot[cc][pdg][0]->Fill(true_Ev->at(inu),weight);
	EnuPlot[0][0][mode]->Fill(true_Ev->at(inu),weight);
	EnuPlot[0][pdg][0]->Fill(true_Ev->at(inu),weight);
	EnuPlot[cc][0][0]->Fill(true_Ev->at(inu),weight);
	EnuPlot[cc][pdg][mode]->Fill(true_Ev->at(inu),weight);

	LepPPlot[0][0][0]->Fill(true_particle_p->at(inu).at(lepton_index),weight);
	LepPPlot[0][pdg][mode]->Fill(true_particle_p->at(inu).at(lepton_index),weight);
	LepPPlot[cc][0][mode]->Fill(true_particle_p->at(inu).at(lepton_index),weight);
	LepPPlot[cc][pdg][0]->Fill(true_particle_p->at(inu).at(lepton_index),weight);
	LepPPlot[0][0][mode]->Fill(true_particle_p->at(inu).at(lepton_index),weight);
	LepPPlot[0][pdg][0]->Fill(true_particle_p->at(inu).at(lepton_index),weight);
	LepPPlot[cc][0][0]->Fill(true_particle_p->at(inu).at(lepton_index),weight);
	LepPPlot[cc][pdg][mode]->Fill(true_particle_p->at(inu).at(lepton_index),weight);

	LepCosPlot[0][0][0]->Fill(true_particle_costheta->at(inu).at(lepton_index),weight);
	LepCosPlot[0][pdg][mode]->Fill(true_particle_costheta->at(inu).at(lepton_index),weight);
	LepCosPlot[cc][0][mode]->Fill(true_particle_costheta->at(inu).at(lepton_index),weight);
	LepCosPlot[cc][pdg][0]->Fill(true_particle_costheta->at(inu).at(lepton_index),weight);
	LepCosPlot[0][0][mode]->Fill(true_particle_costheta->at(inu).at(lepton_index),weight);
	LepCosPlot[0][pdg][0]->Fill(true_particle_costheta->at(inu).at(lepton_index),weight);
	LepCosPlot[cc][0][0]->Fill(true_particle_costheta->at(inu).at(lepton_index),weight);
	LepCosPlot[cc][pdg][mode]->Fill(true_particle_costheta->at(inu).at(lepton_index),weight);

	if (n_protons > 0) {

	  ProtPPlot[0][0][0]->Fill(true_particle_p->at(inu).at(proton_index),weight);
	  ProtPPlot[0][pdg][mode]->Fill(true_particle_p->at(inu).at(proton_index),weight);
	  ProtPPlot[cc][0][mode]->Fill(true_particle_p->at(inu).at(proton_index),weight);
	  ProtPPlot[cc][pdg][0]->Fill(true_particle_p->at(inu).at(proton_index),weight);
	  ProtPPlot[0][0][mode]->Fill(true_particle_p->at(inu).at(proton_index),weight);
	  ProtPPlot[0][pdg][0]->Fill(true_particle_p->at(inu).at(proton_index),weight);
	  ProtPPlot[cc][0][0]->Fill(true_particle_p->at(inu).at(proton_index),weight);
	  ProtPPlot[cc][pdg][mode]->Fill(true_particle_p->at(inu).at(proton_index),weight);

	  ProtCosPlot[0][0][0]->Fill(true_particle_costheta->at(inu).at(proton_index),weight);
	  ProtCosPlot[0][pdg][mode]->Fill(true_particle_costheta->at(inu).at(proton_index),weight);
	  ProtCosPlot[cc][0][mode]->Fill(true_particle_costheta->at(inu).at(proton_index),weight);
	  ProtCosPlot[cc][pdg][0]->Fill(true_particle_costheta->at(inu).at(proton_index),weight);
	  ProtCosPlot[0][0][mode]->Fill(true_particle_costheta->at(inu).at(proton_index),weight);
	  ProtCosPlot[0][pdg][0]->Fill(true_particle_costheta->at(inu).at(proton_index),weight);
	  ProtCosPlot[cc][0][0]->Fill(true_particle_costheta->at(inu).at(proton_index),weight);
	  ProtCosPlot[cc][pdg][mode]->Fill(true_particle_costheta->at(inu).at(proton_index),weight);

	}

	//----------------------------------------//

      } // end of the loop over the neutrinos in the event

   } // end of the loop over the events

  //----------------------------------------//

   fout->cd();
   fout->Write();
   fout->Close();
   cout << "File " << file_name << " has been created!" << endl;

   //----------------------------------------//

} // end of teh program
