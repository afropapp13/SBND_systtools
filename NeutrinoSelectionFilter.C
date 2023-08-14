#define NeutrinoSelectionFilter_cxx
#include "NeutrinoSelectionFilter.h"
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TH1D.h>
#include <TH2D.h>

#include <fstream>
#include <iomanip>
#include <vector>
#include <map>

#include "Tools.h"

#include "sbnanaobj/StandardRecord/SRGlobal.h"

using namespace std;

//----------------------------------------//

void NeutrinoSelectionFilter::Loop() {

	//----------------------------------------//

	TH1D::SetDefaultSumw2();
	TH2D::SetDefaultSumw2();

	//----------------------------------------//

	// Output Files

	TString filename = "sbnd_true_preselection.root";
	TFile* outfile = new TFile(filename,"recreate");
	cout << endl << "File " << filename << " to be created"<< endl << endl;

	//----------------------------------------//

	if (fChain == 0) return;
	Long64_t nentries = fChain->GetEntriesFast();
	Long64_t nbytes = 0, nb = 0;

	//----------------------------------------//

	// TTree Declarations

	TTree* truetree = new TTree("truepreselection","truepreselection");

	//----------------------------------------//

	// Common True & Reco branches	

	int ismc = -1;
	double pot = -1.;
	int run = -1;
	int subrun = -1;
	int event = -1;

	truetree->Branch("ismc",&ismc);
	truetree->Branch("pot",&pot);
	truetree->Branch("run",&run);
	truetree->Branch("subrun",&subrun);
	truetree->Branch("event",&event);

	//----------------------------------------//

	// True Branches

	int mc_nnu = -1;                                   truetree->Branch("mc_nnu",&mc_nnu);
	// number of dials for systematic variations
	//int mc_wgt_totarraysize = -1;                      truetree->Branch("mc_wgt_totarraysize",&mc_wgt_totarraysize);
	// total number of univeses
	//int mc_wgt_univ_totarraysize = -1;                 truetree->Branch("mc_wgt_univ_totarraysize",&mc_wgt_univ_totarraysize);

	vector<double> true_cc;                            truetree->Branch("true_cc",&true_cc);
	vector<double> true_nu_pdg;                        truetree->Branch("true_nu_pdg",&true_nu_pdg);
	vector<double> true_Ev;                            truetree->Branch("true_Ev",&true_Ev);
	vector<double> true_vx;                            truetree->Branch("true_vx",&true_vx);
	vector<double> true_vy;                            truetree->Branch("true_vy",&true_vy);
	vector<double> true_vz;                            truetree->Branch("true_vz",&true_vz);
	vector<double> true_pvx;                           truetree->Branch("true_pvx",&true_pvx);
	vector<double> true_pvy;                           truetree->Branch("true_pvy",&true_pvy);
	vector<double> true_pvz;                           truetree->Branch("true_pvz",&true_pvz);	
	vector<int>    true_mode;                          truetree->Branch("true_mode",&true_mode);
	vector<int>    true_hitnuc;                        truetree->Branch("true_hitnuc",&true_hitnuc);		

	// number of dials for each neutrino in a spill (117)
	//vector<int>    mc_wgt_length;                      truetree->Branch("mc_wgt_length",&mc_wgt_length);		
	// index where the weights for each neutrino start in a spill (0,117,234, ...)
	//vector<int>    mc_wgt_idx;                        truetree->Branch("mc_wgt_idx",&mc_wgt_idx);		

	// number of universes for a given dial
	//vector<int>    mc_wgt_univ_length;                      truetree->Branch("mc_wgt_univ_length",&mc_wgt_univ_length);		
	// index where a given dial universe starts from
	//vector<int>    mc_wgt_univ_idx;                        truetree->Branch("mc_wgt_univ_idx",&mc_wgt_univ_idx);		

	vector< vector<int> >    true_particle_pdg; 	   truetree->Branch("true_particle_pdg",&true_particle_pdg);
	vector< vector<double> > true_particle_E; 	   truetree->Branch("true_particle_E",&true_particle_E);
	vector< vector<double> > true_particle_p; 	   truetree->Branch("true_particle_p",&true_particle_p);
	vector< vector<double> > true_particle_px; 	   truetree->Branch("true_particle_px",&true_particle_px);
	vector< vector<double> > true_particle_py; 	   truetree->Branch("true_particle_py",&true_particle_py);
	vector< vector<double> > true_particle_pz; 	   truetree->Branch("true_particle_pz",&true_particle_pz);			
	vector< vector<double> > true_particle_phi;	   truetree->Branch("true_particle_phi",&true_particle_phi);
	vector< vector<double> > true_particle_costheta;   truetree->Branch("true_particle_costheta",&true_particle_costheta);
	vector< vector<double> > true_particle_startx;	   truetree->Branch("true_particle_startx",&true_particle_startx);
	vector< vector<double> > true_particle_starty;	   truetree->Branch("true_particle_starty",&true_particle_starty);
	vector< vector<double> > true_particle_startz;	   truetree->Branch("true_particle_startz",&true_particle_startz);

	// neutrinos in spill x dials x universes for each dial
	vector< std::map<TString,vector<double> > > syst_weights;     truetree->Branch("syst_weight",&syst_weights);

	Tools tools;

	//----------------------------------------//

	// Counters

	int totalcounter = 0;
	int true_ccincl = 0;	

	//----------------------------------------//

	// TTree for weights

	TChain* WeightChain = new TChain("globalTree","globalNeutrinoSelectionFilter");
	WeightChain->Add(fFile);

	if (WeightChain == 0) return;
	Long64_t Weightnentries = WeightChain->GetEntriesFast();
	Long64_t Weightnbytes = 0, Weightnb = 0;

	//caf::SRGlobal* srglobal;
	//WeightChain->SetBranchAddress( "global", &srglobal);

	caf::SRGlobal global;
	caf::SRGlobal* pglobal = &global;
	WeightChain->SetBranchAddress("global", &pglobal);
	outfile->cd();
	TTree globalOut("globalTree", "globalTree");
	globalOut.Branch("global", "caf::SRGlobal", &global);
	//assert(WeightChain->GetEntries() == 1);
	WeightChain->GetEntry(0);
	globalOut.Fill();
	globalOut.Write();

	// Loop over the WeightChain entries

	/*for (Long64_t Weightjentry = 0; Weightjentry < Weightnentries; Weightjentry++) {

	  Long64_t Weightientry = WeightChain->LoadTree(Weightjentry);
	  if (Weightientry < 0) break;
	  Weightnb = WeightChain->GetEntry(Weightjentry);   
	  Weightnbytes += Weightnb;
	  
	  std::vector<caf::SRWeightPSet> wgts = pglobal->wgts;
	  int wgts_size = wgts.size();

	  // Loop over the weight universes

	  for (int iw = 0; iw < wgts_size; iw++) {

	    std::string name = wgts[iw].name;
	    caf::ReweightType_t type = wgts[iw].type;
	    int nuniv = wgts[iw].nuniv;
	    std::vector<float> covmx = wgts[iw].covmx;
	    std::vector<caf::SRWeightMapEntry> map = wgts[iw].map;

	    cout << "iw = " << iw << ", name = " << name << ", type = " << type << ", nuniv = " << nuniv << endl;
                           
	    int nmap = map.size();
	    cout << "nmap = " << nmap << endl;
	    
	    // Loop over the map entries

	    for (auto const &pair: map) {

	      std::cout << "name = " << pair.param.name << endl;

	      int nknobs = pair.vals.size();

	      // Loop over the knobs

	      for (int iknob = 0; iknob < nknobs; iknob++) {

		cout << pair.vals[iknob] << " sigma, ";

	      } // end of the loop over the knobs
	      
	      cout << endl;

	    } // End of the loop of the map entries

	  } // End of the loop over the weight universes

	  cout << endl;
	                                                                                              
	  }*/ // End of the loop over the WeightChain entries, should be 1


	//----------------------------------------//

	// Loop over the events

	for (Long64_t jentry=0; jentry<nentries;jentry++) {

		//----------------------------------------//
      
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   
		nbytes += nb;

		//----------------------------------------//

		if (jentry%1000 == 0) cout << jentry/1000 << " k " << setprecision(3) << double(jentry)/double(fChain->GetEntries())*100. << " %"<< endl;

		//----------------------------------------//

		// Common branches

		ismc = rec_hdr_ismc;
		pot = rec_hdr_pot;
		run = rec_hdr_run;
		subrun = rec_hdr_subrun;
		event = rec_hdr_subevt;

		//mc_wgt_totarraysize = rec_mc_nu_wgt__totarraysize;
		//mc_wgt_univ_totarraysize = rec_mc_nu_wgt_univ__totarraysize;

		//----------------------------------------//

		// True Branches

		true_cc.clear();
		true_nu_pdg.clear();
		true_Ev.clear();
		true_vx.clear();
		true_vy.clear();
		true_vz.clear();
		true_pvx.clear();
		true_pvy.clear();
		true_pvz.clear();		
		true_mode.clear();
		true_hitnuc.clear();									

		//mc_wgt_length.clear();									
		//mc_wgt_idx.clear();									

		true_particle_pdg.clear();
		true_particle_E.clear();
		true_particle_p.clear();
		true_particle_px.clear();
		true_particle_py.clear();
		true_particle_pz.clear();
		true_particle_phi.clear();
		true_particle_costheta.clear();
		true_particle_startx.clear();
		true_particle_starty.clear();
		true_particle_startz.clear();

		// resize vectors of vectors containing the true primary particle information 

		true_particle_pdg.resize(rec_mc_nnu);
		true_particle_E.resize(rec_mc_nnu);
		true_particle_p.resize(rec_mc_nnu);
		true_particle_px.resize(rec_mc_nnu);	
		true_particle_py.resize(rec_mc_nnu);
		true_particle_pz.resize(rec_mc_nnu);					
		true_particle_phi.resize(rec_mc_nnu);
		true_particle_costheta.resize(rec_mc_nnu);
		true_particle_startx.resize(rec_mc_nnu);
		true_particle_starty.resize(rec_mc_nnu);
		true_particle_startz.resize(rec_mc_nnu);					

		//for (int i=0; i < 312;i++) {

		  //cout << "rec_mc_nu_wgt_univ__length[" << i << "] = " << rec_mc_nu_wgt_univ__length[i] << ", rec_mc_nu_wgt_univ__idx[" << i << "] = " << rec_mc_nu_wgt_univ__idx[i] << endl;

		//}


		//for (int i=0; i < 85368;i++) {

		//cout << "rec_mc_nu_wgt_univ[" << i << "] = " << rec_mc_nu_wgt_univ[i] << endl;

		//}

		//cout << endl;

		//----------------------------------------//

		// Loop over the true neutrino interactions in a given event		
		
		mc_nnu = rec_mc_nnu;
		//cout << "mc_nnu = " << mc_nnu << endl;

		//cout << "rec_mc_nu_wgt__totarraysize = " << rec_mc_nu_wgt__totarraysize << endl;
		//cout << "rec_mc_nu_wgt_univ__totarraysize = " << rec_mc_nu_wgt_univ__totarraysize << endl;

		for (int nnu = 0; nnu < (int)rec_mc_nnu; nnu++) {

		  //cout << "rec_mc_nu_wgt__length["<< nnu <<"] = " << rec_mc_nu_wgt__length[nnu] << endl;
		  //cout << "rec_mc_nu_wgt__idx["<< nnu <<"] = " << rec_mc_nu_wgt__idx[nnu] << endl;

			totalcounter++;			

			//----------------------------------------//	

			// True fully contained CC numu inclusive selection
			// Later, at an event selection stage,
			// we will demand only true CC νμ + n -> μ- + p + π0 primary interactions
			// Demand only numu CC events contained in the fiducial volume		

			/*			if (
				ismc && rec_mc_nu_iscc[nnu] == 1 && rec_mc_nu_pdg[nnu] == 14 && 
				tools.inFV(rec_mc_nu_position_x[nnu], rec_mc_nu_position_y[nnu], rec_mc_nu_position_z[nnu]) == 1
				) {*/

				//----------------------------------------//	

				true_ccincl++;			
				
				// Store the true neutrino info

				TVector3 true_vertex(rec_mc_nu_position_x[nnu],rec_mc_nu_position_y[nnu],rec_mc_nu_position_z[nnu]);

				true_nu_pdg.push_back(rec_mc_nu_pdg[nnu]);
				true_cc.push_back(rec_mc_nu_iscc[nnu]);
				true_Ev.push_back(rec_mc_nu_E[nnu]);
				true_vx.push_back(rec_mc_nu_position_x[nnu]);
				true_vy.push_back(rec_mc_nu_position_y[nnu]);
				true_vz.push_back(rec_mc_nu_position_z[nnu]);
				true_pvx.push_back(rec_mc_nu_momentum_x[nnu]);
				true_pvy.push_back(rec_mc_nu_momentum_y[nnu]);
				true_pvz.push_back(rec_mc_nu_momentum_z[nnu]);				
				true_mode.push_back(rec_mc_nu_genie_mode[nnu]);
				true_hitnuc.push_back(rec_mc_nu_hitnuc[nnu]);				

				//mc_wgt_length.push_back(rec_mc_nu_wgt__length[nnu]);				
				//mc_wgt_idx.push_back(rec_mc_nu_wgt__idx[nnu]);				

				//----------------------------------------//

				int NPrimaries = rec_mc_nu_nprim[nnu];
				//cout << "nnu = " << nnu << ", NPrimaries = " << NPrimaries << endl;

				// Loop over the primary MC particles
				// after the primary neutrino interaction

				//for (int iprim = 0; iprim < NPrimaries; iprim++) {
				for (int iprim = 0; iprim < 59; iprim++) {

					// Only fully contained true start points
				  /*
					if ( tools.inFV(rec_mc_nu_prim_gen_x[iprim], rec_mc_nu_prim_gen_y[iprim], rec_mc_nu_prim_gen_z[iprim]) ) {
				  */
				                TVector3 prim_vertex(rec_mc_nu_prim_gen_x[iprim],rec_mc_nu_prim_gen_y[iprim],rec_mc_nu_prim_gen_z[iprim]);
						TVector3 p(rec_mc_nu_prim_genp_x[iprim],rec_mc_nu_prim_genp_y[iprim],rec_mc_nu_prim_genp_z[iprim]);
						double distance_prim_vertex = (true_vertex-prim_vertex).Mag();
						//cout << "nnu =" << nnu << ", distance_prim_vertex = " << distance_prim_vertex << endl;
						
						// Feel in the vectors of vectors only when the true neutrino vertex 
						// and the primary particle vertex are the same
						if (distance_prim_vertex == 0) {
 
						  double phi = p.Phi(); // rad
						  double costheta = p.CosTheta(); // (-1,1)
						  double mag = p.Mag(); // GeV/c

						  true_particle_pdg.at(nnu).push_back(rec_mc_nu_prim_pdg[iprim]);
						  true_particle_E.at(nnu).push_back(rec_mc_nu_prim_genE[iprim]);
						  true_particle_p.at(nnu).push_back(mag);
						  true_particle_px.at(nnu).push_back(rec_mc_nu_prim_genp_x[iprim]);
						  true_particle_py.at(nnu).push_back(rec_mc_nu_prim_genp_y[iprim]);
						  true_particle_pz.at(nnu).push_back(rec_mc_nu_prim_genp_z[iprim]);																								
						  true_particle_phi.at(nnu).push_back(phi);
						  true_particle_costheta.at(nnu).push_back(costheta);	
						  true_particle_startx.at(nnu).push_back(rec_mc_nu_prim_gen_x[iprim]);
						  true_particle_starty.at(nnu).push_back(rec_mc_nu_prim_gen_y[iprim]);	
						  true_particle_startz.at(nnu).push_back(rec_mc_nu_prim_gen_z[iprim]);	

						} // end of primary vertex distance of 0 required

						/*
					} // End of only fully contained true start points
						*/																										
				

				} // End of the loop over the primary particles

				//cout << "nnu = " << nnu << ", true_particle_startz.at(nnu).size() = " << true_particle_startz.at(nnu).size() << endl;

				/*}*/ // End of the if-statement that the numu CC inclusive requirement is satified

		} // End of the loop over the true neutrino interactions in a given event

		//cout << "true_nu_pdg.size() = " << true_nu_pdg.size() << endl;
		//cout << "true_particle_startz.size() = " << true_particle_startz.size() << endl;

		//cout << endl;
		truetree->Fill();		

		//----------------------------------------//

	} // End of the loop over the events

	//----------------------------------------//

	cout << endl << "The file initially contains " << totalcounter << " true neutrino interactions" << endl;
	cout << true_ccincl << " are passing the true inclusive selection" << endl;	

	//----------------------------------------//

	outfile->cd();
	truetree->Write();
	outfile->Close();
	cout << endl << "File " << filename << " has been created"<< endl << endl;

	//----------------------------------------//

} // End of the program
