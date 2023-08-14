#ifndef CONSTANTS_H
#define CONSTANTS_H

#include "TString.h"
#include "TMath.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <map>

using namespace std;

namespace Constants {

	//----------------------------------------//

	// User  
  
	TString UserID = "apapadop";
	TString Experiment = "sbnd";
	
	// Constants

	static const double Units = 1E38; // so that the extracted cross-section is 10^{-38} cm^{2}
	static const double NA = 6.02214 * TMath::Power(10.,23.); // Avogadro's number, mol^-1 
	
	static const int NuMuPdg = 14, MuonPdg = 13, NuEPdg = 12, ElectronPdg = 11;
	static const int ProtonPdg = 2212, NeutronPdg = 2112;
	static const int AbsChargedPionPdg = 211, NeutralPionPdg = 111, KaonPdg = 321;
	static const int DeuteriumPdg = 1000010020, HeliumPdg = 1000020040, ArgonPdg = 1000180400;

	static const double MuonMass = 106, ProtonMass = 938.272, NeutronMass = 939.565; // MeV
	static const double MuonMass_GeV = 0.106, ProtonMass_GeV = 0.938272, NeutronMass_GeV = 0.939565; // GeV
	static const double DeltaM2 = TMath::Power(NeutronMass_GeV,2.) - TMath::Power(ProtonMass_GeV,2.); // GeV^2	
	
	// Argon 

	static const double A = 40.;
	static const double Z = 18.;
	static const double ArgonDensity = 1.3836; // g/cm^3
	static const double ArgonMolMass = 39.95; // g/mol
	static const double BE = 0.03; // Argon binding energy, GeV	

	// Detector

        double FVx = 199.15; //cm                                                                                                                                                                        
        double FVy = 200; //cm                                                                                                                                                                           
        double FVz = 500; //cm                                                                                                                                                                          

        double borderx = 10; //cm                                                                                                                                                                        
        double bordery = 10; //cm                                                                                                                                                                        
        double borderz = 10; //cm    

	//----------------------------------------//

	// Paths

	const TString PathToFiles          = "/"+Experiment+"/data/users/"+UserID+"/myEvents/OutputFiles/";

	//----------------------------------------//

	std::vector<TString> ccnc_label = {"All","CC","NC"};
	std::vector<TString> nu_label = {"All","numu","nue","numubar","nuebar"};
	std::vector<TString> inte_label = {"All","QE","MEC","RES","DIS","COH","Other"};

	//----------------------------------------//		

}
#endif
