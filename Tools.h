#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <vector>

#include "TString.h"
#include "TMath.h"
#include <TVector3.h>
#include "TH1D.h"
#include "TH2D.h"
#include <TVectorD.h>
#include <TMatrixD.h>

#include "Constants.h"

using namespace Constants;

class Tools {

//	private:

	public:

		// Default constructor
		Tools(){ fFVx = FVx; fFVy = FVy; fFVz = FVz; fborderx = borderx; fbordery = bordery; fborderz = borderz;}

		// Default destructor
		~Tools(){}

		bool    inFV(double x, double y, double z);
		bool    inFVVector(TVector3 vector);	

		double  fFVx; // cm
		double  fFVy; // cm
		double  fFVz; // cm
		double  fborderx; // cm
		double  fbordery; // cm
		double  fborderz; // cm

};


//----------------------------------------//

bool Tools::inFVVector(TVector3 vector) {

	bool contained = true;

	if( !(TMath::Abs(vector.X()) < (fFVx - fborderx) && ( TMath::Abs(vector.X()) > fborderx)) ) { contained = false; }
	if( !(vector.Y() < (fFVy - fbordery) && ( vector.Y() > fbordery)) ) { contained = false; }
	if( !(vector.Z() < (fFVz - fborderz) && ( vector.Z() > fborderz)) ) { contained = false; }			 

	return contained;
};

//----------------------------------------//

bool Tools::inFV(double x, double y, double z) {

	bool contained = true;

	if( !(TMath::Abs(x) < (fFVx - fborderx) && ( TMath::Abs(x) > fborderx)) ) { contained = false; }
	if( !(y < (fFVy - fbordery) && ( y > fbordery)) ) { contained = false; }
	if( !(z < (fFVz - fborderz) && ( z > fborderz)) ) { contained = false; }		 

	return contained;
};
