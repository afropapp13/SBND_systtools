//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Aug  8 15:10:44 2023 by ROOT version 6.26/07
// from TTree truepreselection/truepreselection
// found on file: sbnd_true_preselection_v09_56_00.root
//////////////////////////////////////////////////////////

#ifndef truepreselection_h
#define truepreselection_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"

#include "Constants.h"

class truepreselection {

 private:
  TString fdial;
  int funiverse;


public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           ismc;
   Double_t        pot;
   Int_t           run;
   Int_t           subrun;
   Int_t           event;
   Int_t           mc_nnu;
   vector<double>  *true_cc;
   vector<double>  *true_nu_pdg;
   vector<double>  *true_Ev;
   vector<double>  *true_vx;
   vector<double>  *true_vy;
   vector<double>  *true_vz;
   vector<double>  *true_pvx;
   vector<double>  *true_pvy;
   vector<double>  *true_pvz;
   vector<int>     *true_mode;
   vector<int>     *true_hitnuc;
   vector<vector<int> > *true_particle_pdg;
   vector<vector<double> > *true_particle_E;
   vector<vector<double> > *true_particle_p;
   vector<vector<double> > *true_particle_px;
   vector<vector<double> > *true_particle_py;
   vector<vector<double> > *true_particle_pz;
   vector<vector<double> > *true_particle_phi;
   vector<vector<double> > *true_particle_costheta;
   vector<vector<double> > *true_particle_startx;
   vector<vector<double> > *true_particle_starty;
   vector<vector<double> > *true_particle_startz;

   // List of branches
   TBranch        *b_ismc;   //!
   TBranch        *b_pot;   //!
   TBranch        *b_run;   //!
   TBranch        *b_subrun;   //!
   TBranch        *b_event;   //!
   TBranch        *b_mc_nnu;   //!
   TBranch        *b_true_cc;   //!
   TBranch        *b_true_nu_pdg;   //!
   TBranch        *b_true_Ev;   //!
   TBranch        *b_true_vx;   //!
   TBranch        *b_true_vy;   //!
   TBranch        *b_true_vz;   //!
   TBranch        *b_true_pvx;   //!
   TBranch        *b_true_pvy;   //!
   TBranch        *b_true_pvz;   //!
   TBranch        *b_true_mode;   //!
   TBranch        *b_true_hitnuc;   //!
   TBranch        *b_true_particle_pdg;   //!
   TBranch        *b_true_particle_E;   //!
   TBranch        *b_true_particle_p;   //!
   TBranch        *b_true_particle_px;   //!
   TBranch        *b_true_particle_py;   //!
   TBranch        *b_true_particle_pz;   //!
   TBranch        *b_true_particle_phi;   //!
   TBranch        *b_true_particle_costheta;   //!
   TBranch        *b_true_particle_startx;   //!
   TBranch        *b_true_particle_starty;   //!
   TBranch        *b_true_particle_startz;   //!

   truepreselection(TString dial = "CV", int universe = -1,TTree *tree=0);
   virtual ~truepreselection();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef truepreselection_cxx
truepreselection::truepreselection(TString dial, int universe, TTree *tree) : fChain(0) 
{

    fdial = dial; 
    funiverse = universe;

// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("sbnd_true_preselection.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("sbnd_true_preselection.root");
      }
      f->GetObject("truepreselection",tree);

   }
   Init(tree);
}

truepreselection::~truepreselection()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t truepreselection::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t truepreselection::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void truepreselection::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   true_cc = 0;
   true_nu_pdg = 0;
   true_Ev = 0;
   true_vx = 0;
   true_vy = 0;
   true_vz = 0;
   true_pvx = 0;
   true_pvy = 0;
   true_pvz = 0;
   true_mode = 0;
   true_hitnuc = 0;
   true_particle_pdg = 0;
   true_particle_E = 0;
   true_particle_p = 0;
   true_particle_px = 0;
   true_particle_py = 0;
   true_particle_pz = 0;
   true_particle_phi = 0;
   true_particle_costheta = 0;
   true_particle_startx = 0;
   true_particle_starty = 0;
   true_particle_startz = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("ismc", &ismc, &b_ismc);
   fChain->SetBranchAddress("pot", &pot, &b_pot);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("subrun", &subrun, &b_subrun);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("mc_nnu", &mc_nnu, &b_mc_nnu);
   fChain->SetBranchAddress("true_cc", &true_cc, &b_true_cc);
   fChain->SetBranchAddress("true_nu_pdg", &true_nu_pdg, &b_true_nu_pdg);
   fChain->SetBranchAddress("true_Ev", &true_Ev, &b_true_Ev);
   fChain->SetBranchAddress("true_vx", &true_vx, &b_true_vx);
   fChain->SetBranchAddress("true_vy", &true_vy, &b_true_vy);
   fChain->SetBranchAddress("true_vz", &true_vz, &b_true_vz);
   fChain->SetBranchAddress("true_pvx", &true_pvx, &b_true_pvx);
   fChain->SetBranchAddress("true_pvy", &true_pvy, &b_true_pvy);
   fChain->SetBranchAddress("true_pvz", &true_pvz, &b_true_pvz);
   fChain->SetBranchAddress("true_mode", &true_mode, &b_true_mode);
   fChain->SetBranchAddress("true_hitnuc", &true_hitnuc, &b_true_hitnuc);
   fChain->SetBranchAddress("true_particle_pdg", &true_particle_pdg, &b_true_particle_pdg);
   fChain->SetBranchAddress("true_particle_E", &true_particle_E, &b_true_particle_E);
   fChain->SetBranchAddress("true_particle_p", &true_particle_p, &b_true_particle_p);
   fChain->SetBranchAddress("true_particle_px", &true_particle_px, &b_true_particle_px);
   fChain->SetBranchAddress("true_particle_py", &true_particle_py, &b_true_particle_py);
   fChain->SetBranchAddress("true_particle_pz", &true_particle_pz, &b_true_particle_pz);
   fChain->SetBranchAddress("true_particle_phi", &true_particle_phi, &b_true_particle_phi);
   fChain->SetBranchAddress("true_particle_costheta", &true_particle_costheta, &b_true_particle_costheta);
   fChain->SetBranchAddress("true_particle_startx", &true_particle_startx, &b_true_particle_startx);
   fChain->SetBranchAddress("true_particle_starty", &true_particle_starty, &b_true_particle_starty);
   fChain->SetBranchAddress("true_particle_startz", &true_particle_startz, &b_true_particle_startz);
   Notify();
}

Bool_t truepreselection::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void truepreselection::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t truepreselection::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef truepreselection_cxx
