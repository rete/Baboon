//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jan  8 18:43:23 2013 by ROOT version 5.34/01
// from TTree tree/Shower variables
// found on file: single_calorimeterhit_pi-_80GeV_I0.root
//////////////////////////////////////////////////////////

#ifndef Analyze_h
#define Analyze_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <TH1.h> 
#include <TH2.h>

#include <iostream>
#include <string.h>
#include <vector>

// Fixed size dimensions of array or collections stored in the TTree if any.

typedef struct{
  std::vector<int> nhit;
  float radius;
  int end;
  int begin;
  float cog;
  int nlayer;
  int ebeam;
  float NhitPerLay;
  int coresize;
  int corebegin;
  int coreend;
  int trackMultiplicity;
  std::vector<int> trackLength;
  vector<int> density;
} Variable_t;

TH1 *hist0;
TH1 *hist1;
TH1 *hist2;
TH1 *hist3;
TH1 *hist4;
TH1 *hist5;
TH1 *hist6;
TH1 *hist7;
TH1 *hist8;
TH1 *hist9;
TH1 *hist10;
TH1 *hist11;

//fstream txtnhit;
//fstream txtnhit1;
//fstream txtnhit2;
//fstream txtnhit3;
//fstream txtbegin;
//fstream txtcog;
//fstream txtradius;
//fstream txttrackmulti;
//fstream txtcoresize;
//fstream txtcorelength;

class Analyze {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           Nhit;
   Int_t           Nhit1;
   Int_t           Nhit2;
   Int_t           Nhit3;
   Int_t           Nlayer;
   Int_t           Begin;
   Int_t           End;
   Float_t         Radius;
   Int_t           Ebeam;
   Int_t           Lasthit;
   vector<int>     *Density;
   Float_t         EHratio;
   Int_t           Hole;
   Int_t           Hole1;
   Float_t         CoG[3];
   Float_t         LongitudinalCut;
   Float_t         TransversalCut;
   Int_t           Ncores;
   vector<int>     *CoreSize;
   vector<int>     *CoreBegin;
   vector<int>     *CoreEnd;
   vector<double>  *CorePos;
   vector<double>  *CoreRadius;
   vector<int>     *TrackBegin;
   vector<int>     *TrackEnd;
   vector<int>     *TrackLength;
   Int_t           TrackMultiplicity;

   // List of branches
   TBranch        *b_Nhit;   //!
   TBranch        *b_Nhit1;   //!
   TBranch        *b_Nhit2;   //!
   TBranch        *b_Nhit3;   //!
   TBranch        *b_Nlayer;   //!
   TBranch        *b_Begin;   //!
   TBranch        *b_End;   //!
   TBranch        *b_Radius;   //!
   TBranch        *b_Ebeam;   //!
   TBranch        *b_Lasthit;   //!
   TBranch        *b_Density;   //!
   TBranch        *b_EHratio;   //!
   TBranch        *b_Hole;   //!
   TBranch        *b_Hole1;   //!
   TBranch        *b_LongitudinalCut;   //!
   TBranch        *b_TransversalCut;   //!
   TBranch        *b_CoG;   //!
   TBranch        *b_Ncores;   //!
   TBranch        *b_CoreSize;   //!
   TBranch        *b_CoreBegin;   //!
   TBranch        *b_CoreEnd;   //!
   TBranch        *b_CorePos;   //!
   TBranch        *b_CoreRadius;   //!
   TBranch        *b_TrackBegin;   //!
   TBranch        *b_TrackEnd;   //!
   TBranch        *b_TrackLength;   //!
   TBranch        *b_TrackMultiplicity;   //!

   Analyze(const char* input);
   virtual ~Analyze();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   std::vector<Variable_t> Loop(int E);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
 private:
   TFile *file;
};

#endif

#ifdef Analyze_cxx
Analyze::Analyze(const char* input) : fChain(0) 
{
  file = new TFile(input);  
  TTree *tree = (TTree*)file->Get("tree");

  Init(tree);
}

Analyze::~Analyze()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Analyze::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Analyze::LoadTree(Long64_t entry)
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

void Analyze::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   Density = 0;
   CoreSize = 0;
   CoreBegin = 0;
   CoreEnd = 0;
   CorePos = 0;
   CoreRadius = 0;
   TrackBegin = 0;
   TrackEnd = 0;
   TrackLength = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Nhit", &Nhit, &b_Nhit);
   fChain->SetBranchAddress("Nhit1", &Nhit1, &b_Nhit1);
   fChain->SetBranchAddress("Nhit2", &Nhit2, &b_Nhit2);
   fChain->SetBranchAddress("Nhit3", &Nhit3, &b_Nhit3);
   fChain->SetBranchAddress("Nlayer", &Nlayer, &b_Nlayer);
   fChain->SetBranchAddress("Begin", &Begin, &b_Begin);
   fChain->SetBranchAddress("End", &End, &b_End);
   fChain->SetBranchAddress("Radius", &Radius, &b_Radius);
   fChain->SetBranchAddress("Ebeam", &Ebeam, &b_Ebeam);
   fChain->SetBranchAddress("Lasthit", &Lasthit, &b_Lasthit);
   fChain->SetBranchAddress("Density", &Density, &b_Density);
   fChain->SetBranchAddress("EHratio", &EHratio, &b_EHratio);
   fChain->SetBranchAddress("Hole", &Hole, &b_Hole);
   fChain->SetBranchAddress("Hole1", &Hole1, &b_Hole1);
   fChain->SetBranchAddress("LongitudinalCut", &LongitudinalCut, &b_LongitudinalCut);
   fChain->SetBranchAddress("TransversalCut", &TransversalCut, &b_TransversalCut);
   fChain->SetBranchAddress("CoG", CoG, &b_CoG);
   fChain->SetBranchAddress("Ncores", &Ncores, &b_Ncores);
   fChain->SetBranchAddress("CoreSize", &CoreSize, &b_CoreSize);
   fChain->SetBranchAddress("CoreBegin", &CoreBegin, &b_CoreBegin);
   fChain->SetBranchAddress("CoreEnd", &CoreEnd, &b_CoreEnd);
   fChain->SetBranchAddress("CorePos", &CorePos, &b_CorePos);
   fChain->SetBranchAddress("CoreRadius", &CoreRadius, &b_CoreRadius);
   fChain->SetBranchAddress("TrackBegin", &TrackBegin, &b_TrackBegin);
   fChain->SetBranchAddress("TrackEnd", &TrackEnd, &b_TrackEnd);
   fChain->SetBranchAddress("TrackLength", &TrackLength, &b_TrackLength);
   fChain->SetBranchAddress("TrackMultiplicity", &TrackMultiplicity, &b_TrackMultiplicity);
   Notify();
}

Bool_t Analyze::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Analyze::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Analyze::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Analyze_cxx
