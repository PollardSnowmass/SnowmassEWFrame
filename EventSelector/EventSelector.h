//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Feb 18 13:18:06 2013 by ROOT version 5.34/04
// from TTree SmEWObjNtup/SmEWObjNtup
// found on file: test.out.root
//////////////////////////////////////////////////////////

#ifndef EventSelector_h
#define EventSelector_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>

// Header file for the classes stored in the TTree if any.
#include <vector>

using namespace std;

// Fixed size dimensions of array or collections stored in the TTree if any.

class EventSelector : public TSelector {
    public :
        TTree          *fChain;   //!pointer to the analyzed TTree or TChain

        // Declaration of leaf types
        Int_t           jet_n;
        vector<float>   *jet_eta;
        vector<float>   *jet_phi;
        vector<float>   *jet_m;
        vector<float>   *jet_pt;
        vector<float>   *jet_E;
        Int_t           jet_vv_n;
        vector<float>   *jet_vv_eta;
        vector<float>   *jet_vv_phi;
        vector<float>   *jet_vv_m;
        vector<float>   *jet_vv_pt;
        vector<float>   *jet_vv_E;
        Int_t           jet_akt10_n;
        vector<float>   *jet_akt10_eta;
        vector<float>   *jet_akt10_phi;
        vector<float>   *jet_akt10_m;
        vector<float>   *jet_akt10_pt;
        vector<float>   *jet_akt10_E;
        Int_t           el_n;
        vector<int>     *el_charge;
        vector<float>   *el_eta;
        vector<float>   *el_phi;
        vector<float>   *el_m;
        vector<float>   *el_pt;
        vector<float>   *el_E;
        vector<int>     *el_trigger;
        Int_t           mu_n;
        vector<float>   *mu_charge;
        vector<float>   *mu_eta;
        vector<float>   *mu_phi;
        vector<float>   *mu_m;
        vector<float>   *mu_pt;
        vector<float>   *mu_E;
        vector<int>     *mu_trigger;
        Float_t         met;
        Float_t         met_phi;
        Float_t         sumet;
        Float_t         truth_mass;

        float ht;
        float ht2j;
        float ht_met;
        float ht2j_met;
        float ht_lep;
        float ht_lep_met;

        float m_ll;
        float vv_max_dijet_mass;
        float vv_leading_dijet_mass;
        float vv_leading_dijet_deta;
        float vv_4body_mass;

        bool pass_cuts;

        void Reset();
        void FillMll();
        void FillVVMasses();
        void FillHT();

        TFile *fout;
        TTree *tout;


        // List of branches
        TBranch        *b_jet_n;   //!
        TBranch        *b_jet_eta;   //!
        TBranch        *b_jet_phi;   //!
        TBranch        *b_jet_m;   //!
        TBranch        *b_jet_pt;   //!
        TBranch        *b_jet_E;   //!
        TBranch        *b_jet_vv_n;   //!
        TBranch        *b_jet_vv_eta;   //!
        TBranch        *b_jet_vv_phi;   //!
        TBranch        *b_jet_vv_m;   //!
        TBranch        *b_jet_vv_pt;   //!
        TBranch        *b_jet_vv_E;   //!
        TBranch        *b_jet_akt10_n;   //!
        TBranch        *b_jet_akt10_eta;   //!
        TBranch        *b_jet_akt10_phi;   //!
        TBranch        *b_jet_akt10_m;   //!
        TBranch        *b_jet_akt10_pt;   //!
        TBranch        *b_jet_akt10_E;   //!
        TBranch        *b_el_n;   //!
        TBranch        *b_el_charge;   //!
        TBranch        *b_el_eta;   //!
        TBranch        *b_el_phi;   //!
        TBranch        *b_el_m;   //!
        TBranch        *b_el_pt;   //!
        TBranch        *b_el_E;   //!
        TBranch        *b_el_trigger;   //!
        TBranch        *b_mu_n;   //!
        TBranch        *b_mu_charge;   //!
        TBranch        *b_mu_eta;   //!
        TBranch        *b_mu_phi;   //!
        TBranch        *b_mu_m;   //!
        TBranch        *b_mu_pt;   //!
        TBranch        *b_mu_E;   //!
        TBranch        *b_mu_trigger;   //!
        TBranch        *b_met;   //!
        TBranch        *b_met_phi;   //!
        TBranch        *b_sumet;   //!
        TBranch        *b_truth_mass;   //!

        EventSelector(TTree * /*tree*/ =0) : fChain(0) { }
        virtual ~EventSelector() { }
        virtual Int_t   Version() const { return 2; }
        virtual void    Begin(TTree *tree);
        virtual void    SlaveBegin(TTree *tree);
        virtual void    Init(TTree *tree);
        virtual Bool_t  Notify();
        virtual Bool_t  Process(Long64_t entry);
        virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
        virtual void    SetOption(const char *option) { fOption = option; }
        virtual void    SetObject(TObject *obj) { fObject = obj; }
        virtual void    SetInputList(TList *input) { fInput = input; }
        virtual TList  *GetOutputList() const { return fOutput; }
        virtual void    SlaveTerminate();
        virtual void    Terminate();

        ClassDef(EventSelector,0);
};

#endif

#ifdef EventSelector_cxx
void EventSelector::Init(TTree *tree)
{
    // The Init() function is called when the selector needs to initialize
    // a new tree or chain. Typically here the branch addresses and branch
    // pointers of the tree will be set.
    // It is normally not necessary to make changes to the generated
    // code, but the routine can be extended by the user if needed.
    // Init() will be called many times when running on PROOF
    // (once per file to be processed).

    // Set object pointer
    jet_eta = 0;
    jet_phi = 0;
    jet_m = 0;
    jet_pt = 0;
    jet_E = 0;
    jet_vv_eta = 0;
    jet_vv_phi = 0;
    jet_vv_m = 0;
    jet_vv_pt = 0;
    jet_vv_E = 0;
    jet_akt10_eta = 0;
    jet_akt10_phi = 0;
    jet_akt10_m = 0;
    jet_akt10_pt = 0;
    jet_akt10_E = 0;
    el_charge = 0;
    el_eta = 0;
    el_phi = 0;
    el_m = 0;
    el_pt = 0;
    el_E = 0;
    el_trigger = 0;
    mu_charge = 0;
    mu_eta = 0;
    mu_phi = 0;
    mu_m = 0;
    mu_pt = 0;
    mu_E = 0;
    mu_trigger = 0;
    // Set branch addresses and branch pointers
    if (!tree) return;
    fChain = tree;
    fChain->SetMakeClass(1);

    fChain->SetBranchAddress("jet_n", &jet_n, &b_jet_n);
    fChain->SetBranchAddress("jet_eta", &jet_eta, &b_jet_eta);
    fChain->SetBranchAddress("jet_phi", &jet_phi, &b_jet_phi);
    fChain->SetBranchAddress("jet_m", &jet_m, &b_jet_m);
    fChain->SetBranchAddress("jet_pt", &jet_pt, &b_jet_pt);
    fChain->SetBranchAddress("jet_E", &jet_E, &b_jet_E);
    fChain->SetBranchAddress("jet_vv_n", &jet_vv_n, &b_jet_vv_n);
    fChain->SetBranchAddress("jet_vv_eta", &jet_vv_eta, &b_jet_vv_eta);
    fChain->SetBranchAddress("jet_vv_phi", &jet_vv_phi, &b_jet_vv_phi);
    fChain->SetBranchAddress("jet_vv_m", &jet_vv_m, &b_jet_vv_m);
    fChain->SetBranchAddress("jet_vv_pt", &jet_vv_pt, &b_jet_vv_pt);
    fChain->SetBranchAddress("jet_vv_E", &jet_vv_E, &b_jet_vv_E);
    fChain->SetBranchAddress("jet_akt10_n", &jet_akt10_n, &b_jet_akt10_n);
    fChain->SetBranchAddress("jet_akt10_eta", &jet_akt10_eta, &b_jet_akt10_eta);
    fChain->SetBranchAddress("jet_akt10_phi", &jet_akt10_phi, &b_jet_akt10_phi);
    fChain->SetBranchAddress("jet_akt10_m", &jet_akt10_m, &b_jet_akt10_m);
    fChain->SetBranchAddress("jet_akt10_pt", &jet_akt10_pt, &b_jet_akt10_pt);
    fChain->SetBranchAddress("jet_akt10_E", &jet_akt10_E, &b_jet_akt10_E);
    fChain->SetBranchAddress("el_n", &el_n, &b_el_n);
    fChain->SetBranchAddress("el_charge", &el_charge, &b_el_charge);
    fChain->SetBranchAddress("el_eta", &el_eta, &b_el_eta);
    fChain->SetBranchAddress("el_phi", &el_phi, &b_el_phi);
    fChain->SetBranchAddress("el_m", &el_m, &b_el_m);
    fChain->SetBranchAddress("el_pt", &el_pt, &b_el_pt);
    fChain->SetBranchAddress("el_E", &el_E, &b_el_E);
    fChain->SetBranchAddress("el_trigger", &el_trigger, &b_el_trigger);
    fChain->SetBranchAddress("mu_n", &mu_n, &b_mu_n);
    fChain->SetBranchAddress("mu_charge", &mu_charge, &b_mu_charge);
    fChain->SetBranchAddress("mu_eta", &mu_eta, &b_mu_eta);
    fChain->SetBranchAddress("mu_phi", &mu_phi, &b_mu_phi);
    fChain->SetBranchAddress("mu_m", &mu_m, &b_mu_m);
    fChain->SetBranchAddress("mu_pt", &mu_pt, &b_mu_pt);
    fChain->SetBranchAddress("mu_E", &mu_E, &b_mu_E);
    fChain->SetBranchAddress("mu_trigger", &mu_trigger, &b_mu_trigger);
    fChain->SetBranchAddress("met", &met, &b_met);
    fChain->SetBranchAddress("met_phi", &met_phi, &b_met_phi);
    fChain->SetBranchAddress("sumet", &sumet, &b_sumet);
    fChain->SetBranchAddress("truth_mass", &truth_mass, &b_truth_mass);
}

Bool_t EventSelector::Notify()
{
    // The Notify() function is called when a new file is opened. This
    // can be either for a new TTree in a TChain or when when a new TTree
    // is started when using PROOF. It is normally not necessary to make changes
    // to the generated code, but the routine can be extended by the
    // user if needed. The return value is currently not used.

    return kTRUE;
}

#endif // #ifdef EventSelector_cxx
