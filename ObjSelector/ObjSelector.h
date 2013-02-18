//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Apr  3 10:29:07 2012 by ROOT version 5.32/01
// from TTree truth_physics/truth_physics
// found on file: ../ObjPythia/out/kkg5000.root
//////////////////////////////////////////////////////////

#ifndef ObjSelector_h
#define ObjSelector_h

#include <TROOT.h>
#include <TChain.h>
#include <TTree.h>
#include <TFile.h>
#include <TSelector.h>
#include "TRandom3.h"
#include <TH1F.h>
#include <TH2F.h>

// Header file for the classes stored in the TTree if any.
#include <vector>

using namespace std;

class ObjSelector : public TSelector {
    public :
        TTree          *fChain;   //!pointer to the analyzed TTree or TChain

        // Declaration of leaf types
        Int_t           mc_n;
        vector<float>   *mc_E;
        vector<float>   *mc_pt;
        vector<float>   *mc_m;
        vector<float>   *mc_eta;
        vector<float>   *mc_phi;
        vector<int>     *mc_status;
        vector<int>     *mc_barcode;
        vector<vector<int> > *mc_parents;
        vector<vector<int> > *mc_children;
        vector<int>     *mc_pdgId;
        vector<float>   *mc_charge;
        vector<vector<int> > *mc_child_index;
        vector<vector<int> > *mc_parent_index;
        Int_t           mcevt_n;
        vector<int>     *mcevt_signal_process_id;
        vector<int>     *mcevt_event_number;
        vector<double>  *mcevt_event_scale;
        vector<double>  *mcevt_alphaQCD;
        vector<double>  *mcevt_alphaQED;
        vector<int>     *mcevt_pdf_id1;
        vector<int>     *mcevt_pdf_id2;
        vector<double>  *mcevt_pdf_x1;
        vector<double>  *mcevt_pdf_x2;
        vector<double>  *mcevt_pdf_scale;
        vector<double>  *mcevt_pdf1;
        vector<double>  *mcevt_pdf2;
        vector<vector<double> > *mcevt_weight;
        Int_t           jet_n;
        vector<float>   *jet_eta;
        vector<float>   *jet_phi;
        vector<float>   *jet_E;
        vector<float>   *jet_pt;
        vector<float>   *jet_m;
        vector<int>     *jet_bjet;
        Int_t           jet_akt10_n;
        vector<float>   *jet_akt10_eta;
        vector<float>   *jet_akt10_phi;
        vector<float>   *jet_akt10_E;
        vector<float>   *jet_akt10_pt;
        vector<float>   *jet_akt10_m;
        vector<int>     *el_indices;
        vector<int>     *mu_indices;

        // analysis variables.
        int mu;
        bool do_reco;
        bool dobtag_veto;

        float el_pt_cut;
        float el_eta_cut;
        float mu_pt_cut;
        float mu_eta_cut;
        bool mu_eff_corr; 
        float jet_pt_cut;
        float jet_eta_cut;
        float jet_vv_pt_cut;
        float jet_vv_eta_cut;
        float jet_akt10_pt_cut;
        float jet_akt10_eta_cut;
        float jet_akt10_m_cut;
        float jet_bjet_pt_thresh;
        bool has_bjet;

        int jet_n_out;
        vector<float> jet_eta_out;
        vector<float> jet_phi_out;
        vector<float> jet_m_out;
        vector<float> jet_pt_out;
        vector<float> jet_E_out;

        int jet_vv_n_out;
        vector<float> jet_vv_eta_out;
        vector<float> jet_vv_phi_out;
        vector<float> jet_vv_m_out;
        vector<float> jet_vv_pt_out;
        vector<float> jet_vv_E_out;

        int jet_akt10_n_out;
        vector<float> jet_akt10_eta_out;
        vector<float> jet_akt10_phi_out;
        vector<float> jet_akt10_m_out;
        vector<float> jet_akt10_pt_out;
        vector<float> jet_akt10_E_out;

        int el_n_out;
        vector<int> el_charge_out;
        vector<float> el_eta_out;
        vector<float> el_phi_out;
        vector<float> el_m_out;
        vector<float> el_pt_out;
        vector<float> el_E_out;
        vector<int> el_trigger_out;

        int mu_n_out;
        vector<float> mu_charge_out;
        vector<float> mu_eta_out;
        vector<float> mu_phi_out;
        vector<float> mu_m_out;
        vector<float> mu_pt_out;
        vector<float> mu_E_out;
        vector<int> mu_trigger_out;

        float nu_pt_out;
        float nu_eta_out;
        float nu_phi_out;
        float nu_E_out;

        float sumet_out;

        float m_ll_out;
        float m_ll01_out;
        float m_ll02_out;
        float m_ll03_out;
        float m_ll12_out;
        float m_ll13_out;
        float m_ll23_out;
        float m_llll_out;

        float met_out;
        float met_phi_out;

        float ht_out;
        float ht_met_out;
        float ht2j_out;
        float ht2j_met_out;
        float ht_lep_out;
        float ht_lep_met_out;

        int ljets_lepjet_idx_out;
        float ljets_leptop_mass_out;
        int ljets_hadjet_idx_out;
        float ljets_hadtop_mass_out;
        float ljets_wmass_out;
        float ljets_mass_out;

        float vv_max_dijet_mass_out;
        float vv_leading_dijet_mass_out;
        float vv_leading_dijet_deta_out;
        float vv_4body_mass_out;

        float wz_m_wz_out;
        float wz_m_w_out;
        float wz_m_z_out;

        float truth_mass;

        TFile *fout;
        TTree *tout;


        // List of branches
        TBranch        *b_mc_n;   //!
        TBranch        *b_mc_E;   //!
        TBranch        *b_mc_pt;   //!
        TBranch        *b_mc_m;   //!
        TBranch        *b_mc_eta;   //!
        TBranch        *b_mc_phi;   //!
        TBranch        *b_mc_status;   //!
        TBranch        *b_mc_barcode;   //!
        TBranch        *b_mc_parents;   //!
        TBranch        *b_mc_children;   //!
        TBranch        *b_mc_pdgId;   //!
        TBranch        *b_mc_charge;   //!
        TBranch        *b_mc_child_index;   //!
        TBranch        *b_mc_parent_index;   //!
        TBranch        *b_mcevt_n;   //!
        TBranch        *b_mcevt_signal_process_id;   //!
        TBranch        *b_mcevt_event_number;   //!
        TBranch        *b_mcevt_event_scale;   //!
        TBranch        *b_mcevt_alphaQCD;   //!
        TBranch        *b_mcevt_alphaQED;   //!
        TBranch        *b_mcevt_pdf_id1;   //!
        TBranch        *b_mcevt_pdf_id2;   //!
        TBranch        *b_mcevt_pdf_x1;   //!
        TBranch        *b_mcevt_pdf_x2;   //!
        TBranch        *b_mcevt_pdf_scale;   //!
        TBranch        *b_mcevt_pdf1;   //!
        TBranch        *b_mcevt_pdf2;   //!
        TBranch        *b_mcevt_weight;   //!
        TBranch        *b_jet_n;   //!
        TBranch        *b_jet_eta;   //!
        TBranch        *b_jet_phi;   //!
        TBranch        *b_jet_E;   //!
        TBranch        *b_jet_pt;   //!
        TBranch        *b_jet_m;   //!
        TBranch        *b_jet_bjet;   //!
        TBranch        *b_jet_akt10_n;   //!
        TBranch        *b_jet_akt10_eta;   //!
        TBranch        *b_jet_akt10_phi;   //!
        TBranch        *b_jet_akt10_E;   //!
        TBranch        *b_jet_akt10_pt;   //!
        TBranch        *b_jet_akt10_m;   //!
        TBranch        *b_el_indices;   //!
        TBranch        *b_mu_indices;   //!

        ObjSelector(TTree * /*tree*/ =0) : fChain(0) { }
        virtual ~ObjSelector() { }
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

        // added by chris.
        TRandom3 *myRand;
        void Reset();
        void ReadConfig();
        float Smear(TRandom3 *rand1, float width);
        bool PassEfficiency(TRandom3 *rand1, float efficiency);
        float quadsum(const float &x, const float &y);

        void FillGoodElectrons();
        void FillGoodMuons();
        void FillGoodJets();
        void FillGoodVVJets();
        void FillGoodAkt10Jets();

        bool IsGoodElectron(int iEl);
        bool IsGoodMuon(int iMu);
        bool IsGoodJet(int iJet);
        bool IsGoodVVJet(int iJet);
        bool IsGoodAkt10Jet(int iJet);

        void ElJetOverlap();
        void JetElOverlap();

        void FillTrigger();
        void FillSumEt();
        void FillMET();
        void FillTruthMass();

        void SmearElectrons();
        void SmearMuons();
        void SmearJets();
        void SmearAkt10Jets();

        ClassDef(ObjSelector,0);
};

#endif

#ifdef ObjSelector_cxx
void ObjSelector::Init(TTree *tree)
{
    // The Init() function is called when the selector needs to initialize
    // a new tree or chain. Typically here the branch addresses and branch
    // pointers of the tree will be set.
    // It is normally not necessary to make changes to the generated
    // code, but the routine can be extended by the user if needed.
    // Init() will be called many times when running on PROOF
    // (once per file to be processed).

    // Set object pointer
    mc_E = 0;
    mc_pt = 0;
    mc_m = 0;
    mc_eta = 0;
    mc_phi = 0;
    mc_status = 0;
    mc_barcode = 0;
    mc_parents = 0;
    mc_children = 0;
    mc_pdgId = 0;
    mc_charge = 0;
    mc_child_index = 0;
    mc_parent_index = 0;
    mcevt_signal_process_id = 0;
    mcevt_event_number = 0;
    mcevt_event_scale = 0;
    mcevt_alphaQCD = 0;
    mcevt_alphaQED = 0;
    mcevt_pdf_id1 = 0;
    mcevt_pdf_id2 = 0;
    mcevt_pdf_x1 = 0;
    mcevt_pdf_x2 = 0;
    mcevt_pdf_scale = 0;
    mcevt_pdf1 = 0;
    mcevt_pdf2 = 0;
    mcevt_weight = 0;
    jet_eta = 0;
    jet_phi = 0;
    jet_E = 0;
    jet_pt = 0;
    jet_m = 0;
    jet_bjet = 0;
    jet_akt10_eta = 0;
    jet_akt10_phi = 0;
    jet_akt10_E = 0;
    jet_akt10_pt = 0;
    jet_akt10_m = 0;
    el_indices = 0;
    mu_indices = 0;

    // Set branch addresses and branch pointers
    if (!tree) return;
    fChain = tree;
    fChain->SetMakeClass(1);

    fChain->SetBranchAddress("mc_n", &mc_n, &b_mc_n);
    fChain->SetBranchAddress("mc_E", &mc_E, &b_mc_E);
    fChain->SetBranchAddress("mc_pt", &mc_pt, &b_mc_pt);
    fChain->SetBranchAddress("mc_m", &mc_m, &b_mc_m);
    fChain->SetBranchAddress("mc_eta", &mc_eta, &b_mc_eta);
    fChain->SetBranchAddress("mc_phi", &mc_phi, &b_mc_phi);
    fChain->SetBranchAddress("mc_status", &mc_status, &b_mc_status);
    fChain->SetBranchAddress("mc_barcode", &mc_barcode, &b_mc_barcode);
    fChain->SetBranchAddress("mc_parents", &mc_parents, &b_mc_parents);
    fChain->SetBranchAddress("mc_children", &mc_children, &b_mc_children);
    fChain->SetBranchAddress("mc_pdgId", &mc_pdgId, &b_mc_pdgId);
    fChain->SetBranchAddress("mc_charge", &mc_charge, &b_mc_charge);
    fChain->SetBranchAddress("mc_child_index", &mc_child_index, &b_mc_child_index);
    fChain->SetBranchAddress("mc_parent_index", &mc_parent_index, &b_mc_parent_index);
    fChain->SetBranchAddress("mcevt_n", &mcevt_n, &b_mcevt_n);
    fChain->SetBranchAddress("mcevt_signal_process_id", &mcevt_signal_process_id, &b_mcevt_signal_process_id);
    fChain->SetBranchAddress("mcevt_event_number", &mcevt_event_number, &b_mcevt_event_number);
    fChain->SetBranchAddress("mcevt_event_scale", &mcevt_event_scale, &b_mcevt_event_scale);
    fChain->SetBranchAddress("mcevt_alphaQCD", &mcevt_alphaQCD, &b_mcevt_alphaQCD);
    fChain->SetBranchAddress("mcevt_alphaQED", &mcevt_alphaQED, &b_mcevt_alphaQED);
    fChain->SetBranchAddress("mcevt_pdf_id1", &mcevt_pdf_id1, &b_mcevt_pdf_id1);
    fChain->SetBranchAddress("mcevt_pdf_id2", &mcevt_pdf_id2, &b_mcevt_pdf_id2);
    fChain->SetBranchAddress("mcevt_pdf_x1", &mcevt_pdf_x1, &b_mcevt_pdf_x1);
    fChain->SetBranchAddress("mcevt_pdf_x2", &mcevt_pdf_x2, &b_mcevt_pdf_x2);
    fChain->SetBranchAddress("mcevt_pdf_scale", &mcevt_pdf_scale, &b_mcevt_pdf_scale);
    fChain->SetBranchAddress("mcevt_pdf1", &mcevt_pdf1, &b_mcevt_pdf1);
    fChain->SetBranchAddress("mcevt_pdf2", &mcevt_pdf2, &b_mcevt_pdf2);
    fChain->SetBranchAddress("mcevt_weight", &mcevt_weight, &b_mcevt_weight);
    fChain->SetBranchAddress("jet_n", &jet_n, &b_jet_n);
    fChain->SetBranchAddress("jet_eta", &jet_eta, &b_jet_eta);
    fChain->SetBranchAddress("jet_phi", &jet_phi, &b_jet_phi);
    fChain->SetBranchAddress("jet_E", &jet_E, &b_jet_E);
    fChain->SetBranchAddress("jet_pt", &jet_pt, &b_jet_pt);
    fChain->SetBranchAddress("jet_m", &jet_m, &b_jet_m);
    fChain->SetBranchAddress("jet_bjet", &jet_bjet, &b_jet_bjet);
    fChain->SetBranchAddress("jet_akt10_n", &jet_akt10_n, &b_jet_akt10_n);
    fChain->SetBranchAddress("jet_akt10_eta", &jet_akt10_eta, &b_jet_akt10_eta);
    fChain->SetBranchAddress("jet_akt10_phi", &jet_akt10_phi, &b_jet_akt10_phi);
    fChain->SetBranchAddress("jet_akt10_E", &jet_akt10_E, &b_jet_akt10_E);
    fChain->SetBranchAddress("jet_akt10_pt", &jet_akt10_pt, &b_jet_akt10_pt);
    fChain->SetBranchAddress("jet_akt10_m", &jet_akt10_m, &b_jet_akt10_m);
    fChain->SetBranchAddress("el_indices", &el_indices, &b_el_indices);
    fChain->SetBranchAddress("mu_indices", &mu_indices, &b_mu_indices);
}

Bool_t ObjSelector::Notify()
{
    // The Notify() function is called when a new file is opened. This
    // can be either for a new TTree in a TChain or when when a new TTree
    // is started when using PROOF. It is normally not necessary to make changes
    // to the generated code, but the routine can be extended by the
    // user if needed. The return value is currently not used.

    return kTRUE;
}

#ifdef __CINT__
#pragma link C++ class ObjSelector+;
#endif

#endif // #ifdef ObjSelector_cxx
