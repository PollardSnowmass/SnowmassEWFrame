#define EventSelector_cxx

#include "EventSelector.h"
#include "LoadVector.h"
#include <iostream>
#include <TLorentzVector.h>

void EventSelector::Begin(TTree * /*tree*/) {
    // The Begin() function is called at the start of the query.
    // When running with PROOF Begin() is only called on the client.
    // The tree argument is deprecated (on PROOF 0 is passed).

    return;
}

void EventSelector::SlaveBegin(TTree * /*tree*/) {
    // The SlaveBegin() function is called after the Begin() function.
    // When running with PROOF SlaveBegin() is called on each slave server.
    // The tree argument is deprecated (on PROOF 0 is passed).

    TString fname = GetOption();
    fout = new TFile(fname, "recreate");
    tout = new TTree("SmEWEventNtup", "SmEWEventNtup");

    // central jets
    tout->Branch("jet_n", &jet_n, "jet_n/I");
    tout->Branch("jet_eta", &jet_eta);
    tout->Branch("jet_phi", &jet_phi);
    tout->Branch("jet_m", &jet_m);
    tout->Branch("jet_pt", &jet_pt);
    tout->Branch("jet_E", &jet_E);

    // central + forward jets (eta < 5.0)
    tout->Branch("jet_vv_n", &jet_vv_n, "jet_vv_n/I");
    tout->Branch("jet_vv_eta", &jet_vv_eta);
    tout->Branch("jet_vv_phi", &jet_vv_phi);
    tout->Branch("jet_vv_m", &jet_vv_m);
    tout->Branch("jet_vv_pt", &jet_vv_pt);
    tout->Branch("jet_vv_E", &jet_vv_E);

    // "fat" jets
    tout->Branch("jet_akt10_n", &jet_akt10_n, "jet_akt10_n/I");
    tout->Branch("jet_akt10_eta", &jet_akt10_eta);
    tout->Branch("jet_akt10_phi", &jet_akt10_phi);
    tout->Branch("jet_akt10_m", &jet_akt10_m);
    tout->Branch("jet_akt10_pt", &jet_akt10_pt);
    tout->Branch("jet_akt10_E", &jet_akt10_E);

    // electrons
    tout->Branch("el_n", &el_n, "el_n/I");
    tout->Branch("el_charge", &el_charge);
    tout->Branch("el_eta", &el_eta);
    tout->Branch("el_phi", &el_phi);
    tout->Branch("el_m", &el_m);
    tout->Branch("el_pt", &el_pt);
    tout->Branch("el_E", &el_E);
    tout->Branch("el_trigger", &el_trigger);

    // muons
    tout->Branch("mu_n", &mu_n, "mu_n/I");
    tout->Branch("mu_charge", &mu_charge);
    tout->Branch("mu_eta", &mu_eta);
    tout->Branch("mu_phi", &mu_phi);
    tout->Branch("mu_m", &mu_m);
    tout->Branch("mu_pt", &mu_pt);
    tout->Branch("mu_E", &mu_E);
    tout->Branch("mu_trigger", &mu_trigger);

    // met & sumet
    tout->Branch("met", &met, "met/F");
    tout->Branch("met_phi", &met_phi, "met_phi/F");
    tout->Branch("sumet", &sumet, "sumet/F");

    // invariant mass of hard scatter
    tout->Branch("truth_mass", &truth_mass, "truth_mass/F");

    tout->Branch("ht", &ht, "ht/F");
    tout->Branch("ht_met", &ht_met, "ht_met/F");
    tout->Branch("ht2j", &ht2j, "ht2j/F");
    tout->Branch("ht2j_met", &ht2j_met, "ht2j_met/F");
    tout->Branch("ht_lep", &ht_lep, "ht_lep/F");
    tout->Branch("ht_lep_met", &ht_lep_met, "ht_lep_met/F");

    tout->Branch("m_ll", &m_ll, "m_ll/F");
    tout->Branch("vv_max_dijet_mass", &vv_max_dijet_mass, "vv_max_dijet_mass/F");
    tout->Branch("vv_leading_dijet_mass", &vv_leading_dijet_mass, "vv_leading_dijet_mass/F");
    tout->Branch("vv_leading_dijet_deta", &vv_leading_dijet_deta, "vv_leading_dijet_deta/F");
    tout->Branch("vv_4body_mass", &vv_4body_mass, "vv_4body_mass/F");

    tout->Branch("truth_mass", &truth_mass, "truth_mass/F");

    return;
}

Bool_t EventSelector::Process(Long64_t entry) {
    // The Process() function is called for each entry in the tree (or possibly
    // keyed object in the case of PROOF) to be processed. The entry argument
    // specifies which entry in the currently loaded tree is to be processed.
    // It can be passed to either EventSelector::GetEntry() or TBranch::GetEntry()
    // to read either all or the required parts of the data. When processing
    // keyed objects with PROOF, the object is already loaded and is available
    // via the fObject pointer.
    //
    // This function should contain the "body" of the analysis. It can contain
    // simple or elaborate selection criteria, run algorithms on the data
    // of the event and typically fill histograms.
    //
    // The processing can be stopped by calling Abort().
    //
    // Use fStatus to set the return value of TTree::Process().
    //
    // The return value is currently not used.

    GetEntry(entry);

    if (!(entry % 100))
        cout << "Processing event " << entry << "." << endl;

    Reset();

    // fill event-wide variables.
    FillMll();
    FillVVMasses();
    FillHT();

    pass_cuts = 1;

    if (pass_cuts)
        tout->Fill();

    return kTRUE;
}

void EventSelector::SlaveTerminate() {
    // The SlaveTerminate() function is called after all entries or objects
    // have been processed. When running with PROOF SlaveTerminate() is called
    // on each slave server.

    GetOutputList()->Add(tout);

    return;
}

void EventSelector::Terminate() {
    // The Terminate() function is the last function to be called during
    // a query. It always runs on the client, it can be used to present
    // the results graphically or save the results to file.

    TList *l = GetOutputList();
    size_t imax = l->GetSize();
    for (size_t i = 0; i < imax; i++)
        l->At(i)->Write();

    fout->Close();

    return;
}

void EventSelector::FillHT() {

    for (int i = 0; i < el_n; i++)
        ht_lep += el_pt->at(i);
    for (int i = 0; i < mu_n; i++)
        ht_lep += mu_pt->at(i);

    ht = ht_lep;
    ht2j = ht;
    for (int i = 0; i < jet_n; i++) {
        ht += jet_pt->at(i);
        if (i < 2)
            ht2j += jet_pt->at(i);
    }

    ht_lep_met = met + ht;
    ht_met = met + ht;
    ht2j_met = met + ht2j;
}

void EventSelector::FillMll() {
    TLorentzVector tlvl0, tlvl1;

    if (el_n > 1 && mu_n == 0) {
        tlvl0.SetPtEtaPhiE(el_pt->at(0), el_eta->at(0),
                el_phi->at(0), el_E->at(0));
        tlvl1.SetPtEtaPhiE(el_pt->at(1), el_eta->at(1),
                el_phi->at(1), el_E->at(1));
        m_ll = (tlvl0 + tlvl1).M();
    } else if (el_n && mu_n) {
        tlvl0.SetPtEtaPhiE(el_pt->at(0), el_eta->at(0),
                el_phi->at(0), el_E->at(0));
        tlvl1.SetPtEtaPhiE(mu_pt->at(0), mu_eta->at(0),
                mu_phi->at(0), mu_E->at(0));
        m_ll = (tlvl0 + tlvl1).M();
    } else if (mu_n > 1 && el_n == 0) {
        tlvl0.SetPtEtaPhiE(mu_pt->at(0), mu_eta->at(0),
                mu_phi->at(0), mu_E->at(0));
        tlvl1.SetPtEtaPhiE(mu_pt->at(1), mu_eta->at(1),
                mu_phi->at(1), mu_E->at(1));
        m_ll = (tlvl0 + tlvl1).M();
    } else {
        m_ll = 0;
    }

    return;
}


void EventSelector::FillVVMasses() {
    if (jet_vv_n < 2)
        return;

    TLorentzVector tlvj1, tlvj2;
    float mjj;
    int j1_idx = -1, j2_idx = -1;

    int njets = jet_vv_n;
    for (int iJet = 0; iJet < njets; iJet++) {

        if (j1_idx < 0) {
            j1_idx = iJet;
        } else if (jet_vv_pt->at(iJet) >
                jet_vv_pt->at(j1_idx)) {
            j2_idx = j1_idx;
            j1_idx = iJet;
        } else if (j2_idx < 0) {
            j2_idx = iJet;
        } else if (jet_vv_pt->at(iJet) >
                jet_vv_pt->at(j2_idx)) {
            j2_idx = iJet;
        }

        tlvj1.SetPtEtaPhiE(jet_vv_pt->at(iJet),
                jet_vv_eta->at(iJet), jet_vv_phi->at(iJet),
                jet_vv_E->at(iJet));

        for (int jJet = iJet+1; jJet < njets; jJet++) {
            tlvj2.SetPtEtaPhiE(jet_vv_pt->at(jJet),
                    jet_vv_eta->at(jJet), jet_vv_phi->at(jJet),
                    jet_vv_E->at(jJet));

            mjj = (tlvj1 + tlvj2).M();
            if (mjj > vv_max_dijet_mass)
                vv_max_dijet_mass = mjj;
        }
    }

    tlvj1.SetPtEtaPhiE(jet_vv_pt->at(j1_idx),
            jet_vv_eta->at(j1_idx), jet_vv_phi->at(j1_idx),
            jet_vv_E->at(j1_idx));

    tlvj2.SetPtEtaPhiE(jet_vv_pt->at(j2_idx),
            jet_vv_eta->at(j2_idx), jet_vv_phi->at(j2_idx),
            jet_vv_E->at(j2_idx));

    TLorentzVector tlvjets = tlvj1 + tlvj2;
    vv_leading_dijet_mass = tlvjets.M();
    vv_leading_dijet_deta = TMath::Abs(tlvj1.Eta() - tlvj2.Eta());

    TLorentzVector tlvl0, tlvl1;
    if (el_n > 1 && mu_n == 0) {
        tlvl0.SetPtEtaPhiE(el_pt->at(0), el_eta->at(0),
                el_phi->at(0), el_E->at(0));
        tlvl1.SetPtEtaPhiE(el_pt->at(1), el_eta->at(1),
                el_phi->at(1), el_E->at(1));

        vv_4body_mass = (tlvl0 + tlvl1 + tlvjets).M();
    } else if (el_n && mu_n) {
        tlvl0.SetPtEtaPhiE(el_pt->at(0), el_eta->at(0),
                el_phi->at(0), el_E->at(0));
        tlvl1.SetPtEtaPhiE(mu_pt->at(0), mu_eta->at(0),
                mu_phi->at(0), mu_E->at(0));

        vv_4body_mass = (tlvl0 + tlvl1 + tlvjets).M();
    } else if (mu_n > 1 && el_n == 0) {
        tlvl0.SetPtEtaPhiE(mu_pt->at(0), mu_eta->at(0),
                mu_phi->at(0), mu_E->at(0));
        tlvl1.SetPtEtaPhiE(mu_pt->at(1), mu_eta->at(1),
                mu_phi->at(1), mu_E->at(1));

        vv_4body_mass = (tlvl0 + tlvl1 + tlvjets).M();
    } else {
        vv_4body_mass = 0.0;
    }

    
    return;
}

void EventSelector::Reset() {
    ht = 0;
    ht2j = 0;
    ht_met = 0;
    ht2j_met = 0;
    ht_lep = 0;
    ht_lep_met = 0;

    vv_max_dijet_mass = 0;
    vv_leading_dijet_mass = 0;
    vv_leading_dijet_deta = 0;
    vv_4body_mass = 0;

    return;
}
