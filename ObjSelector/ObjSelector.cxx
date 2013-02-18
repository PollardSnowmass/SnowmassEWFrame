#define ObjSelector_cxx

#include "ObjSelector.h"
#include "LoadVector.h"
#include <iostream>
#include <TLorentzVector.h>
#include "config.h"

void ObjSelector::Begin(TTree * /*tree*/) {
    // The Begin() function is called at the start of the query.
    // When running with PROOF Begin() is only called on the client.
    // The tree argument is deprecated (on PROOF 0 is passed).

    return;
}

void ObjSelector::SlaveBegin(TTree * /*tree*/) {
    // The SlaveBegin() function is called after the Begin() function.
    // When running with PROOF SlaveBegin() is called on each slave server.
    // The tree argument is deprecated (on PROOF 0 is passed).

    TString fname = GetOption();
    fout = new TFile(fname, "recreate");
    tout = new TTree("SmEWObjNtup", "SmEWObjNtup");

    myRand = new TRandom3();

    ReadConfig();

    // central jets
    tout->Branch("jet_n", &jet_n_out, "jet_n/I");
    tout->Branch("jet_eta", &jet_eta_out);
    tout->Branch("jet_phi", &jet_phi_out);
    tout->Branch("jet_m", &jet_m_out);
    tout->Branch("jet_pt", &jet_pt_out);
    tout->Branch("jet_E", &jet_E_out);

    // central + forward jets (eta < 5.0)
    tout->Branch("jet_vv_n", &jet_vv_n_out, "jet_vv_n/I");
    tout->Branch("jet_vv_eta", &jet_vv_eta_out);
    tout->Branch("jet_vv_phi", &jet_vv_phi_out);
    tout->Branch("jet_vv_m", &jet_vv_m_out);
    tout->Branch("jet_vv_pt", &jet_vv_pt_out);
    tout->Branch("jet_vv_E", &jet_vv_E_out);

    // "fat" jets
    tout->Branch("jet_akt10_n", &jet_akt10_n_out, "jet_akt10_n/I");
    tout->Branch("jet_akt10_eta", &jet_akt10_eta_out);
    tout->Branch("jet_akt10_phi", &jet_akt10_phi_out);
    tout->Branch("jet_akt10_m", &jet_akt10_m_out);
    tout->Branch("jet_akt10_pt", &jet_akt10_pt_out);
    tout->Branch("jet_akt10_E", &jet_akt10_E_out);

    // electrons
    tout->Branch("el_n", &el_n_out, "el_n/I");
    tout->Branch("el_charge", &el_charge_out);
    tout->Branch("el_eta", &el_eta_out);
    tout->Branch("el_phi", &el_phi_out);
    tout->Branch("el_m", &el_m_out);
    tout->Branch("el_pt", &el_pt_out);
    tout->Branch("el_E", &el_E_out);
    tout->Branch("el_trigger", &el_trigger_out);

    // muons
    tout->Branch("mu_n", &mu_n_out, "mu_n/I");
    tout->Branch("mu_charge", &mu_charge_out);
    tout->Branch("mu_eta", &mu_eta_out);
    tout->Branch("mu_phi", &mu_phi_out);
    tout->Branch("mu_m", &mu_m_out);
    tout->Branch("mu_pt", &mu_pt_out);
    tout->Branch("mu_E", &mu_E_out);
    tout->Branch("mu_trigger", &mu_trigger_out);

    // met & sumet
    tout->Branch("met", &met_out, "met/F");
    tout->Branch("met_phi", &met_phi_out, "met_phi/F");
    tout->Branch("sumet", &sumet_out, "sumet/F");

    // invariant mass of hard scatter
    tout->Branch("truth_mass", &truth_mass, "truth_mass/F");

    return;
}

Bool_t ObjSelector::Process(Long64_t entry) {
    // The Process() function is called for each entry in the tree (or possibly
    // keyed object in the case of PROOF) to be processed. The entry argument
    // specifies which entry in the currently loaded tree is to be processed.
    // It can be passed to either ObjSelector::GetEntry() or TBranch::GetEntry()
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

    if (do_reco) {
        SmearElectrons();
        SmearMuons();
        SmearJets();
        SmearAkt10Jets();
    }

    // fill good particles.
    FillGoodElectrons();
    FillGoodMuons();
    FillGoodJets();
    FillGoodVVJets();
    FillGoodAkt10Jets();

    // deal with jet/electron overlap.
    ElJetOverlap();
    JetElOverlap();

    // get trigger information.
    FillTrigger();

    // fill event-wide variables.
    FillSumEt();
    FillMET();

    // get the mass of the resonant particle.
    FillTruthMass();

    tout->Fill();

    return kTRUE;
}

void ObjSelector::SlaveTerminate() {
    // The SlaveTerminate() function is called after all entries or objects
    // have been processed. When running with PROOF SlaveTerminate() is called
    // on each slave server.

    GetOutputList()->Add(tout);

    return;
}

void ObjSelector::Terminate() {
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

void ObjSelector::FillGoodElectrons() {
    int iEl;
    for (int i = 0; i < el_indices->size(); i++) {
        iEl = el_indices->at(i);
        if (IsGoodElectron(iEl)) {
            el_n_out++;
            el_charge_out.push_back(mc_charge->at(iEl));
            el_eta_out.push_back(mc_eta->at(iEl));
            el_phi_out.push_back(mc_phi->at(iEl));
            el_pt_out.push_back(mc_pt->at(iEl));
            el_E_out.push_back(mc_E->at(iEl));
            el_m_out.push_back(mc_m->at(iEl));
        }
    }

    return;
}

void ObjSelector::FillGoodMuons() {
    int iMu;
    for (int i = 0; i < mu_indices->size(); i++) {
        iMu = mu_indices->at(i);
        if (IsGoodMuon(iMu)) {
            mu_n_out++;
            mu_charge_out.push_back(mc_charge->at(iMu));
            mu_eta_out.push_back(mc_eta->at(iMu));
            mu_phi_out.push_back(mc_phi->at(iMu));
            mu_pt_out.push_back(mc_pt->at(iMu));
            mu_E_out.push_back(mc_E->at(iMu));
            mu_m_out.push_back(mc_m->at(iMu));
        }
    }

    return;
}

void ObjSelector::FillGoodJets() {
    for (int iJet = 0; iJet < jet_n; iJet++) {
        if (IsGoodJet(iJet)) {
            jet_n_out++;
            jet_eta_out.push_back(jet_eta->at(iJet));
            jet_phi_out.push_back(jet_phi->at(iJet));
            jet_pt_out.push_back(jet_pt->at(iJet));
            jet_E_out.push_back(jet_E->at(iJet));
            jet_m_out.push_back(jet_m->at(iJet));
        }
    }

    return;
}

void ObjSelector::FillGoodVVJets() {
    has_bjet = 0;
    for (int iJet = 0; iJet < jet_n; iJet++) {
        if (IsGoodVVJet(iJet)) {
            jet_vv_n_out++;
            jet_vv_eta_out.push_back(jet_eta->at(iJet));
            jet_vv_phi_out.push_back(jet_phi->at(iJet));
            jet_vv_pt_out.push_back(jet_pt->at(iJet));
            jet_vv_E_out.push_back(jet_E->at(iJet));
            jet_vv_m_out.push_back(jet_m->at(iJet));
        }
    }

    return;
}

void ObjSelector::FillGoodAkt10Jets() {
    for (int iJet = 0; iJet < jet_akt10_n; iJet++) {
        if (IsGoodAkt10Jet(iJet)) {
            jet_akt10_n_out++;
            jet_akt10_eta_out.push_back(jet_akt10_eta->at(iJet));
            jet_akt10_phi_out.push_back(jet_akt10_phi->at(iJet));
            jet_akt10_pt_out.push_back(jet_akt10_pt->at(iJet));
            jet_akt10_E_out.push_back(jet_akt10_E->at(iJet));
            jet_akt10_m_out.push_back(jet_akt10_m->at(iJet));
        }
    }

    return;
}

bool ObjSelector::IsGoodElectron(int iEl) {
    if (mc_status->at(iEl) <= 0)
        return 0;

    int abspdgId = TMath::Abs(mc_pdgId->at(iEl));
    if (abspdgId != 11)
        return 0;

    if (mc_pt->at(iEl) < el_pt_cut)
        return 0;

    float abseta = TMath::Abs(mc_eta->at(iEl));
    if (abseta > el_eta_cut)
        return 0;

    if (do_reco &&
            !PassEfficiency(myRand, 0.85-0.191*TMath::Exp(1-mc_pt->at(iEl)/20000.0)))
        return 0;

    return 1;
}

bool ObjSelector::IsGoodMuon(int iMu) {
    if (mc_status->at(iMu) <= 0)
        return 0;

    int abspdgId = TMath::Abs(mc_pdgId->at(iMu));
    if (abspdgId != 13)
        return 0;

    float abseta = TMath::Abs(mc_eta->at(iMu));
    if (abseta > mu_eta_cut)
        return 0;

    if (mc_pt->at(iMu) < mu_pt_cut)
        return 0;

    float mu_eff = (mu_eff_corr ? sqrt(0.45) : 0.97)
        - 0.03*mc_pt->at(iMu)/1000000.0;

    if (do_reco && !PassEfficiency(myRand, mu_eff))
        return 0;

    return 1;
}

bool ObjSelector::IsGoodJet(int iJet) {
    float abseta = TMath::Abs(jet_eta->at(iJet));
    if (abseta > jet_eta_cut)
        return 0;

    if (jet_pt->at(iJet) < jet_pt_cut)
        return 0;

    return 1;
}

bool ObjSelector::IsGoodVVJet(int iJet) {
    float abseta = TMath::Abs(jet_eta->at(iJet));
    if (abseta > jet_vv_eta_cut)
        return 0;

    float btag_eff;
    if (jet_bjet->at(iJet) && jet_pt->at(iJet) > jet_bjet_pt_thresh) {
        btag_eff = 0.0;
        if (dobtag_veto) {
            btag_eff = 0.5+0.03*jet_pt->at(iJet)/1000000.0;
            if (abseta > 2.5) btag_eff =- 0.3*abseta;
            if (btag_eff<0.0) btag_eff = 0.0;
        }

        if (do_reco && PassEfficiency(myRand, btag_eff))
            has_bjet = 1;
    }

    if (jet_pt->at(iJet) < jet_vv_pt_cut)
        return 0;

    return 1;
}

bool ObjSelector::IsGoodAkt10Jet(int iJet) {
    float abseta = TMath::Abs(jet_akt10_eta->at(iJet));
    if (abseta > jet_akt10_eta_cut)
        return 0;

    if (jet_akt10_pt->at(iJet) < jet_akt10_pt_cut)
        return 0;

    if (jet_akt10_m->at(iJet) < jet_akt10_m_cut)
        return 0;

    return 1;
}

void ObjSelector::FillSumEt() {
    int abspdgId;
    float abseta;
    float sumet = 0;
    for (int iPart = 0; iPart < mc_n; iPart++) {
        if (mc_status->at(iPart) <= 0)
            continue;

        abspdgId = TMath::Abs(mc_pdgId->at(iPart));
        if (abspdgId == 12 || abspdgId == 14 || abspdgId == 16)
            continue;

        abseta = TMath::Abs(mc_eta->at(iPart));
        if (abseta > 5)
            continue;

        sumet += mc_pt->at(iPart);
    }

    sumet_out = sumet;

    return;
}

void ObjSelector::FillTrigger()  {
    if (do_reco) {
        for (int iEl = 0; iEl < el_n_out; iEl++) {
            // check if it triggers...
            el_trigger_out.push_back(PassEfficiency(myRand, 0.88));
        }

        for (int iMu = 0; iMu < mu_n_out; iMu++) {
            if (TMath::Abs(mu_eta_out[iMu]) < 1.0) {
                mu_trigger_out.push_back(PassEfficiency(myRand, 0.64));
            } else {
                mu_trigger_out.push_back(PassEfficiency(myRand, 0.86));
            }
        }
    } else {
        el_trigger_out.resize(el_n_out, 1);
        mu_trigger_out.resize(mu_n_out, 1);
    }

    return;
}

void ObjSelector::FillMET() {
    TLorentzVector tlv_met, tlv_tmp;

    tlv_met = TLorentzVector(0, 0, 0, 0);
    for (size_t i = 0; i < jet_eta_out.size(); i++) {
        tlv_tmp.SetPtEtaPhiE(jet_pt_out[i], jet_eta_out[i], jet_phi_out[i], jet_E_out[i]);
        tlv_met -= tlv_tmp;
    }

    for (size_t i = 0; i < el_eta_out.size(); i++) {
        tlv_tmp.SetPtEtaPhiE(el_pt_out[i], el_eta_out[i], el_phi_out[i], el_E_out[i]);
        tlv_met -= tlv_tmp;
    }

    for (size_t i = 0; i < mu_eta_out.size(); i++) {
        tlv_tmp.SetPtEtaPhiE(mu_pt_out[i], mu_eta_out[i], mu_phi_out[i], mu_E_out[i]);
        tlv_met -= tlv_tmp;
    }

    float met_x = tlv_met.X(), met_y = tlv_met.Y();
    float met_res, smear;
    if (do_reco) {
        // smear met.
        met_res = (0.40+0.09*sqrt(mu))*sqrt(sumet_out/1000.0 + mu*20);

        smear = Smear(myRand, met_res/(met_x/1000.0));
        met_x *= smear;
        smear = Smear(myRand, met_res/(met_y/1000.0));
        met_y *= smear;
    }

    met_out = sqrt(met_x*met_x + met_y*met_y);
    met_phi_out = TMath::ATan2(met_y, met_x);

    return;
}

void ObjSelector::FillTruthMass() {
    if (mc_m->size() > 5)
        truth_mass = mc_m->at(5);
    else
        truth_mass = 0;

    return;
}

void ObjSelector::ElJetOverlap() {
    // remove jets on top of electrons.

    TLorentzVector tlvel, tlvjet;

    float drmin, dr;
    int ijetmin = -1;
    for (int iEl = 0; iEl < el_n_out; iEl++) {
        tlvel.SetPtEtaPhiE(el_pt_out[iEl], el_eta_out[iEl], el_phi_out[iEl], el_E_out[iEl]);
        drmin = 999;
        for (int iJet = 0; iJet < jet_n_out; iJet++) {
            tlvjet.SetPtEtaPhiE(jet_pt_out[iJet], jet_eta_out[iJet], jet_phi_out[iJet], jet_E_out[iJet]);
            dr = tlvjet.DeltaR(tlvel);
            if (dr < drmin) {
                drmin = dr;
                ijetmin = iJet;
            }
        }

        if (drmin < 0.2) {
            jet_n_out--;
            jet_eta_out.erase(jet_eta_out.begin()+ijetmin);
            jet_phi_out.erase(jet_phi_out.begin()+ijetmin);
            jet_m_out.erase(jet_m_out.begin()+ijetmin);
            jet_pt_out.erase(jet_pt_out.begin()+ijetmin);
            jet_E_out.erase(jet_E_out.begin()+ijetmin);
        }
    }

    // same for vv jets.
    for (int iEl = 0; iEl < el_n_out; iEl++) {
        tlvel.SetPtEtaPhiE(el_pt_out[iEl], el_eta_out[iEl], el_phi_out[iEl], el_E_out[iEl]);
        drmin = 999;
        for (int iJet = 0; iJet < jet_vv_n_out; iJet++) {
            tlvjet.SetPtEtaPhiE(jet_vv_pt_out[iJet], jet_vv_eta_out[iJet], jet_vv_phi_out[iJet], jet_vv_E_out[iJet]);
            dr = tlvjet.DeltaR(tlvel);
            if (dr < drmin) {
                drmin = dr;
                ijetmin = iJet;
            }
        }

        if (drmin < 0.2) {
            jet_vv_n_out--;
            jet_vv_eta_out.erase(jet_vv_eta_out.begin()+ijetmin);
            jet_vv_phi_out.erase(jet_vv_phi_out.begin()+ijetmin);
            jet_vv_m_out.erase(jet_vv_m_out.begin()+ijetmin);
            jet_vv_pt_out.erase(jet_vv_pt_out.begin()+ijetmin);
            jet_vv_E_out.erase(jet_vv_E_out.begin()+ijetmin);
        }
    }

    return;
}

void ObjSelector::JetElOverlap() {
    // remove electrons in jets.

    TLorentzVector tlvel, tlvjet;

    float drmin, dr;
    for (int iJet = 0; iJet < jet_n_out; iJet++) {
        tlvjet.SetPtEtaPhiE(jet_pt_out[iJet], jet_eta_out[iJet], jet_phi_out[iJet], jet_E_out[iJet]);
        for (int iEl = 0; iEl < el_n_out; iEl++) {
            tlvel.SetPtEtaPhiE(el_pt_out[iEl], el_eta_out[iEl], el_phi_out[iEl], el_E_out[iEl]);
            dr = tlvjet.DeltaR(tlvel);
            if (dr < 0.2) {
                el_n_out--;
                el_charge_out.erase(el_charge_out.begin()+iEl);
                el_eta_out.erase(el_eta_out.begin()+iEl);
                el_phi_out.erase(el_phi_out.begin()+iEl);
                el_pt_out.erase(el_pt_out.begin()+iEl);
                el_E_out.erase(el_E_out.begin()+iEl);
                iEl--;
            }
        }
    }

    return;
}

float ObjSelector::quadsum(const float &x, const float &y) {
    return sqrt(x*x + y*y);
}

void ObjSelector::SmearElectrons() {
    if (!do_reco)
        return;

    float smear;
    int iEl;
    for (int i = 0; i < el_indices->size(); i++) {
        iEl = el_indices->at(i);
        if (mc_status->at(iEl) < 0)
            continue;

        // calorimeter saturation
        if (mc_E->at(iEl) > 6000000)
            mc_E->at(iEl) = 6000000;

        if (TMath::Abs(mc_pdgId->at(iEl)) == 11) {
            smear = Smear(myRand,
                    quadsum(0.125/TMath::Sqrt(mc_E->at(iEl)/1000.), 0.0125));

            mc_E->at(iEl) = mc_E->at(iEl)*smear;
            mc_pt->at(iEl) = mc_pt->at(iEl)*smear;
        }
    }

    return;
}

void ObjSelector::SmearMuons() {
    if (!do_reco)
        return;

    float smear;
    int iMu;
    float q, pt, qbypt, sigma_qbypt;
    for (int i = 0; i < mu_indices->size(); i++) {
        iMu = mu_indices->at(i);
        if (mc_status->at(iMu) < 0)
            continue;

        if (TMath::Abs(mc_pdgId->at(iMu)) == 13) {
            q = mc_charge->at(iMu);
            pt = mc_pt->at(iMu);
            qbypt = q/pt;

            // get pt resolution from ES muons tool.
            sigma_qbypt = qbypt*(0.35*pt/3000000);

            // get value by which to smear q/pt.
            smear = Smear(myRand, TMath::Abs(sigma_qbypt/qbypt));

            // smear q/pt.
            qbypt *= smear;

            // re-evaluate pt, charge.
            pt = 1.0/TMath::Abs(qbypt);
            q = qbypt > 0 ? 1 : -1;

            // store new charge.
            mc_charge->at(iMu) = q;

            // smear pt and energy.
            smear = pt/mc_pt->at(iMu);
            mc_pt->at(iMu) = mc_pt->at(iMu)*smear;
            mc_E->at(iMu) = mc_E->at(iMu)*smear;
        }
    }

    return;
}

void ObjSelector::SmearJets() {
    if (!do_reco)
        return;

    float smear;
    float jet_res;
    for (int iJet = 0; iJet < jet_n; iJet++) {
        if (TMath::Abs(jet_eta->at(iJet)) > 5.0)
            continue;

        jet_res = quadsum(1.1/TMath::Sqrt(jet_E->at(iJet)/1000.0), 0.05);
        smear = Smear(myRand, jet_res);

        jet_E->at(iJet) = jet_E->at(iJet)*smear;
        jet_pt->at(iJet) = jet_pt->at(iJet)*smear;
    }
}

void ObjSelector::SmearAkt10Jets() {
    if (!do_reco)
        return;

    float smear;
    float jet_akt10_res;
    TLorentzVector tlv;
    for (int iJet = 0; iJet < jet_akt10_n; iJet++) {
        if (TMath::Abs(jet_akt10_eta->at(iJet)) > 5.0)
            continue;

        jet_akt10_res = quadsum(1.1/TMath::Sqrt(jet_akt10_E->at(iJet)/1000.0), 0.05);
        smear = Smear(myRand, jet_akt10_res);

        jet_akt10_E->at(iJet) = jet_akt10_E->at(iJet)*smear;
        jet_akt10_pt->at(iJet) = jet_akt10_pt->at(iJet)*smear;

        tlv.SetPtEtaPhiE(jet_akt10_pt->at(iJet), jet_akt10_eta->at(iJet),
                jet_akt10_phi->at(iJet), jet_akt10_E->at(iJet));

        jet_akt10_m->at(iJet) = tlv.M();
    }
}

float ObjSelector::Smear(TRandom3 *myRand1, float width) {
    return myRand1->Gaus(1, width);
}

bool ObjSelector::PassEfficiency(TRandom3 *myRand1, float efficiency) {
    return myRand1->Uniform() < efficiency;
}

void ObjSelector::Reset() {
    sumet_out = 0;

    jet_n_out = 0;
    jet_eta_out.clear();
    jet_phi_out.clear();
    jet_m_out.clear();
    jet_pt_out.clear();
    jet_E_out.clear();

    jet_vv_n_out = 0;
    jet_vv_eta_out.clear();
    jet_vv_phi_out.clear();
    jet_vv_m_out.clear();
    jet_vv_pt_out.clear();
    jet_vv_E_out.clear();

    jet_akt10_n_out = 0;
    jet_akt10_eta_out.clear();
    jet_akt10_phi_out.clear();
    jet_akt10_m_out.clear();
    jet_akt10_pt_out.clear();
    jet_akt10_E_out.clear();

    el_n_out = 0;
    el_charge_out.clear();
    el_eta_out.clear();
    el_phi_out.clear();
    el_m_out.clear();
    el_pt_out.clear();
    el_E_out.clear();
    el_trigger_out.clear();

    mu_n_out = 0;
    mu_charge_out.clear();
    mu_eta_out.clear();
    mu_phi_out.clear();
    mu_m_out.clear();
    mu_pt_out.clear();
    mu_E_out.clear();
    mu_trigger_out.clear();

    met_out = 0;
    met_phi_out = 0;

    truth_mass = 0;

    return;
}

void ObjSelector::ReadConfig() {
    mu = USConfig::mu;
    do_reco = USConfig::do_reco;
    dobtag_veto = USConfig::dobtag_veto;

    el_pt_cut = USConfig::el_pt_cut;
    el_eta_cut = USConfig::el_eta_cut;

    mu_pt_cut = USConfig::mu_pt_cut;
    mu_eta_cut = USConfig::mu_eta_cut;
    mu_eff_corr = USConfig::mu_eff_corr;

    jet_pt_cut = USConfig::jet_pt_cut;
    jet_eta_cut = USConfig::jet_eta_cut;
    jet_bjet_pt_thresh = USConfig::jet_bjet_pt_thresh;

    jet_vv_pt_cut = USConfig::jet_vv_pt_cut;
    jet_vv_eta_cut = USConfig::jet_vv_eta_cut;

    jet_akt10_pt_cut = USConfig::jet_akt10_pt_cut;
    jet_akt10_eta_cut = USConfig::jet_akt10_eta_cut;
    jet_akt10_m_cut = USConfig::jet_akt10_m_cut;

    return;
}
