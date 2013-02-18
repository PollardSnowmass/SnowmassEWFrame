#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "Pythia.h"
#include <vector>
#include <iostream>
#include "fastjet/ClusterSequence.hh"

using namespace fastjet;
using namespace Pythia8;

bool endsWith (std::string const &fullString, std::string const &ending) {
    if (fullString.length() >= ending.length()) {
        return (0 == fullString.compare (fullString.length() - ending.length(), ending.length(), ending));
    } else {
        return false;
    }
}

int main(int argc, char *argv[]) {

    gROOT->ProcessLine(".L LoadVector.h++");

    if (argc != 3) {
        cout << "Usage: " << argv[0] << " cmnd_or_lhe_file outfile.root" << endl;
        exit(1);
    }

    cout << "Trying to open input file " << argv[argc-1] << endl;

    TFile *outf = TFile::Open(argv[argc-1], "recreate");
    if (!outf)
        return 0;

    TTree *t = new TTree("SmEWTruthNtup", "SmEWTruthNtup");

    std::vector<TLorentzVector> *tlv_bquarks = new std::vector<TLorentzVector>();

    int mc_n = 0;
    std::vector<float> *mc_E = new std::vector<float>();
    std::vector<float> *mc_pt = new std::vector<float>();
    std::vector<float> *mc_m = new std::vector<float>();
    std::vector<float> *mc_eta = new std::vector<float>();
    std::vector<float> *mc_phi = new std::vector<float>();
    std::vector<int> *mc_status = new std::vector<int>();
    std::vector<int> *mc_barcode = new std::vector<int>();
    std::vector<std::vector<int> > *mc_parents = new std::vector<std::vector<int> >();
    std::vector<std::vector<int> > *mc_children = new std::vector<std::vector<int> >();
    std::vector<int> *mc_pdgId = new std::vector<int>();
    std::vector<float> *mc_charge = new std::vector<float>();
    std::vector<std::vector<int> > *mc_child_index = new std::vector<std::vector<int> >();
    std::vector<std::vector<int> > *mc_parent_index = new std::vector<std::vector<int> >();

    int mcevt_n = 0;
    std::vector<int> *mcevt_signal_process_id = new std::vector<int>();
    std::vector<int> *mcevt_event_number = new std::vector<int>();
    std::vector<double> *mcevt_event_scale = new std::vector<double>();
    std::vector<double> *mcevt_alphaQCD = new std::vector<double>();
    std::vector<double> *mcevt_alphaQED = new std::vector<double>();
    std::vector<int> *mcevt_pdf_id1 = new std::vector<int>();
    std::vector<int> *mcevt_pdf_id2 = new std::vector<int>();
    std::vector<double> *mcevt_pdf_x1 = new std::vector<double>();
    std::vector<double> *mcevt_pdf_x2 = new std::vector<double>();
    std::vector<double> *mcevt_pdf_scale = new std::vector<double>();
    std::vector<double> *mcevt_pdf1 = new std::vector<double>();
    std::vector<double> *mcevt_pdf2 = new std::vector<double>();
    std::vector<std::vector<double> > *mcevt_weight = new std::vector<std::vector<double> >();

    std::vector<int> *el_indices = new std::vector<int>();
    std::vector<int> *mu_indices = new std::vector<int>();

    int jet_n = 0;
    std::vector<float> *jet_eta = new std::vector<float>();
    std::vector<float> *jet_phi = new std::vector<float>();
    std::vector<float> *jet_E = new std::vector<float>();
    std::vector<float> *jet_pt = new std::vector<float>();
    std::vector<float> *jet_m = new std::vector<float>();
    std::vector<int> *jet_bjet = new std::vector<int>();

    int jet_akt10_n = 0;
    std::vector<float> *jet_akt10_eta = new std::vector<float>();
    std::vector<float> *jet_akt10_phi = new std::vector<float>();
    std::vector<float> *jet_akt10_E = new std::vector<float>();
    std::vector<float> *jet_akt10_pt = new std::vector<float>();
    std::vector<float> *jet_akt10_m = new std::vector<float>();

    t->Branch("mc_n", &mc_n, "mc_n/I");
    t->Branch("mc_E", &mc_E);
    t->Branch("mc_pt", &mc_pt);
    t->Branch("mc_m", &mc_m);
    t->Branch("mc_eta", &mc_eta);
    t->Branch("mc_phi", &mc_phi);
    t->Branch("mc_status", &mc_status);
    t->Branch("mc_barcode", &mc_barcode);
    t->Branch("mc_parents", mc_parents);
    t->Branch("mc_children", &mc_children);
    t->Branch("mc_pdgId", &mc_pdgId);
    t->Branch("mc_charge", &mc_charge);
    t->Branch("mc_child_index", &mc_child_index);
    t->Branch("mc_parent_index", &mc_parent_index);

    t->Branch("mcevt_n", &mcevt_n, "mcevt_n/I");
    t->Branch("mcevt_signal_process_id", &mcevt_signal_process_id);
    t->Branch("mcevt_event_number", &mcevt_event_number);
    t->Branch("mcevt_event_scale", &mcevt_event_scale);
    t->Branch("mcevt_alphaQCD", &mcevt_alphaQCD);
    t->Branch("mcevt_alphaQED", &mcevt_alphaQED);
    t->Branch("mcevt_pdf_id1", &mcevt_pdf_id1);
    t->Branch("mcevt_pdf_id2", &mcevt_pdf_id2);
    t->Branch("mcevt_pdf_x1", &mcevt_pdf_x1);
    t->Branch("mcevt_pdf_x2", &mcevt_pdf_x2);
    t->Branch("mcevt_pdf_scale", &mcevt_pdf_scale);
    t->Branch("mcevt_pdf1", &mcevt_pdf1);
    t->Branch("mcevt_pdf2", &mcevt_pdf2);
    t->Branch("mcevt_weight", &mcevt_weight);

    t->Branch("el_indices", &el_indices);
    t->Branch("mu_indices", &mu_indices);

    t->Branch("jet_n", &jet_n, "jet_n/I");
    t->Branch("jet_eta", &jet_eta);
    t->Branch("jet_phi", &jet_phi);
    t->Branch("jet_E", &jet_E);
    t->Branch("jet_pt", &jet_pt);
    t->Branch("jet_m", &jet_m);
    t->Branch("jet_bjet", &jet_bjet);

    t->Branch("jet_akt10_n", &jet_akt10_n, "jet_akt10_n/I");
    t->Branch("jet_akt10_eta", &jet_akt10_eta);
    t->Branch("jet_akt10_phi", &jet_akt10_phi);
    t->Branch("jet_akt10_E", &jet_akt10_E);
    t->Branch("jet_akt10_pt", &jet_akt10_pt);
    t->Branch("jet_akt10_m", &jet_akt10_m);

    vector<PseudoJet> particles;
    TLorentzVector tlv_tmp;

    int pdgId;
    JetDefinition akt4_def(antikt_algorithm, 0.4);
    vector<PseudoJet> akt4_jets;

    JetDefinition akt10_def(antikt_algorithm, 1.0);
    vector<PseudoJet> akt10_jets;

    Pythia pythia;
    if (endsWith(string(argv[1]), string(".lhe")) ||
            endsWith(string(argv[1]), string(".lhef"))) {
        pythia.readFile("default.cmnd");
        pythia.init(argv[1]);
    } else {
        pythia.readFile(argv[1]);
        pythia.init();
    }

    int nEvent = pythia.mode("Main:numberOfEvents");

    // begin event loop.
    for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

        if (!(iEvent % 100))
            cout << "Generating event " << iEvent << "." << endl;

        // skip if error
        if (!pythia.next()) continue;

        mc_n = pythia.event.size();

        mc_E->resize(mc_n);
        mc_pt->resize(mc_n);
        mc_m->resize(mc_n);
        mc_eta->resize(mc_n);
        mc_phi->resize(mc_n);
        mc_status->resize(mc_n);
        mc_barcode->resize(mc_n);
        mc_parents->resize(mc_n);
        mc_children->resize(mc_n);
        mc_pdgId->resize(mc_n);
        mc_charge->resize(mc_n);
        mc_child_index->resize(mc_n);
        mc_parent_index->resize(mc_n);

        mcevt_n = 0;
        mcevt_signal_process_id->clear();
        mcevt_event_number->clear();
        mcevt_event_scale->clear();
        mcevt_alphaQCD->clear();
        mcevt_alphaQED->clear();
        mcevt_pdf_id1->clear();
        mcevt_pdf_id2->clear();
        mcevt_pdf_x1->clear();
        mcevt_pdf_x2->clear();
        mcevt_pdf_scale->clear();
        mcevt_pdf1->clear();
        mcevt_pdf2->clear();
        mcevt_weight->clear();

        el_indices->clear();
        mu_indices->clear();

        tlv_bquarks->clear();

        jet_n = 0;
        jet_eta->clear();
        jet_phi->clear();
        jet_E->clear();
        jet_pt->clear();
        jet_m->clear();
        jet_bjet->clear();

        jet_akt10_n = 0;
        jet_akt10_eta->clear();
        jet_akt10_phi->clear();
        jet_akt10_E->clear();
        jet_akt10_pt->clear();
        jet_akt10_m->clear();

        particles.clear();

        for (int iPart = 0; iPart < mc_n; iPart++) {
            mc_E->at(iPart) = (pythia.event[iPart].e()*1000);
            mc_pt->at(iPart) = (pythia.event[iPart].pT()*1000);
            mc_m->at(iPart) = (pythia.event[iPart].m()*1000);
            mc_eta->at(iPart) = (pythia.event[iPart].eta());
            mc_phi->at(iPart) = (pythia.event[iPart].phi());
            mc_status->at(iPart) = (pythia.event[iPart].status());
            mc_barcode->at(iPart) = (-999);
            mc_parents->at(iPart) = (pythia.event.motherList(iPart));
            mc_children->at(iPart) = (pythia.event.daughterList(iPart));
            mc_pdgId->at(iPart) = (pythia.event[iPart].id());
            mc_charge->at(iPart) = (pythia.event[iPart].charge());
            mc_child_index->at(iPart) = (pythia.event.daughterList(iPart));
            mc_parent_index->at(iPart) = (pythia.event.motherList(iPart));

            // set up jet clustering over outgoing partons.
            // straight from fastjet quickstart guide.

            pdgId = abs(pythia.event[iPart].id());
            if (pdgId == 11)
                el_indices->push_back(iPart);
            else if (pdgId == 13)
                mu_indices->push_back(iPart);

            // cluster outgoing particles except neutrinos, muons.
            // ignore eta > 3.0 particles.
            if (pythia.event[iPart].status() > 0 &&
                    fabs(pythia.event[iPart].eta()) < 7.0 &&
                    pdgId != 12 && pdgId != 13 &&
                    pdgId != 14 && pdgId != 16) {

                tlv_tmp.SetPtEtaPhiE(pythia.event[iPart].pT(),
                        pythia.event[iPart].eta(),
                        pythia.event[iPart].phi(),
                        pythia.event[iPart].e());

                particles.push_back(PseudoJet(tlv_tmp));
            }

            // check for bjets.
            if (pdgId == 5) {
                tlv_tmp.SetPtEtaPhiE(pythia.event[iPart].pT(),
                        pythia.event[iPart].eta(),
                        pythia.event[iPart].phi(),
                        pythia.event[iPart].e());

                tlv_bquarks->push_back(TLorentzVector(tlv_tmp));
            }
        }

        ClusterSequence akt4_cs(particles, akt4_def);
        akt4_jets = sorted_by_pt(akt4_cs.inclusive_jets());

        jet_n = akt4_jets.size();
        for (size_t i = 0; i < jet_n; i++) {
            jet_eta->push_back(akt4_jets[i].eta());
            jet_phi->push_back(akt4_jets[i].phi());
            jet_E->push_back(akt4_jets[i].E()*1000);
            jet_pt->push_back(akt4_jets[i].pt()*1000);
            jet_m->push_back(akt4_jets[i].m()*1000);

            tlv_tmp.SetPtEtaPhiE(akt4_jets[i].pt(),
                    akt4_jets[i].eta(),
                    akt4_jets[i].phi(),
                    akt4_jets[i].E());

            jet_bjet->push_back(0);
            for (size_t j = 0; j < tlv_bquarks->size(); j++) {
                if (tlv_tmp.DeltaR(tlv_bquarks->at(j)) < 0.4) {
                    jet_bjet->at(i) = 1;
                    break;
                }
            }
        }


        ClusterSequence akt10_cs(particles, akt10_def);
        akt10_jets = sorted_by_pt(akt10_cs.inclusive_jets());

        jet_akt10_n = akt10_jets.size();
        for (size_t i = 0; i < jet_akt10_n; i++) {
            jet_akt10_eta->push_back(akt10_jets[i].eta());
            jet_akt10_phi->push_back(akt10_jets[i].phi());
            jet_akt10_E->push_back(akt10_jets[i].E()*1000);
            jet_akt10_pt->push_back(akt10_jets[i].pt()*1000);
            jet_akt10_m->push_back(akt10_jets[i].m()*1000);
        }

        // done with jet clustering.

        mcevt_event_number->push_back(iEvent);
        mcevt_event_scale->push_back(pythia.event.scale());
        mcevt_pdf_id1->push_back(pythia.info.id1());
        mcevt_pdf_id2->push_back(pythia.info.id2());
        mcevt_pdf_x1->push_back(pythia.info.x1());
        mcevt_pdf_x2->push_back(pythia.info.x2());
        mcevt_pdf_scale->push_back(pythia.event.scale());
        mcevt_pdf1->push_back(pythia.info.pdf1());
        mcevt_pdf2->push_back(pythia.info.pdf2());
        mcevt_weight->push_back(std::vector<double>(1, 1.0));

        t->Fill();
    }

    t->Write();
    outf->Close();

    pythia.statistics(true);

    return 0;
}
