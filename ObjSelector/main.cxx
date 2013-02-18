#include <vector>
#include <iostream>
#include <sstream>
#include <TChain.h>
#include <TProof.h>
#include <TSystem.h>
#include <TROOT.h>

using namespace std;

void split(const string& s, char c, vector<string>& v) {
    string::size_type i = 0;
    string::size_type j = s.find(c);
    while (j != string::npos) {
        v.push_back(s.substr(i, j-i));
        i = ++j;
        j = s.find(c, j);
    }
    v.push_back(s.substr(i, s.length()));
}

int main(int argc, char *argv[]) {
    if (argc != 3) {
        cout << "Usage: " << argv[0] << " comma,separated,infiles outfile." << endl;
        return 1;
    }

    TChain *c = new TChain("SmEWTruthNtup");
    string line;
    vector<string> v;
    int l;

    split((string) argv[1], ',', v);
    l = v.size();

    cout << l << " File" << ((l==1)?"": "s") << " to Analyze." << endl;

    // build the chain.
    for (int i = 0; i < l; i++) {
        c->AddFile(v.at(i).c_str());
    }

    // instantiate new ObjSelector object and loop.
    // TSystem *gSystem = new TSystem();
    // TProof::Open(gSystem->GetFromPipe("pod-info -c"));
    c->Process("ObjSelector.cxx+", argv[2]);

    delete c;
    // delete p;

    return 0;
}
