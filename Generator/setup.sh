export LHAPDF=~/Programming/lhapdf
export FASTJET=~/Programming/fastjet
export PYTHIA8=~/Programming/pythia
export PYTHIA8DATA=$PYTHIA8/xmldoc

export LD_LIBRARY_PATH=$LHAPDF/lib:$PYTHIA8/lib:$FASTJET/lib:$LD_LIBRARY_PATH
# for osx
# export DYLD_LIBRARY_PATH=$LHAPDF/lib:$PYTHIA8/lib:$FASTJET/lib:$LD_LIBRARY_PATH
