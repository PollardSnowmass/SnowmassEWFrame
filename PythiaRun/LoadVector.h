#include<vector>
#ifdef __CINT__
#pragma link C++ class vector<vector<int> >;
#pragma link C++ class vector<float>;
#else
template class std::vector<std::vector<int> >;
template class std::vector<float>;
#endif
