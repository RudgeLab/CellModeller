
%module SPP
%{
  #include "SPP.h"
%}

%naturalvar SPP::cell_centers;
%naturalvar SPP::cell_directions;

%include "pyabc.i"
%include "std_vector.i"


namespace std {
  %template(VecFloat) vector<float>;
  %template(VecVecFloat) vector< vector<float> >;
}

%include "SPP.h"

