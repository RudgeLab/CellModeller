
%module SPP
%{
  #include "SPP.h"
%}

%naturalvar SPP::cell_centers;
%naturalvar SPP::cell_polarization;

%include "pyabc.i"
%include "std_vector.i"


%include "SPP.h"
%template (floatVector) std::vector<float>;
