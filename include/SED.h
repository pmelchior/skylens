#ifndef SKYLENS_SED_H
#define SKYLENS_SED_H

#include "Filter.h"

namespace skylens {
  
  /// Class for spectral energy distributions.
  /// SEDs are maps from frequency to emission used to describe object 
  /// spectral flux.
  /// 
  /// Derived from Matthias Bartelmann's libastro.
  class SED : public Filter {
  public:
    SED();
    SED(const std::string& filename, double threshold=1e-5);
    void shift(double z);
    double kCorrection(const Filter& f, double z) const;
    double color(const Filter& f1, const Filter& f2) const;
  };

} // end namespace

#endif
