#include "../include/SED.h"
#include <math.h>

namespace skylens {

  SED::SED() : Filter() { }
  SED::SED(const std::string& filename, double threshold) :  Filter(filename, threshold) { }

  void SED::shift(double z) {
    Filter::z = z;
  }

  double SED::kCorrection(const Filter& f, double z) const {
    SED sa(*this), sb(*this);
    sa.shift(z);
    sa *= f;
    sb *= f;
    return -2.5*log10((1+z)*sa.computeNorm()/sb.computeNorm());
  }

  double SED::color(const Filter& f1, const Filter& f2) const {
    SED s1(*this), s2(*this);
    s1 *= f1;
    s2 *= f2;
    return -2.5*log10(s1.computeNorm()/s2.computeNorm()*f2.computeNorm()/f1.computeNorm());
  }
  

} // end namespace
