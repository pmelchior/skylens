#include "../include/SED.h"
#include <shapelens/FITS.h>
#include <math.h>

namespace skylens {

  using shapelens::FITS;

  SED::SED() : Filter() { }
  SED::SED(const std::string& filename, double threshold) : Filter() {
  fitsfile* fptr = FITS::openTable(filename);
    long nrows = FITS::getTableRows(fptr);
    int nu_col = FITS::getTableColumnNumber(fptr, "frequency");
    int fn_col = FITS::getTableColumnNumber(fptr, "spectral_flux");
    double nu, fn, max_fn = 0;
    for (long i=0; i < nrows; i++) {
      FITS::readTableValue(fptr, i, nu_col, nu);
      FITS::readTableValue(fptr, i, fn_col, fn);
      Filter::insert(std::pair<double, double>(nu, fn));
      if (fn > max_fn)
	max_fn = fn;
    }
    removeZeros(max_fn*threshold);
  }

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
