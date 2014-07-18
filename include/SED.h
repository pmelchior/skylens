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
    /// Default constructor.
    SED();
    /// Argumented constructor to read a FITS table from \p filename.
    /// The table should have a \p frequency column (in units of \f$10^{15}\f$ 
    /// Hz) and a \p spectral_flux column 
    /// (\f$f\f$: typically in units of erg/Hz). 
    /// To speed up computation, only entries with \f[
    /// f > \mathtt{threshold} * \max(f)
    /// \f]
    /// are considered; no ordering of the frequency array is needed.
    SED(const std::string& filename, double threshold=1e-3);
    /// Redshift the emission.
    /// Subsequent calls to this function do not retain memory.
    void shift(double z);
    /// Compute k-correction of the SED redshifted by \p z and seen in the
    /// rest-frame filter \p f. 
    double kCorrection(const Filter& f, double z) const;
    /// Returns the color of the SED between the filter curves \p f1 and \p f2.
    double color(const Filter& f1, const Filter& f2) const;
  };

} // end namespace

#endif
