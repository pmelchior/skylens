#ifndef SKYLENS_FILTER_H
#define SKYLENS_FILTER_H

#include <map>
#include <string>

namespace skylens {
  
  /// Class for filter curves.
  /// Filters are maps from frequency to transmission efficiency 
  /// used to describe observation bands. 
  /// 
  /// Derived from Matthias Bartelmann's libastro,
  /// substantially modified to work as a \p std::map by Peter Melchior
  class Filter : public std::map<double, double> {
  public:
    /// Default constructor.
    Filter();
    /// Argumented constructor to read a FITS table from \p filename.
    /// The table should have a \p frequency column (in units of \f$10^{15}\f$ 
    /// Hz) and a \p transmission column (\f$\tau\f$: typically from 0 to 1). 
    /// To speed up computation, only entries with \f[
    /// \tau > \mathtt{threshold} * \max(\tau)
    /// \f]
    /// are considered; no ordering of the frequency array is needed.  
    Filter(const std::string& filename, double threshold=1e-3);
    /// Multiply transmission with a constant \p c.
    Filter& operator*= (double c);
    /// Divide transmission by a constant \p c.
    Filter& operator/= (double c);
    /// Multiply transmission with the one from \p f.
    Filter& operator*= (const Filter& f);
    /// Add filter curve from \p f.
    /// Currently only extending the transmission by concatenating \p f is
    /// implemented, overlapping filters are not.
    Filter& operator+= (const Filter& f);
    /// Get transmission at frequency \p nu.
    /// Linear interpolation across filter elements bracketing \p nu.
    double operator() (double nu) const;
    /// Save filter curve in the format required by constructor.
    void save(const std::string& filename) const;
    /// Erase consecutive elemenst with tranmission below \p threshold.
    void removeZeros(double threshold);
    /// Get minimum frequency with non-zero transmission.
    double getNuMin() const;
    /// Get maximum frequency with non-zero transmission.
    double getNuMax() const;
    /// Compute the integral \f$\int d\nu\ \tau/\nu\f$.
    /// Uses trapezoid rule.
    double computeNorm() const;
    /// Compute the average transmission between \p numin and \p numax.
    /// Uses linear interpolation and trapezoid rule.
    double avg(double numin, double numax) const;
    /// Compute effective wavelength (in Angstrom).
    double computeLambdaEff() const;
    /// Compute width of filter (in Angstrom).
    double computeWidth() const;
  protected:
    double prefactor, z;
  };

} // end namespace

#endif
