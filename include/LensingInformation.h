#ifndef SKYLENS_LENSINGINFORMATION_H
#define SKYLENS_LENSINGINFORMATION_H

#include <shapelens/utils/Singleton.h>
#include <map>

namespace skylens {
  /// Stores global lensing information.
  /// Access through a Singleton:
  /// \code
  /// LensingInformation& li = SingleLensingInformation::getInstance();
  /// \endcode
  class LensingInformation {
  public:
    /// Constructor.
    LensingInformation();
    /// Redshift of the first lensing layer.
    double z_first_lens;
    /// Angular diameter distance to source layers
    std::map<double, double> Ds;
    /// Iterator pointing on the entry of \p Ds of the current source layer.
    std::map<double, double>::iterator current_source;
    /// Angular diameter distance from lens to source layer.
    std::map<double, std::map<double, double> > Dls;
    /// Value of \f$c/H_0\f$ in Mpc (required to rescale Ds).
    double c_H0;
  };
  typedef shapelens::Singleton<LensingInformation> SingleLensingInformation;
} // end namespace

#endif
