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
    /// Angular diameter distance to all source layers
    std::map<double, double> Ds;
    /// Iterator pointing on the entry of \p Ds of the current source layer.
    std::map<double, double>::iterator current_source;
  };
  typedef shapelens::Singleton<LensingInformation> SingleLensingInformation;
} // end namespace

#endif
