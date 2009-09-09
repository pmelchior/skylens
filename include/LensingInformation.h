#ifndef SKYLENS_LENSINGINFORMATION_H
#define SKYLENS_LENSINGINFORMATION_H

#include <shapelens/utils/Singleton.h>

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
    /// Redshift of the first lensed source layer.
    double z_first_lensed_source;
    /// Redshift of the currently active source layer.
    double z_current_source;
  };
  typedef shapelens::Singleton<LensingInformation> SingleLensingInformation;
} // end namespace

#endif
