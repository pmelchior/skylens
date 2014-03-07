#ifndef SKYLENS_RNG_H
#define SKYLENS_RNG_H

#include <gsl/gsl_randist.h>
#include "Singleton.h"

namespace skylens {
  /// Common class for Random Number Generators.
  /// The main purpose of the class is to provide a single RNG with a single 
  /// seed such that a simulation can be exactly reproduced once the 
  /// seed is known.\n\n
  /// The preferred use is via a Singleton:
  /// \code
  /// RNG& rng = Singleton<RNG>::getInstance();
  /// const gsl_rng * r = rng.getRNG();
  /// \endcode
  class RNG {
  public:
    /// Constructor.
    RNG();
    /// Destructor
    ~RNG();
    /// Access to underlying RNG from GSL
    const gsl_rng * getRNG() const;
  private:
    gsl_rng * r;
  };
} // end namespace

#endif
