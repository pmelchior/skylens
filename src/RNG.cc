#include "../include/RNG.h"

namespace skylens {
  RNG::RNG() {
    const gsl_rng_type * T;
    // uses env varibale to set up RNG
    // see GSL manual
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
  }

  RNG::~RNG() {
    gsl_rng_free (r);
  }

  const gsl_rng * RNG::getRNG() const {
    return r;
  }
} // end namespace
