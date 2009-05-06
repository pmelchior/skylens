#include "../include/Layer.h"
#include <gsl/gsl_randist.h>

using namespace skylens;

NoiseLayer::NoiseLayer(double ron, double flat_field) :
  ron(ron), flat_field(flat_field), 
  // automatically creates a single instance of LayerStack
  ls(SingleLayerStack::getInstance())
{
  Layer::z = -2;
  Layer::transparent = false;
  me = ls.insert(std::pair<double,Layer*>(z,this));

  // set up RNG
  const gsl_rng_type * T;
  gsl_rng_env_setup(); // read env variables to set seed and RNG type
  T = gsl_rng_default;
  r = gsl_rng_alloc(T);
}

NoiseLayer::~NoiseLayer() {
  gsl_rng_free (r);
}

// sum all fluxes until the next transformation layer is found
// and then add noise to it
double NoiseLayer::getFlux(double x, double y) const {
  double flux = 0;
  LayerStack::iterator iter = me;
  iter++; // next layer
  for (iter; iter != ls.end(); iter++) {
    std::string type = iter->second->getType();
    flux += iter->second->getFlux(x,y);
    if (type[0] == 'T')
      break;
  }
  if (!transparent)
    // eq. (31)
    flux += gsl_ran_gaussian_ziggurat (r,sqrt(ron + fabs(flux) + flat_field*flux*flux));
  return flux;
}

std::string NoiseLayer::getType() const {
  return "TN";
}
