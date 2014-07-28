#include "../include/Layer.h"

using namespace skylens;

SkyFluxLayer::SkyFluxLayer(double flux) :
  // automatically creates a single instance of LayerStack
  ls(SingleLayerStack::getInstance()),
  flux(flux)
{
  Layer::z = -1;
  Layer::transparent = false;
  ls.insert(std::pair<double,Layer*>(z,this));
}

double SkyFluxLayer::getFlux(const shapelens::Point<double>& P, double* z_) const {
  if (z_ != NULL)
    if (*z_ != z)
      return 0;
  
  if (!transparent)
    return flux;
  else
    return 0;
}

std::string SkyFluxLayer::getType() const {
  return "SS";
}

void SkyFluxLayer::setFlux(double flux_) {
  flux = flux_;
}
