#include <Layer.h>

using namespace skylens;

ConstFluxLayer::ConstFluxLayer(double z, double flux) :
  // automatically creates a single instance of LayerStack
  ls(SingleLayerStack::getInstance()),
  flux(flux)
{
  Layer::z = z;
  // try to insert this layer
  // fails if layer at that redshift already exists
  std::pair<LayerStack::iterator, bool> p = ls.insert(std::pair<double,Layer*>(z,this));
  if (p.second == false) { // insertation failed
    std::ostringstream s;
    s << z;
    throw std::invalid_argument("ConstFluxLayer: A layer at redshift " + s.str() + " already exists!");
  }
}

double ConstFluxLayer::getFlux(double x, double y) const {
  return flux;
}

std::string ConstFluxLayer::getType() const {
  return "S=";
}
