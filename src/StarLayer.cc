#include <Layer.h>

using namespace skylens;

StarLayer::StarLayer(double z, const PSF& psf) :
  // automatically creates a single instance of LayerStack
  ls(SingleLayerStack::getInstance()),
  psf(psf)
{
  Layer::z = z;
  // try to insert this layer
  // fails if layer at that redshift already exists
  std::pair<LayerStack::iterator, bool> p = ls.insert(std::pair<double,Layer*>(z,this));
  if (p.second == false) { // insertation failed
    std::ostringstream s;
    s << z;
    throw std::invalid_argument("StarLayer: A layer at redshift " + s.str() + " already exists!");
  }
}

// check for object at given position and return its flux
double StarLayer::getFlux(double x, double y) const {
  // FIXME: howto use r-tree, how to specify stellar ensemble
  double flux = 0;
  return flux;
}

std::string StarLayer::getType() const {
  return "S*";
}