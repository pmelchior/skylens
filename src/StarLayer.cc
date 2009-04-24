#include <Layer.h>

using namespace skylens;

StarLayer::StarLayer(const PSF& psf) :
  // automatically creates a single instance of LayerStack
  ls(SingleLayerStack::getInstance()),
  psf(psf)
{
  Layer::z = 0;
  Layer::transparent = false;
  ls.insert(std::pair<double,Layer*>(z,this));
}

// check for object at given position and return its flux
double StarLayer::getFlux(double x, double y) const {
  // FIXME: howto use r-tree, how to specify stellar ensemble
  double flux = 0;
  if (!transparent) {
  }
  return flux;
}

std::string StarLayer::getType() const {
  return "S*";
}
