#include "../include/Layer.h"

using namespace skylens;

NullLayer::NullLayer():
  // automatically creates a single instance of LayerStack
  ls(SingleLayerStack::getInstance())
{
  Layer::z = -1000;
  Layer::transparent = false;
  me = ls.insert(std::pair<double,Layer*>(z,this));
}

// sum all fluxes until the next transformation layer is found
double NullLayer::getFlux(double x, double y) const {
  double flux = 0;
  LayerStack::iterator iter = me;
  iter++; // next layer
  for (iter; iter != ls.end(); iter++) {
    std::string type = iter->second->getType();
    flux += iter->second->getFlux(x,y);
    if (type[0] == 'T')
      break;
  }
  return flux;
}

std::string NullLayer::getType() const {
  return "T0";
}
