#include "../include/Layer.h"

using namespace skylens;

ShearLayer::ShearLayer(double z, std::complex<double> gamma) :
  gamma(gamma),
  // automatically creates a single instance of LayerStack
  ls(SingleLayerStack::getInstance())
{
  Layer::z = z;
  Layer::transparent = false;
  me = ls.insert(std::pair<double,Layer*>(z,this));
}

// sum all fluxes until the next transformation layer is found
double ShearLayer::getFlux(const shapelens::Point<double>& P) const {
  double flux = 0, x_, y_;
  // apply lens equation
  if (!transparent) {
    x_ = (1-real(gamma))*P(0) - imag(gamma)*P(1);
    y_ = -imag(gamma)*P(0) + (1+real(gamma))*P(1);
  } else {
    x_ = P(0);
    y_ = P(1);
  }
  LayerStack::iterator iter = me;
  iter++; // next layer
  for (iter; iter != ls.end(); iter++) {
    std::string type = iter->second->getType();
    flux += iter->second->getFlux(shapelens::Point<double>(x_,y_));
    if (type[0] == 'T')
      break;
  }
  return flux;
}

std::string ShearLayer::getType() const {
  return "TS";
}
