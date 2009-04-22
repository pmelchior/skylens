#include <Layer.h>

using namespace skylens;

ShearLayer::ShearLayer(double z, complex<double> gamma) :
  gamma(gamma),
  // automatically creates a single instance of LayerStack
  ls(SingleLayerStack::getInstance())
{
  Layer::z = z;
  me = ls.insert(std::pair<double,Layer*>(z,this));
}

// sum all fluxes until the next transformation layer is found
double ShearLayer::getFlux(double x, double y) const {
  double flux = 0;
  // apply lens equation
  double x_ = (1-real(gamma))*x - imag(gamma)*y;
  double y_ = -imag(gamma)*x + (1+real(gamma))*y;
  LayerStack::iterator iter = me;
  iter++; // next layer
  for (iter; iter != ls.end(); iter++) {
    std::string type = iter->second->getType();
    flux += iter->second->getFlux(x_,y_);
    if (type[0] == 'T')
      break;
  }
  return flux;
}

std::string ShearLayer::getType() const {
  return "TS";
}
