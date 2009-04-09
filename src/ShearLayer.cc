#include <Layer.h>

using namespace skylens;

ShearLayer::ShearLayer(double z, complex<double> gamma) :
  gamma(gamma),
  // automatically creates a single instance of LayerStack
  ls(SingleLayerStack::getInstance())
{
  Layer::z = z;
  // try to insert this layer
  // fails if layer at that redshift already exists
  std::pair<LayerStack::iterator, bool> p = ls.insert(std::pair<double,Layer*>(z,this));
  if (p.second == true) // insertation took place
    me = p.first;
  else {
    std::ostringstream s;
    s << z;
    throw std::invalid_argument("ShearLayer: A layer at redshift " + s.str() + " already exists!");
  }
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
