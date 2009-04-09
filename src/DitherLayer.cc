#include <Layer.h>

using namespace skylens;

DitherLayer::DitherLayer(double z, double dx, double dy) :
  dx(dx), dy(dy),
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
    throw std::invalid_argument("DitherLayer: A layer at redshift " + s.str() + " already exists!");
  }
}

// sum all fluxes until the next transformation layer is found
double DitherLayer::getFlux(double x, double y) const {
  double flux = 0;
  LayerStack::iterator iter = me;
  iter++; // next layer
  for (iter; iter != ls.end(); iter++) {
    std::string type = iter->second->getType();
    double x_, y_;
    x_ = x - dx;
    y_ = y - dy;
    flux += iter->second->getFlux(x_,y_);
    if (type[0] == 'T')
      break;
  }
  return flux;
}

std::string DitherLayer::getType() const {
  return "TD";
}

void DitherLayer::setDisplacement(double dx_, double dy_) {
  dx = dx_;
  dy = dy_;
}
