#include <Layer.h>

using namespace skylens;

DitherLayer::DitherLayer(double dx, double dy) :
  dx(dx), dy(dy),
  // automatically creates a single instance of LayerStack
  ls(SingleLayerStack::getInstance())
{
  Layer::z = -4;
  Layer::transparent = false;
  me = ls.insert(std::pair<double,Layer*>(z,this));
}

// sum all fluxes until the next transformation layer is found
double DitherLayer::getFlux(double x, double y) const {
  double flux = 0, x_, y_;
  if (!transparent) {
    x_ = x - dx;
    y_ = y - dy;
  } else {
    x_ = x;
    y_ = y;
  }
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

std::string DitherLayer::getType() const {
  return "TD";
}

void DitherLayer::setDisplacement(double dx_, double dy_) {
  dx = dx_;
  dy = dy_;
  // it might be necessary to recompute some coordindates
  // or source layers contents of the displacement is to large...
}
