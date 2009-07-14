#include "../include/Layer.h"

using namespace skylens;

DitherLayer::DitherLayer(double dx, double dy) :
  dx(dx), dy(dy),
  // automatically creates a single instance of LayerStack
  ls(SingleLayerStack::getInstance())
{
  Layer::z = -3;
  Layer::transparent = false;
  me = ls.insert(std::pair<double,Layer*>(z,this));
}

// sum all fluxes until the next transformation layer is found
double DitherLayer::getFlux(const shapelens::Point<double>& P) const {
  double flux = 0;
  shapelens::Point<double> P_ = P;
  if (!transparent) {
    P_(0) -= dx;
    P_(1) -= dy;
  }
  LayerStack::iterator iter = me;
  iter++; // next layer
  for (iter; iter != ls.end(); iter++) {
    std::string type = iter->second->getType();
    flux += iter->second->getFlux(P_);
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
  // it might be necessary to recompute some coordinates
  // or source layers contents if the displacement is to large...
}
