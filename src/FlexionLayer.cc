#include "../include/Layer.h"

using namespace skylens;

FlexionLayer::FlexionLayer(double z, complex<double> F, complex<double> G) :
  // automatically creates a single instance of LayerStack
  ls(SingleLayerStack::getInstance())
{
  Layer::z = z;
  Layer::transparent = false;
  me = ls.insert(std::pair<double,Layer*>(z,this));

  // Bacon et al. (2006), eq. (14)
  double dGamma11 = 0.5*(real(F) + real(G));
  double dGamma12 = 0.5*(imag(G) - imag(F));
  double dGamma21 = 0.5*(imag(F) + imag(G));
  double dGamma22 = 0.5*(real(F) - real(G));
  // 1/2 of eq. (5)
  D111 = -dGamma11-0.5*dGamma22;
  D121 = -0.5*dGamma21;
  D211 = -0.5*dGamma21;
  D221 = -0.5*dGamma22;
  D112 = -0.5*dGamma21;
  D122 = -0.5*dGamma22;
  D212 = -0.5*dGamma22;
  D222 = dGamma12-0.5*dGamma21;
}

// sum all fluxes until the next transformation layer is found
double FlexionLayer::getFlux(const shapelens::Point<double>& P) const {
  double flux = 0, x_, y_;
  // apply eq. (3) without shear
  if (!transparent) {
    x_ = D111*P(0)*P(0) + (D112 + D121)*P(0)*P(1) + D122*P(1)*P(1);
    y_ = D211*P(0)*P(0) + (D212 + D221)*P(0)*P(1) + D222*P(1)*P(1);
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

std::string FlexionLayer::getType() const {
  return "TF";
}
