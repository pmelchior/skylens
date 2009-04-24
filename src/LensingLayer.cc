#include <Layer.h>
#include <utils/Interpolation.h>
#include <utils/IO.h>

using namespace skylens;

LensingLayer::LensingLayer(double z, std::string angle_file) :
  // automatically creates a single instance of LayerStack
  ls(SingleLayerStack::getInstance())
{
  Layer::z = z;
  Layer::transparent = false;
  me = ls.insert(std::pair<double,Layer*>(z,this));
  // FIXME: needs to open two real-valued images of first and second component
  //a = shapelens::Image<complex<double> >(angle_file);
}

// sum all fluxes until the next transformation layer is found
double LensingLayer::getFlux(double x, double y) const {
  double flux = 0;
  // apply lens equation
  complex<double> p(x,y);
  if (!transparent)
    p -= a.interpolate(x,y); // FIXME: what type of interpolation
  LayerStack::iterator iter = me;
  iter++; // next layer
  for (iter; iter != ls.end(); iter++) {
    std::string type = iter->second->getType();
    flux += iter->second->getFlux(real(p),imag(p));
    if (type[0] == 'T')
      break;
  }
  return flux;
}

std::string LensingLayer::getType() const {
  return "TL";
}
