#include "../include/Layer.h"

using namespace skylens;

MaskLayer::MaskLayer(const std::list<shapelens::Polygon<double> >& masks_) :
  masks(masks_),
  // automatically creates a single instance of LayerStack
  ls(SingleLayerStack::getInstance())
{
  Layer::z = -2;
  Layer::transparent = false;
  me = ls.insert(std::pair<double,Layer*>(Layer::z,this));
}

MaskLayer::MaskLayer(std::string maskfile) :
  // automatically creates a single instance of LayerStack
  ls(SingleLayerStack::getInstance())
{
  Layer::z = -3;
  Layer::transparent = false;
  me = ls.insert(std::pair<double,Layer*>(Layer::z,this));

  // open maskfile and read (possibly many) polygons
  std::ifstream ifs(maskfile.c_str());
  if (ifs.bad())
    throw std::ios_base::failure("MaskLayer: Cannot read from istream");

  shapelens::Polygon<double> poly;
  while (ifs >> poly) {
    masks.push_back(poly);
  }
  ifs.close();
}

// sum all fluxes until the next transformation layer is found
double MaskLayer::getFlux(const shapelens::Point<double>& P) const {
  double flux = 0;
  // check whether this position is masked
  bool masked = false;
  if (!transparent) {
    for (std::list<shapelens::Polygon<double> >::const_iterator piter = masks.begin(); piter != masks.end(); piter++) {
      if (piter->contains(P)) {
	masked = true;
	break;
      }
    }
  }
  if (masked == false) {
    LayerStack::iterator iter = me;
    iter++; // next layer
    for (iter; iter != ls.end(); iter++) {
      std::string type = iter->second->getType();
      flux += iter->second->getFlux(P);
      if (type[0] == 'T')
	break;
    }
  }
  return flux;
}

std::string MaskLayer::getType() const {
  return "TM";
}
