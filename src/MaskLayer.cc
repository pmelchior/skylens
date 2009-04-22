#include <Layer.h>

using namespace skylens;

MaskLayer::MaskLayer(double FoV, const std::list<shapelens::Polygon<double> >& masks_) :
  masks(masks_),
  // automatically creates a single instance of LayerStack
  ls(SingleLayerStack::getInstance())
{
  Layer::z = -4;
  me = ls.insert(std::pair<double,Layer*>(Layer::z,this));

  // rescale masks coordinates to arcsec via FoV
  for (std::list<shapelens::Polygon<double> >::iterator iter = masks.begin(); iter != masks.end(); iter++)
    iter->rescale(FoV);
}

MaskLayer::MaskLayer(double FoV, std::string maskfile) :
  // automatically creates a single instance of LayerStack
  ls(SingleLayerStack::getInstance())
{
  Layer::z = -4;
  me = ls.insert(std::pair<double,Layer*>(Layer::z,this));

  // open maskfile and read (possibly many) polygons
  std::ifstream ifs(maskfile.c_str());
  if (ifs.bad())
    throw std::ios_base::failure("MaskLayer: Cannot read from istream");

  shapelens::Polygon<double> poly;
  while (ifs >> poly) {
    poly.rescale(FoV);
    masks.push_back(poly);
  }
  ifs.close();
}

// sum all fluxes until the next transformation layer is found
double MaskLayer::getFlux(double x, double y) const {
  double flux = 0;
  // check whether this position is masked
  bool masked = false;
  shapelens::Point2D<double> p(x,y);
  for (std::list<shapelens::Polygon<double> >::const_iterator piter = masks.begin(); piter != masks.end(); piter++) {
    if (piter->isInside(p)) {
      masked = true;
      break;
    }
  }
  if (masked == false) {
    LayerStack::iterator iter = me;
    iter++; // next layer
    for (iter; iter != ls.end(); iter++) {
      std::string type = iter->second->getType();
      flux += iter->second->getFlux(x,y);
      if (type[0] == 'T')
	break;
    }
  }
  return flux;
}

std::string MaskLayer::getType() const {
  return "TM";
}
