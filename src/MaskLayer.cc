#include <Layer.h>

using namespace skylens;

MaskLayer::MaskLayer(double z, const std::list<shapelens::Polygon<double> >& masks) :
  masks(masks),
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
    throw std::invalid_argument("MaskLayer: A layer at redshift " + s.str() + " already exists!");
  }
}

// sum all fluxes until the next transformation layer is found
double MaskLayer::getFlux(double x, double y) const {
  double flux = 0;
  // check whether this position is masked
  bool masked = false;
  shapelens::Point2D<double> p(x,y);
  for (std::list<shapelens::Polygon<double> >::const_iterator iter = masks.begin(); iter != masks.end(); iter++) {
    if (iter->isInside(p)) {
      masked = true;
      break;
    }
  }
  if (masked == false)
    return 1;
  else return 0;
//   if (masked == false) {
//     LayerStack::iterator iter = me;
//     iter++; // next layer
//     for (iter; iter != ls.end(); iter++) {
//       std::string type = iter->second->getType();
//       flux += iter->second->getFlux(x,y);
//       if (type[0] == 'T')
// 	break;
//     }
//   }
//   return flux;
}

std::string MaskLayer::getType() const {
  return "TM";
}
