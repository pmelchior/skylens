#include "../include/Layer.h"
#include <boost/foreach.hpp>

using namespace skylens;

GalaxyLayer::GalaxyLayer(double z, const SourceModelList& galaxies_) :
  galaxies(galaxies_),
  // automatically creates a single instance of LayerStack
  ls(SingleLayerStack::getInstance())
{
  Layer::z = z;
  Layer::transparent = false;
  ls.insert(std::pair<double,Layer*>(z,this));
  
  // read out support rectangles from galaxies
  // and create RTree from it
  for (size_t i=0; i < galaxies.size(); i++) {
    shapelens::Rectangle<shapelens::data_t> support = (galaxies[i])->getSupport();
    BBoxIndex bbi(BBox(BPoint(support.ll[0], support.ll[1]), 
		       BPoint(support.tr[0], support.tr[1])),
		  i);
    rtree.insert(bbi);
  }
}

// check for object at given position and return its flux
double GalaxyLayer::getFlux(const shapelens::Point<double>& P, double* z_) const {
  double flux = 0;
  if (z_ != NULL)
    if (*z_ != z)
      return flux;

  if (!transparent) {
    std::vector<BBoxIndex> result;
    BPoint bp(P(0), P(1));
    rtree.query(boost::geometry::index::intersects(bp), std::back_inserter(result));
    BOOST_FOREACH(BBoxIndex const& bbi, result)
      flux += std::max(0.,galaxies[bbi.second]->getValue(P));
  }
  return flux;
}

std::string GalaxyLayer::getType() const {
  return "SG";
}
