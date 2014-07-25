#include "../include/Layer.h"
#include <boost/foreach.hpp>
using namespace skylens;

StarLayer::StarLayer(const SourceModelList& stars_) :
  // automatically creates a single instance of LayerStack
  ls(SingleLayerStack::getInstance()),
  stars(stars_)
{
  Layer::z = 0;
  Layer::transparent = false;
  ls.insert(std::pair<double,Layer*>(z,this));

  // read out support rectangles from stars
  // and create RTree from it
  for (size_t i=0; i < stars.size(); i++) {
    shapelens::Rectangle<shapelens::data_t> support = (stars[i])->getSupport();
    BBoxIndex bbi(BBox(BPoint(support.ll[0], support.ll[1]), 
		       BPoint(support.tr[0], support.tr[1])),
		  i);
    rtree.insert(bbi);
  }
}

// check for object at given position and return its flux
double StarLayer::getFlux(const shapelens::Point<double>& P) const {
  double flux = 0;
  if (!transparent) {
    std::vector<BBoxIndex> result;
    BPoint bp(P(0), P(1));
    rtree.query(boost::geometry::index::intersects(bp), std::back_inserter(result));
    BOOST_FOREACH(BBoxIndex const& bbi, result)
      flux += std::max(0.,stars[bbi.second]->getValue(P));
  }
  return flux;
}

std::string StarLayer::getType() const {
  return "S*";
}
