#include "../include/Layer.h"

using namespace skylens;

StarLayer::StarLayer(const shapelens::SourceModelList& stars) :
  // automatically creates a single instance of LayerStack
  ls(SingleLayerStack::getInstance()),
  stars(stars)
{
  Layer::z = 0;
  Layer::transparent = false;
  ls.insert(std::pair<double,Layer*>(z,this));

  // read out support rectangles from stars
  // and create RTree from it
  std::vector<shapelens::Rectangle<double> > patches;
  for (shapelens::SourceModelList::const_iterator iter = stars.begin(); iter != stars.end(); iter++)
    patches.push_back((*iter)->getSupport().getBoundingBox());
  rtree.insertNodes(patches);
}

// check for object at given position and return its flux
double StarLayer::getFlux(const shapelens::Point<double>& P) const {
  double flux = 0;
  if (!transparent) {
    std::list<unsigned long> l = rtree.getMatches(P);
    for(std::list<unsigned long>::const_iterator iter = l.begin(); iter != l.end(); iter++)
      flux += stars[*iter]->getValue(P);
  }
  return flux;
}

std::string StarLayer::getType() const {
  return "S*";
}
