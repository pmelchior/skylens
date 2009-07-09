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
    patches.push_back((*iter)->getSupport());
  rtree.insertNodes(patches);
}

// check for object at given position and return its flux
double StarLayer::getFlux(double x, double y) const {
  double flux = 0;
  if (!transparent) {
    shapelens::Point<double> p(x,y);
    std::list<unsigned long> l = rtree.getMatches(p);
    for(std::list<unsigned long>::const_iterator iter = l.begin(); iter != l.end(); iter++)
      flux += stars[*iter]->getValue(p);
  }
  return flux;
}

std::string StarLayer::getType() const {
  return "S*";
}
