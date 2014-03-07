#include "../include/Layer.h"

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
  std::vector<shapelens::Rectangle<double> > patches;
  for (SourceModelList::const_iterator iter = galaxies.begin(); iter != galaxies.end(); iter++)
    patches.push_back((*iter)->getSupport());
  rtree.insertNodes(patches);
}

// check for object at given position and return its flux
double GalaxyLayer::getFlux(const shapelens::Point<double>& P) const {
  double flux = 0;
  if (!transparent) {
    std::list<unsigned long> l = rtree.getMatches(P);
    for(std::list<unsigned long>::const_iterator iter = l.begin(); iter != l.end(); iter++)
      flux += std::max(0.,galaxies[*iter]->getValue(P));
  }
  return flux;
}

std::string GalaxyLayer::getType() const {
  return "SG";
}
