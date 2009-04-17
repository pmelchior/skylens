#include <Layer.h>

using namespace skylens;

GalaxyLayer::GalaxyLayer(double z, const shapelens::SourceModelList& galaxies) :
  galaxies(galaxies),
  // automatically creates a single instance of LayerStack
  ls(SingleLayerStack::getInstance())
{
  Layer::z = z;
  // try to insert this layer
  // fails if layer at that redshift already exists
  std::pair<LayerStack::iterator, bool> p = ls.insert(std::pair<double,Layer*>(z,this));
  if (p.second == false) { // insertation failed
    std::ostringstream s;
    s << z;
    throw std::invalid_argument("GalaxyLayer: A layer at redshift " + s.str() + " already exists!");
  }

  // read out support rectangles from galaxies
  // and create RTree from it
  std::vector<shapelens::Rectangle<double> > patches;
  for (shapelens::SourceModelList::const_iterator iter = galaxies.begin(); iter != galaxies.end(); iter++)
    patches.push_back((*iter)->getSupport());
  rtree.insertNodes(patches);
}

// check for object at given position and return its flux
double GalaxyLayer::getFlux(double x, double y) const {
  double flux = 0;
  shapelens::Point2D<double> p(x,y);
  std::list<unsigned long> l = rtree.getMatches(p);
  for(std::list<unsigned long>::const_iterator iter = l.begin(); iter != l.end(); iter++)
    flux += galaxies[*iter]->getValue(p);
  return flux;
}

std::string GalaxyLayer::getType() const {
  return "SG";
}
