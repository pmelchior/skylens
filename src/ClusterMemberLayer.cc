#include <Layer.h>

using namespace skylens;

ClusterMemberLayer::ClusterMemberLayer(double z) :
  // automatically creates a single instance of LayerStack
  ls(SingleLayerStack::getInstance())
{
  Layer::z = z;
  ls.insert(std::pair<double,Layer*>(z,this));
}

// check for object at given position and return its flux
double ClusterMemberLayer::getFlux(double x, double y) const {
  // FIXME: howto use r-tree, how to specify galaxy ensemble
  double flux = 0;
  return flux;
}

std::string ClusterMemberLayer::getType() const {
  return "SC";
}
