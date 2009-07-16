#include "../include/Layer.h"

using namespace skylens;

ClusterMemberLayer::ClusterMemberLayer(double z) :
  // automatically creates a single instance of LayerStack
  ls(SingleLayerStack::getInstance())
{
  Layer::z = z;
  Layer::transparent = false;
  ls.insert(std::pair<double,Layer*>(z,this));
}

// check for object at given position and return its flux
double ClusterMemberLayer::getFlux(const shapelens::Point<double>& P) const {
  // FIXME: howto use r-tree, how to specify galaxy ensemble
  double flux = 0;
  if (!transparent) {
  }
  return flux;
}

std::string ClusterMemberLayer::getType() const {
  return "SC";
}
