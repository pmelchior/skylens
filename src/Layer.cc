#include "../include/Layer.h"

using namespace skylens;

Layer::~Layer() {}

double Layer::getRedshift() const {
  return z;
}

LayerStack::LayerStack() : std::multimap<double,Layer*>() {
}

LayerStack::~LayerStack() {
  for (LayerStack::iterator iter = this->begin(); iter != this->end(); iter++)
    delete iter->second;
}
