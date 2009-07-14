#include "../include/Layer.h"

namespace skylens {

Layer::~Layer() {}

double Layer::getRedshift() const {
  return z;
}

LayerStack::LayerStack() : std::multimap<double,Layer*>() {
}

LayerStack::~LayerStack() {
  for (LayerStack::iterator iter = this->begin(); iter != this->end(); iter++) {
    delete iter->second;
  }
}

std::ostream& operator<<(std::ostream& os, const LayerStack& ls) {
  for (LayerStack::const_iterator iter = ls.begin(); iter != ls.end(); iter++)
    os << iter->first << " -> " << iter->second->getType() << std::endl;
  return os;
}

}
