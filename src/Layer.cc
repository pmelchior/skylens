#include <Layer.h>
#include <stdexcept>
#include <sstream>

using namespace skylens;

Layer::~Layer() {}

double Layer::getRedshift() const {
  return z;
}
