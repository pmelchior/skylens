#include <Layer.h>
#include <stdexcept>
#include <sstream>

using namespace skylens;

Layer::~Layer() {}

double Layer::getRedshift() const {
  return z;
}

LensingLayer::LensingLayer(double z, std::string angle_file) :
  // opens deflection angle map
  a(shapelens::Image<double>(angle_file)),
  // automatically creates a single instance of LayerStack
  ls(SingleLayerStack::getInstance())
{
  Layer::z = z;
  // try to insert this layer
  // fails if layer at that redshift already exists
  std::pair<LayerStack::iterator, bool> p = ls.insert(std::pair<double,Layer*>(z,this));
  if (p.second == true) // insertation took place
    me = p.first;
  else {
    std::ostringstream s;
    s << z;
    throw std::invalid_argument("LensingLayer: A layer at redshift " + s.str() + " already exists!");
  }
}

// sum all fluxes until the next transformation layer is found
double LensingLayer::getFlux(double x, double y) const {
  std::cout << getType() << " called at (" << x << "," << y << ")" << std::endl;
  double flux = 0;
  LayerStack::iterator iter = me;
  iter++; // next layer
  for (iter; iter != ls.end(); iter++) {
    std::string type = iter->second->getType();
    double x_, y_;
    // FIXME:
    x_ = x + 0.1;
    y_ = y - 0.2;
    // the following assumes deflection angle map to be of complex type
    //x_ = x - real(a(x,y));
    //y_ = y - imag(a(x,y));
    flux += iter->second->getFlux(x_,y_);
    if (type[0] == 'T')
      break;
  }
  return flux;
}

std::string LensingLayer::getType() const {
  return "TL";
}

DitherLayer::DitherLayer(double z, double dx, double dy) :
  dx(dx), dy(dy),
  // automatically creates a single instance of LayerStack
  ls(SingleLayerStack::getInstance())
{
  Layer::z = z;
  // try to insert this layer
  // fails if layer at that redshift already exists
  std::pair<LayerStack::iterator, bool> p = ls.insert(std::pair<double,Layer*>(z,this));
  if (p.second == true) // insertation took place
    me = p.first;
  else {
    std::ostringstream s;
    s << z;
    throw std::invalid_argument("DitherLayer: A layer at redshift " + s.str() + " already exists!");
  }
}

// sum all fluxes until the next transformation layer is found
double DitherLayer::getFlux(double x, double y) const {
  std::cout << getType() << " called at (" << x << "," << y << ")" << std::endl;
  double flux = 0;
  LayerStack::iterator iter = me;
  iter++; // next layer
  for (iter; iter != ls.end(); iter++) {
    std::string type = iter->second->getType();
    double x_, y_;
    x_ = x - dx;
    y_ = y - dy;
    flux += iter->second->getFlux(x_,y_);
    if (type[0] == 'T')
      break;
  }
  return flux;
}

std::string DitherLayer::getType() const {
  return "TD";
}

GalaxyLayer::GalaxyLayer(double z) :
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
}

// check for object at given position and return its flux
double GalaxyLayer::getFlux(double x, double y) const {
  std::cout << getType() << " called at (" << x << "," << y << ")" << std::endl;
  // FIXME: howto use r-tree, how to specify galaxy ensemble
  double flux = 2;
  return flux;
}

std::string GalaxyLayer::getType() const {
  return "SG";
}

ClusterMemberLayer::ClusterMemberLayer(double z) :
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
    throw std::invalid_argument("ClusterMemberLayer: A layer at redshift " + s.str() + " already exists!");
  }
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


StarLayer::StarLayer(double z, const PSF& psf) :
  // automatically creates a single instance of LayerStack
  ls(SingleLayerStack::getInstance()),
  psf(psf)
{
  Layer::z = z;
  // try to insert this layer
  // fails if layer at that redshift already exists
  std::pair<LayerStack::iterator, bool> p = ls.insert(std::pair<double,Layer*>(z,this));
  if (p.second == false) { // insertation failed
    std::ostringstream s;
    s << z;
    throw std::invalid_argument("StarLayer: A layer at redshift " + s.str() + " already exists!");
  }
}

// check for object at given position and return its flux
double StarLayer::getFlux(double x, double y) const {
  std::cout << getType() << " called at (" << x << "," << y << ")" << std::endl;
  // FIXME: howto use r-tree, how to specify stellar ensemble
  double flux = 0.5;
  return flux;
}

std::string StarLayer::getType() const {
  return "S*";
}
