#include <Layer.h>
#include <iostream>

using namespace skylens;
int main() {
  LensingLayer ll1(0.5,"");
  LensingLayer ll2(1,"");
  GalaxyLayer lg1(0.75);
  GalaxyLayer lg2(2);
  DitherLayer ld(0,0.1,0.2);
  PSF psf("data/SUBARU/psf.fits");
  StarLayer lS(1e-3,psf);
  LayerStack& ls = SingleLayerStack::getInstance();
  for (LayerStack::iterator iter = ls.begin(); iter != ls.end(); iter++) {
    std::cout << iter->second->getRedshift() << "\t" << iter->second->getType() << std::endl;
  }
  
  std::cout << "total flux at (0,0) = " << ls.begin()->second->getFlux(0,0) << std::endl;
}
