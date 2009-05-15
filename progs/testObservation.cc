#include <skylens/Observation.h>
#include <skylens/Layer.h>
#include <skylens/Conventions.h>
#include <skylens/Conversion.h>
#include <shapelens/utils/IO.h>
using namespace skylens;
using namespace shapelens;

int main() {
  SUBARU subaru("I");
  sed sky("sky/moon1.fits",datapath);
  double time = 1000;
  Observation obs(subaru,sky,time);
  LayerStack& ls = SingleLayerStack::getInstance();
  double ADU_sky_pixel; // sky ADUs per pixel
  for (LayerStack::iterator iter = ls.begin(); iter != ls.end(); iter++) {
    if (iter->second->getType() == "SS")
      ADU_sky_pixel = iter->second->getFlux(0,0);
  }
  std::cout << "zeropoint = " << Conversion::zeroPoint(1,subaru) << std::endl;
  std::cout << "sky ADUs = " << ADU_sky_pixel << std::endl;
  std::cout << "sky magnitude 1 = " << Conversion::flux2mag(Conversion::photons2flux(Conversion::ADU2photons(ADU_sky_pixel/(subaru.pixsize*subaru.pixsize),subaru.gain),time,subaru)) << std::endl;
  std::cout << "sky magnitude 2 = " << Conversion::flux2mag(Conversion::ADU2flux(ADU_sky_pixel/(subaru.pixsize*subaru.pixsize),Conversion::zeroPoint(time,subaru))) << std::endl;
  return 0;
}
