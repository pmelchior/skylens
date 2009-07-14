#include <skylens/Observation.h>
#include <skylens/Layer.h>
#include <skylens/Conventions.h>
#include <skylens/Conversion.h>
#include <shapelens/utils/IO.h>
using namespace skylens;
using namespace shapelens;

int main() {
  Telescope tel("SUBARU","B");
  sed sky("sky/moon1.fits",datapath);
  filter atmosphere ("sites/lasilla.fits",datapath);
  double time = 1000;
  Observation obs(tel,time,sky,atmosphere,1,1);
  const filter& transmittance = obs.getTotalTransmittance();
  LayerStack& ls = SingleLayerStack::getInstance();
  double ADU_sky_pixel; // sky ADUs per pixel
  for (LayerStack::iterator iter = ls.begin(); iter != ls.end(); iter++) {
    if (iter->second->getType() == "SS")
      ADU_sky_pixel = iter->second->getFlux(0,0);
  }
  std::cout << "zeropoint = " << Conversion::zeroPoint(tel,transmittance,1) << std::endl;
  std::cout << "sky ADUs = " << ADU_sky_pixel << std::endl;
  std::cout << "sky magnitude = " << Conversion::flux2mag(Conversion::photons2flux(Conversion::ADU2photons(ADU_sky_pixel/(tel.pixsize*tel.pixsize),tel.gain),time,tel,transmittance)) << std::endl;
  return 0;
}
