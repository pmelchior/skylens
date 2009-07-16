#include <skylens/Observation.h>
#include <skylens/Layer.h>
#include <skylens/Conventions.h>
#include <skylens/Conversion.h>
#include <shapelens/utils/IO.h>
#include <tclap/CmdLine.h>

using namespace skylens;
using namespace shapelens;

int main(int argc, char* argv[]) {
  // parse commandline
  TCLAP::CmdLine cmd("Test observation module of SkyLens++", ' ', "0.1");
  TCLAP::ValueArg<std::string> telname("t","telescope","Name of telescope",true,"","string", cmd);
  TCLAP::ValueArg<std::string> bandname("b","band","Name of filter band",true,"","string", cmd);
  TCLAP::ValueArg<double> expt("e","exposure_time","Exposure time (in seconds)",true,0,"double", cmd);
  TCLAP::ValueArg<std::string> skyspectrum("s","sky_spectrum","SED of the sky",true,"","string");
  TCLAP::ValueArg<double> sky_mag("S","sky_mag","Exposure time (in seconds)",true,0,"double");
  cmd.xorAdd(skyspectrum,sky_mag);
  TCLAP::ValueArg<double> airmass("m","air_mass","Airmass",false,1,"double",cmd);
  TCLAP::ValueArg<std::string> atm("a","atmosphere_spectrum","Atmospheric absorption spectrum",false,"","string");
  TCLAP::ValueArg<double> absorption("A","air_absorption","Average atmospheric absorption",false,0,"double");
  cmd.xorAdd(atm,absorption);
  cmd.parse(argc,argv);

  Telescope tel(telname.getValue(),bandname.getValue());
  double exptime = expt.getValue();
  Observation obs(tel,exptime);

  if (skyspectrum.isSet())
    obs.createSkyFluxLayer(sed(skyspectrum.getValue(),datapath));
  else
    obs.createSkyFluxLayer(sky_mag.getValue());
  if (atm.isSet())
    obs.computeTransmittance(filter(atm.getValue(),datapath),airmass.getValue());
  else if (absorption.isSet())
    obs.computeTransmittance(absorption.getValue(),airmass.getValue());
  
  const filter& transmittance = obs.getTotalTransmittance();

  std::cout << "total transmittance = " << transmittance.getQe() << std::endl;
  LayerStack& ls = SingleLayerStack::getInstance();
  double ADU_sky_pixel; // sky ADUs per pixel
  for (LayerStack::iterator iter = ls.begin(); iter != ls.end(); iter++) {
    if (iter->second->getType() == "SS")
      ADU_sky_pixel = iter->second->getFlux(Point<double>(0,0));
  }
  std::cout << "zeropoint = " << Conversion::zeroPoint(tel,transmittance,1) << std::endl;
  std::cout << "sky ADUs = " << ADU_sky_pixel << std::endl;
  std::cout << "sky magnitude = " << Conversion::flux2mag(Conversion::photons2flux(Conversion::ADU2photons(ADU_sky_pixel/(tel.pixsize*tel.pixsize),tel.gain),exptime,tel,transmittance)) << std::endl;
  
  return 0;
}
