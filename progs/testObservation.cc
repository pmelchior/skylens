#include <skylens/SkyLens.h>
#include <shapelens/utils/IO.h>
#include <tclap/CmdLine.h>
#include <fenv.h>

using namespace skylens;
using namespace shapelens;

int main(int argc, char* argv[]) {
  // parse commandline
  TCLAP::CmdLine cmd("Compute fluxes for  SkyLens++ simulator", ' ', "0.2");
  TCLAP::ValueArg<std::string> configfile("c","config","Configuration file",true,"","string", cmd);
  cmd.parse(argc,argv);
  
  // for measuring computation time
  time_t t0,t1;
  t0 = time(NULL);

  // read in global config file
  Property config;
  std::ifstream ifs(configfile.getValue().c_str());
  config.read(ifs);
  ifs.close();

  // get datapath
  std::string datapath = getDatapath();
  
  // set telescope and filter
  Telescope tel(boost::get<std::string>(config["TELESCOPE"]),
		boost::get<std::string>(config["FILTER"]));
  int exptime = boost::get<int>(config["EXPTIME"]);
  Observation obs(tel,exptime);

  // set absorption of the atmosphere
  double airmass = boost::get<data_t>(config["AIRMASS"]);
  try {
    std::string absorption = boost::get<std::string>(config["ATMOSPHERE"]);
    test_open(ifs,datapath,absorption);
    obs.computeTransmittance(absorption,airmass);
  } catch (boost::bad_get) {
    double  absorption = boost::get<data_t>(config["ATMOSPHERE"]);
    obs.computeTransmittance(absorption,airmass);
  }
  const filter& transmittance = obs.getTotalTransmittance();

  // set emission of the sky
  try {
    std::string sky = boost::get<std::string>(config["SKY"]);
    test_open(ifs,datapath,sky);
    obs.createSkyFluxLayer(sed(sky,"/"));
   } catch (boost::bad_get) {
    double sky = boost::get<double>(config["SKY"]);
    obs.createSkyFluxLayer(sky);
  }

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
