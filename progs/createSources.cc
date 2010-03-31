#include <shapelens/ShapeLens.h>
#include <skylens/SkyLens.h>
#include <time.h>
#include <tclap/CmdLine.h>
#include <fstream>

using namespace skylens;
using namespace shapelens;

int main(int argc, char* argv[]) {
  // parse commandline
  TCLAP::CmdLine cmd("Create sources for  SkyLens++ simulator", ' ', "0.1");
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

  // set outfile: no catch, this must be set
  std::string outfile = boost::get<std::string>(config["OUTFILE"]);
  std::string outfileroot = outfile.substr(0,outfile.rfind('.'));

  // get datapath
  std::string datapath = getDatapath();
  
  // set telescope and filter
  Telescope tel(boost::get<std::string>(config["TELESCOPE"]),
		boost::get<std::string>(config["FILTER"]));
  int exptime = boost::get<int>(config["EXPTIME"]);
  Observation obs(tel,exptime);

  // set global cosmology: default is vanilla CDM
  cosmology& cosmo = SingleCosmology::getInstance();
  try {
    std::string cosmofile = boost::get<std::string>(config["COSMOLOGY"]);
    // FIXME: implementation missing
  } catch (std::invalid_argument) {}

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

  // set global RNG seed if demanded
  try {
    int seed =  boost::get<int>(config["RNG_SEED"]);
    RNG& rng = Singleton<RNG>::getInstance();
    const gsl_rng * r = rng.getRNG();
    gsl_rng_set(r,seed);
  } catch (std::invalid_argument) {}

  // get sources from config files
  std::vector<std::string> sourcefiles = boost::get<std::vector<std::string> >(config["SOURCES"]);
  std::ostringstream fileext;
  for (int i=0; i< sourcefiles.size(); i++) {
    test_open(ifs,datapath,sourcefiles[i]);
    SourceCatalog sourcecat(sourcefiles[i]);
    // account for change of FoV from reference to telescope
    sourcecat.adjustNumber(tel);
    std::cout << "Sources: " << sourcecat.size() << std::endl;
    // place them randomly in the FoV 
    // and on the available redshifts of GalaxyLayers
    sourcecat.distribute(tel);
    // find reference bands with overlap to tel
    sourcecat.selectOverlapBands(tel);
    // compute flux of source in each of the remaining bands
    sourcecat.computeADUinBands(tel);
    // save source catalogs
    if (sourcefiles.size() > 1) { // multiple source catalogs
      fileext.str("");
      fileext << outfileroot+".sourcecat" << i+1;
      sourcecat.save(fileext.str());
    } else
      sourcecat.save(outfileroot+".sourcecat");
  }

  t1 = time(NULL);
  std::cout << "Computation time: " << t1-t0 << " seconds" << std::endl;
}
