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

  // connect to application DB
  SQLiteDB db;
  db.connect(outfileroot+".db");

  // get datapath
  std::string datapath = getDatapath();
  
  // set telescope and filter
  Telescope tel(boost::get<std::string>(config["TELESCOPE"]),
		boost::get<std::string>(config["FILTER"]));
  int exptime = boost::get<int>(config["EXPTIME"]);
  Observation obs(tel,exptime);

  // set (top-right corner of ) the global FoV
  // default = telescope's FoV
  Point<double> fov(tel.fov_x,tel.fov_y);
  try {
    fov(0) = boost::get<double>(config["GLOBAL_FOV_X"]);
    fov(1) = boost::get<double>(config["GLOBAL_FOV_Y"]);
  } catch (std::invalid_argument) {}

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
  for (int i=0; i< sourcefiles.size(); i++) {
    test_open(ifs,datapath,sourcefiles[i]);
    SourceCatalog sourcecat(sourcefiles[i]);
    // account for change of FoV from reference to telescope/global
    sourcecat.adjustNumber(fov);
    std::cout << "Sources:\t\t" << sourcecat.size() << std::endl;
    std::cout << "Replication ratio:\t" << sourcecat.getReplicationRatio() << std::endl;
    // place them randomly in the FoV 
    // and on the available redshifts of GalaxyLayers
    sourcecat.distribute(fov);
    // find reference bands with overlap to total transmittance
    sourcecat.selectOverlapBands(transmittance);
    // compute flux of source in each of the remaining bands
    sourcecat.computeADUinBands(tel,transmittance);
    // save source catalogs
    sourcecat.save(db,i);
  }

  t1 = time(NULL);
  std::cout << "Computation time: " << t1-t0 << " seconds" << std::endl;
}
