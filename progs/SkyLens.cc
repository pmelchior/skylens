#include <shapelens/ShapeLens.h>
#include <skylens/SkyLens.h>
#include <time.h>
#include <tclap/CmdLine.h>
#include <fstream>

using namespace skylens;
using namespace shapelens;

int main(int argc, char* argv[]) {
  // parse commandline
  TCLAP::CmdLine cmd("Run SkyLens++ simulator", ' ', "0.4");
  TCLAP::ValueArg<std::string> configfile("c","config","Configuration file",true,"","string", cmd);
  TCLAP::SwitchArg useSources("u","use_sources","Use precomputed sources",cmd, false);
  TCLAP::SwitchArg saveSources("s","save_sources","Save catalog of sources",cmd, false);
  TCLAP::SwitchArg onlySources("o","only_sources","Only prepare source catalog (sets -s)",cmd, false);
  cmd.parse(argc,argv);


  // for measuring computation time
  time_t t0,t1;
  t0 = time(NULL);
  std::cout << "# SkyLens++ v" << cmd.getVersion() << " (git" << STRINGIFY(GITREV) << ")" << std::endl;
  
  // read in global config file
  std::cout << "# reading global configuration from " << configfile.getValue() << std::endl;
  Property config;
  std::ifstream ifs(configfile.getValue().c_str());
  config.read(ifs);
  ifs.close();

  // set outfile: no catch, this must be set
  std::string outfile = boost::get<std::string>(config["OUTFILE"]);
  std::string fileroot = boost::get<std::string>(config["PROJECT"]);

  // connect to application DB
  SQLiteDB& db = ApplicationDB::getInstance();
  db.connect(fileroot+".db");

  // get datapath
  std::string datapath = getDatapath();
  
  // set telescope and filter
  std::cout << "# setting up Telescope and Obervation" << std::endl;
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
  Cosmology& cosmo = SingleCosmology::getInstance();
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
  const Filter& transmittance = obs.getTotalTransmittance();

  // set emission of the sky
  try {
    std::string sky = boost::get<std::string>(config["SKY"]);
    test_open(ifs,datapath,sky);
    obs.createSkyFluxLayer(SED(sky));
   } catch (boost::bad_get) {
    double sky = boost::get<double>(config["SKY"]);
    obs.createSkyFluxLayer(sky);
  }

  // set pixel noise: on or off
  try {
    int noise = boost::get<int>(config["NOISE"]);
    if (noise)
      obs.setNoise();
  } catch (std::invalid_argument) {}

  // set global RNG seed if demanded
  try {
    int seed =  boost::get<int>(config["RNG_SEED"]);
    RNG& rng = Singleton<RNG>::getInstance();
    const gsl_rng * r = rng.getRNG();
    gsl_rng_set(r,seed);
  } catch (std::invalid_argument) {}

  // get sources from config files
  std::vector<std::string> files = boost::get<std::vector<std::string> >(config["SOURCES"]);
  for (int i=0; i< files.size(); i++) {
    // compute source information from DB
    if (!useSources.isSet()) {
      std::cout << "# retrieving source information from " << files[i] << ":" << std::endl;
      test_open(ifs,datapath,files[i]);
      SourceCatalog sourcecat(files[i]);
      // account for change of FoV from reference to telescope/global
      sourcecat.adjustNumber(fov);
      std::cout << "#   number = " << sourcecat.size() << std::endl;
      std::cout << "#   replication ratio = " << sourcecat.getReplicationRatio() << std::endl;
      // place them randomly in the FoV 
      // and on the available redshifts of GalaxyLayers
      sourcecat.distribute(fov);
      // compute flux of source in each of the remaining bands
      std::cout << "# computing source ADU in overlapping bands" << std::endl;
      sourcecat.computeADUinBands(tel,transmittance);
      // save source catalogs, if demanded
      if (saveSources.isSet() || onlySources.isSet()) {
	std::cout << "# saving sources" << std::endl;
	// save source catalogs
	sourcecat.save(db,i);
      }
      // create GalaxyLayers from sources
      if (!onlySources.isSet())
	sourcecat.createGalaxyLayers(exptime);
    }
    else if (!onlySources.isSet()) { // use precomputed sources
      std::cout << "# creating source models from " << files[i] << std::endl;
      SourceCatalog sourcecat(db,i);
      sourcecat.createGalaxyLayers(exptime);
    }
  }

  // stop here: no ray-tracing required
  if (onlySources.isSet()) {
    t1 = time(NULL);
    std::cout << "# computation time: " << t1-t0 << " seconds" << std::endl;
    exit(0);
  }

  // read in lens config
  Point<double> center;
  try {
    files = boost::get<std::vector<std::string> > (config["LENSES"]);
    for (int i=0; i < files.size(); i++) {
      std::cout << "# setting up lens from " << files[i] << std::endl;
      test_open(ifs,datapath,files[i]);
      Property lensconfig;
      lensconfig.read(ifs);
      // create lens layer
      center(0) = boost::get<double>(lensconfig["POS_X"]);
      center(1) = boost::get<double>(lensconfig["POS_Y"]);
      std::string anglefile = boost::get<std::string>(lensconfig["ANGLEFILE"]);
      test_open(ifs,datapath,anglefile);
      new LensingLayer(boost::get<double>(lensconfig["REDSHIFT"]),
		       anglefile,
		       center);
    }
  } catch (std::invalid_argument) {}

  // FIXME: convolution config

  // FIXME: star config

  // do the actual ray tracing
  obs.SUBPIXEL = boost::get<int>(config["OVERSAMPLING"]);
  Image<float> im;
  std::cout << "# ray-tracing ..." << std::endl;
  try {
    center(0) = boost::get<double>(config["POINTING_X"]);
    center(1) = boost::get<double>(config["POINTING_Y"]);
    obs.makeImage(im,&center);
  } catch (std::invalid_argument) {
    obs.makeImage(im);
  }

  // write output
  fitsfile* fptr = FITS::createFile(outfile);
  FITS::writeImage(fptr,im);

  // add elementary WCS parameters
  FITS::updateKeyword(fptr,"WCSAXES",2);
  std::string val = "FK5";
  FITS::updateKeyword(fptr,"RADECSYS",val);
  FITS::updateKeyword(fptr,"EQUINOX",2000.);
  val = "RA---TAN";
  FITS::updateKeyword(fptr,"CTYPE1",val);
  val = "DEC--TAN";
  FITS::updateKeyword(fptr,"CTYPE2",val);
  FITS::updateKeyword(fptr,"CRVAL1",im.grid(0,0)/3600);
  FITS::updateKeyword(fptr,"CRVAL2",im.grid(0,1)/3600);
  FITS::updateKeyword(fptr,"CRPIX1",0.);
  FITS::updateKeyword(fptr,"CRPIX2",0.);
  val = "deg";
  FITS::updateKeyword(fptr,"CUNIT1",val);
  FITS::updateKeyword(fptr,"CUNIT2",val);
  FITS::updateKeyword(fptr,"CDELT1",tel.pixsize/3600);
  FITS::updateKeyword(fptr,"CDELT2",tel.pixsize/3600);
  FITS::closeFile(fptr);


  t1 = time(NULL);
  std::cout << "# computation time: " << t1-t0 << " seconds" << std::endl;
}
