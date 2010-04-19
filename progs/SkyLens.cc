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

  // set emission of the sky
  try {
    std::string sky = boost::get<std::string>(config["SKY"]);
    test_open(ifs,datapath,sky);
    obs.createSkyFluxLayer(sed(sky,"/"));
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
  SourceCatalog sourcecat;
  std::vector<std::string> files = boost::get<std::vector<std::string> >(config["SOURCES"]);
  for (int i=0; i< files.size(); i++) {
    test_open(ifs,datapath,files[i]);
    // compute source information from DB
    if (!useSources.isSet()) {
      sourcecat = SourceCatalog(files[i]);
      // account for change of FoV from reference to telescope/global
      sourcecat.adjustNumber(fov);
      // place them randomly in the FoV 
      // and on the available redshifts of GalaxyLayers
      sourcecat.distribute(fov);
      // find reference bands with overlap to total transmittance
      sourcecat.selectOverlapBands(transmittance);
      // compute flux of source in each of the remaining bands
      sourcecat.computeADUinBands(tel,transmittance);
      // save source catalogs, if demanded
      if (saveSources.isSet()) {
	// save source catalogs
	sourcecat.save(db,i);
      }
    }
    else { // use precomputed sources
      sourcecat = SourceCatalog(db,i);
    }
    // create GalaxyLayers from sources
    sourcecat.createGalaxyLayers(exptime);
  }

  // read in lens config
  Point<double> center;
  try {
    files = boost::get<std::vector<std::string> > (config["LENSES"]);
    for (int i=0; i < files.size(); i++) {
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
  center.clear();
  try {
    center(0) = boost::get<double>(config["POINTING_X"]);
    center(1) = boost::get<double>(config["POINTING_Y"]);
  } catch (std::invalid_argument) {}
  Image<float> im;
  obs.makeImage(im,center);

  // write output
  fitsfile* fptr = IO::createFITSFile(outfile);
  IO::writeFITSImage(fptr,im);

  // add elementary WCS parameters
  IO::updateFITSKeyword(fptr,"WCSAXES",2);
  IO::updateFITSKeywordString(fptr,"RADECSYS","FK5");
  IO::updateFITSKeyword(fptr,"EQUINOX",2000.);
  IO::updateFITSKeywordString(fptr,"CTYPE1","RA---TAN");
  IO::updateFITSKeywordString(fptr,"CTYPE2","DEC--TAN");
  IO::updateFITSKeyword(fptr,"CRVAL1",0.);
  IO::updateFITSKeyword(fptr,"CRVAL2",0.);
  IO::updateFITSKeyword(fptr,"CRPIX1",0.);
  IO::updateFITSKeyword(fptr,"CRPIX2",0.);
  IO::updateFITSKeywordString(fptr,"CUNIT1","deg");
  IO::updateFITSKeywordString(fptr,"CUNIT2","deg");
  IO::updateFITSKeyword(fptr,"CDELT1",tel.pixsize/3600);
  IO::updateFITSKeyword(fptr,"CDELT2",tel.pixsize/3600);
  IO::closeFITSFile(fptr);


  t1 = time(NULL);
  std::cout << "Computation time: " << t1-t0 << " seconds" << std::endl;
}
