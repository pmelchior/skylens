#include <shapelens/ShapeLens.h>
#include <skylens/SkyLens.h>
#include <tclap/CmdLine.h>
#include <fstream>

using namespace skylens;
using namespace shapelens;


int main(int argc, char* argv[]) {
  // parse commandline
  TCLAP::CmdLine cmd("Create simplistic sources for  SkyLens++ simulator", ' ', "0.2");
  TCLAP::ValueArg<std::string> configfile("c","config","Configuration file",true,"","string", cmd);
  TCLAP::ValueArg<data_t> mag("m","mag","Magnitude of sources",true,25,"data_t",cmd);
  TCLAP::ValueArg<data_t> sigma_e("e","sigma_e","Intrinsic ellipticity rms",true,0.3,"data_t",cmd);
  TCLAP::ValueArg<data_t> radius("r","radius","Effective radius [arcsec]",true,0.35,"data_t",cmd);
  TCLAP::SwitchArg regular("R","regular","Regular placement of source", cmd, false);
  TCLAP::ValueArg<data_t> n_sersic("n","n_sersic","Sersic index",true,1.5,"data_t",cmd);
  TCLAP::ValueArg<data_t> density("d","density","Number density [arcmin^-2] of sources",true,100,"data_t",cmd);
  cmd.parse(argc,argv);
  
  // read in global config file
  Property config;
  std::ifstream ifs(configfile.getValue().c_str());
  config.read(ifs);
  ifs.close();

  // set outfile: no catch, this must be set
  std::string outfile = boost::get<std::string>(config["OUTFILE"]);
  std::string fileroot = boost::get<std::string>(config["PROJECT"]);

  // connect to application DB
  SQLiteDB db;
  db.connect(fileroot+".db");

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
  RNG& rng = Singleton<RNG>::getInstance();
  const gsl_rng * r = rng.getRNG();
  try {
    int seed =  boost::get<int>(config["RNG_SEED"]);
    gsl_rng_set(r,seed);
  } catch (std::invalid_argument) {}

  // get sources from config files
  std::vector<std::string> files = boost::get<std::vector<std::string> >(config["SOURCES"]);
  test_open(ifs,datapath,files[0]);
  
  // set up source catalog
  SourceCatalog sourcecat;
  sourcecat.imref.pixsize = 1;
  sourcecat.imref.fov = 0;
  SourceCatalog::Band b;
  b.name = "tel";
  b.overlap = 1;
  sourcecat.imref.bands.insert(b);
  sourcecat.config.read(ifs);

  // get readshifts from config file
  std::vector<double> redshifts = boost::get<std::vector<double> > (sourcecat.config["REDSHIFT"]);
  // adjust numbers
  double Re = radius.getValue() * boost::get<double>(sourcecat.config["PIXSIZE"]);
  // total number of gals in FoV
  double N = fov(0)*fov(1) / 3600 * density.getValue();
  // avg. distance beween N gals in FoV
  int L = (int) floor(sqrt(N));

  // set the galaxy infos
  GalaxyInfo info;
  info.model_type = 0;
  for (unsigned long i=0; i < N; i++) {
    info.object_id = i;
    if (regular.isSet()) {
      info.centroid(0) = (0.5+(i%L))/L * fov(0);
      info.centroid(1) = (0.5+(i/L))/L * fov(1);
    } else {
      info.centroid(0) = gsl_rng_uniform(r) * fov(0);
      info.centroid(1) = gsl_rng_uniform(r) * fov(1);
    }
    info.redshift = info.redshift_layer = redshifts[i%redshifts.size()];
    info.radius = Re;
    info.ellipticity = gsl_ran_rayleigh (r,sigma_e.getValue()/M_SQRT2);
    info.rotation = M_PI * gsl_rng_uniform(r);
    info.n_sersic = n_sersic.getValue();
    info.mag = mag.getValue();
    info.adus["tel"] = Conversion::photons2ADU(Conversion::flux2photons(Conversion::mag2flux(info.mag),1,tel,transmittance),tel.gain);
    info.sed = "none";
    info.sed_norm = 0;
    sourcecat.push_back(info);
  }
  sourcecat.save(db);
}
