#include <shapelens/ShapeLens.h>
#include <skylens/SkyLens.h>
#include <time.h>
#include <tclap/CmdLine.h>
#include <fstream>

using namespace skylens;
using namespace shapelens;

int main(int argc, char* argv[]) {
  // parse commandline
  TCLAP::CmdLine cmd("Create arcs with SkyLens++", ' ', "0.1");
  TCLAP::ValueArg<std::string> configfile("c","config","Configuration file",true,"","string", cmd);
  TCLAP::ValueArg<double> z_s("z","z_s","source redshift",true,1,"double", cmd);
  TCLAP::ValueArg<double> delta("d","delta","displacement from caustic [arcsec]",true,0.5,"double", cmd);
  TCLAP::ValueArg<unsigned int> N("N","number","Number of different source positions",true,100,"unsigned int", cmd);
  TCLAP::ValueArg<double> sigma_e("e","sigma_e","Source ellipticity dispersion",false,0.3,"double",cmd);
  cmd.parse(argc,argv);


  // for measuring computation time
  time_t t0,t1;
  t0 = time(NULL);
  std::cout << "# SkyLens++ v" << cmd.getVersion() << " (svn" << STRINGIFY(SVNREV) << ")" << std::endl;
  
  // read in global config file
  std::cout << "# reading global configuration from " << configfile.getValue() << std::endl;
  Property config;
  std::ifstream ifs(configfile.getValue().c_str());
  config.read(ifs);
  ifs.close();

  // set outfile: no catch, this must be set
  std::string fileroot = boost::get<std::string>(config["PROJECT"]);
  std::string outfile = boost::get<std::string>(config["OUTFILE"]);
  fitsfile* fptr = IO::createFITSFile(outfile);

  /*// connect to application DB
  SQLiteDB& db = ApplicationDB::getInstance();
  db.connect(fileroot+".db");*/

  // get datapath
  std::string datapath = getDatapath();
  
  // set telescope and filter
  std::cout << "# setting up Telescope and Obervation" << std::endl;
  Telescope tel(boost::get<std::string>(config["TELESCOPE"]),
		boost::get<std::string>(config["FILTER"]));
  int exptime = boost::get<int>(config["EXPTIME"]);
  Observation obs(tel,exptime);
  obs.SUBPIXEL = boost::get<int>(config["OVERSAMPLING"]);

  // set (top-right corner of ) the global FoV
  // default = telescope's FoV
  Point<double> fov(tel.fov_x,tel.fov_y);
  try {
    fov(0) = boost::get<double>(config["GLOBAL_FOV_X"]);
    fov(1) = boost::get<double>(config["GLOBAL_FOV_Y"]);
  } catch (std::invalid_argument) {}

  cosmology& cosmo = SingleCosmology::getInstance();

  // set global RNG seed if demanded
  RNG& rng = Singleton<RNG>::getInstance();
  const gsl_rng * r = rng.getRNG();
  try {
    int seed =  boost::get<int>(config["RNG_SEED"]);
    gsl_rng_set(r,seed);
  } catch (std::invalid_argument) {}


  // read in lens config
  std::vector<std::string> files;
  Point<double> center;
  std::map<Point<double>, Point<double> > cpoints;
  std::map<Point<double>, Point<double> >::iterator citer;
  files = boost::get<std::vector<std::string> > (config["LENSES"]);
  if (files.size() != 1)
    throw std::invalid_argument("createArcs: specify exactly one lens layer in config parameter LENSES");
  std::cout << "# setting up lens from " << files[0] << std::endl;
  test_open(ifs,datapath,files[0]);
  Property lensconfig;
  lensconfig.read(ifs);
  // create lens layer
  center(0) = boost::get<double>(lensconfig["POS_X"]);
  center(1) = boost::get<double>(lensconfig["POS_Y"]);
  std::string anglefile = boost::get<std::string>(lensconfig["ANGLEFILE"]);
  test_open(ifs,datapath,anglefile);
  LensingLayer* ll = new LensingLayer(boost::get<double>(lensconfig["REDSHIFT"]), anglefile, center);
  cpoints = ll->findCriticalPoints(z_s.getValue());

  // create a layer with only one source lying close to caustic
  // repeat N times
  shapelens::SourceModelList models;
  LayerStack& ls = SingleLayerStack::getInstance();
  complex<double> I(0,1);
  for (int n=0; n < N.getValue(); n++) {
    // morphological distributions of ellipticit, sersic index, and size
    complex<double> eps(1,0);
    while (abs(eps) >= 0.99) { // restrict to viable ellipticities
      eps = gsl_ran_rayleigh(r,sigma_e.getValue()/M_SQRT2);
      eps *= exp(2.*I*M_PI*gsl_rng_uniform(r));
    }
    data_t ns = 0.5 + 4*gsl_rng_uniform (r);
    data_t Re = 0.1 + 0.9*gsl_rng_uniform (r); // size in arcsec

    // source centroid:
    // randomly select point from caustic, offset slightly via delta
    unsigned int adv = cpoints.size()*gsl_rng_uniform (r);
    citer = cpoints.begin();
    std::advance(citer,adv);
    Point<double> centroid = citer->second;
    centroid(0) += gsl_ran_gaussian(r, delta.getValue());
    centroid(1) += gsl_ran_gaussian(r, delta.getValue());
    shapelens::ShiftTransformation Z(centroid);
    models.push_back(boost::shared_ptr<shapelens::SourceModel>(new shapelens::SersicModel(ns,Re,1,eps,&Z)));
    GalaxyLayer* gl = new GalaxyLayer(z_s.getValue(), models);

      // do the actual ray tracing
    Image<float> im;
    std::cout << "# ray-tracing image " << n << " ..." << std::endl;
    try {
      center(0) = boost::get<double>(config["POINTING_X"]);
      center(1) = boost::get<double>(config["POINTING_Y"]);
      obs.makeImage(im,&center);
    } catch (std::invalid_argument) {
      obs.makeImage(im);
    }
    
    // write output
    IO::writeFITSImage(fptr,im);
    // add source parameters
    complex<double> cent(centroid(0), centroid(1));
    IO::updateFITSKeyword(fptr,"CENTROID",cent);
    IO::updateFITSKeyword(fptr,"SERSIC_N",ns);
    IO::updateFITSKeyword(fptr,"SIZE",Re);
    IO::updateFITSKeyword(fptr,"EPSILON",eps);
    // add elementary WCS parameters
    IO::updateFITSKeyword(fptr,"WCSAXES",2);
    IO::updateFITSKeywordString(fptr,"RADECSYS","FK5");
    IO::updateFITSKeyword(fptr,"EQUINOX",2000.);
    IO::updateFITSKeywordString(fptr,"CTYPE1","RA---TAN");
    IO::updateFITSKeywordString(fptr,"CTYPE2","DEC--TAN");
    IO::updateFITSKeyword(fptr,"CRVAL1",im.grid(0,0)/3600);
    IO::updateFITSKeyword(fptr,"CRVAL2",im.grid(0,1)/3600);
    IO::updateFITSKeyword(fptr,"CRPIX1",0.);
    IO::updateFITSKeyword(fptr,"CRPIX2",0.);
    IO::updateFITSKeywordString(fptr,"CUNIT1","deg");
    IO::updateFITSKeywordString(fptr,"CUNIT2","deg");
    IO::updateFITSKeyword(fptr,"CDELT1",tel.pixsize/3600);
    IO::updateFITSKeyword(fptr,"CDELT2",tel.pixsize/3600);

    // delete GalaxyLayer and remove from LayerStack
    // such that its not present in the following images
    for (LayerStack::iterator iter = ls.begin(); iter != ls.end(); iter++) {
      if (iter->second->getType() == "SG") {
	ls.erase(iter);
	break;
      }
    }
    delete gl;
    
    // delete galaxy from models
    models.clear();
  }

  IO::closeFITSFile(fptr);

  t1 = time(NULL);
  std::cout << "# computation time: " << t1-t0 << " seconds" << std::endl;
}
