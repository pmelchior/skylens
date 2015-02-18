#include <shapelens/ShapeLens.h>
#include <shapelens/MathHelper.h>
#include <shapelens/LensHelper.h>
#include <shapelens/DEIMOSElliptical.h>
#include <skylens/SkyLens.h>
#include <time.h>
#include <tclap/CmdLine.h>
#include <fstream>

using namespace skylens;
using namespace shapelens;

double distance(const Point<double>& P1, const Point<double>& P2) {
  return sqrt(pow2(P1(0)-P2(0)) + pow2(P1(1)-P2(1)));
}

template <class T>
std::vector<T>& append(std::vector<T>& v1, const std::vector<T>& v2) {
  v1.insert(v1.end(), v2.begin(), v2.end());
  return v1;
}

void setTelescopeFoV(const Telescope& tel_orig, const Rectangle<double>& bb, LensingLayer* ll, double z_s, Telescope& tel, Point<double>& center) {
  Rectangle<double> search_area;
  search_area.ll(0) = -tel_orig.fov_x/2;
  search_area.ll(1) = -tel_orig.fov_x/2;
  search_area.tr(0) = tel_orig.fov_x/2;
  search_area.tr(1) = tel_orig.fov_x/2;
  Point<double> P = bb.ll;
  std::vector<Point<data_t> > thetas = ll->findImages(P, z_s, search_area);
  P(1) = bb.tr(1);
  append(thetas, ll->findImages(P, z_s, search_area));
  P(1) = bb.ll(1);
  P(0) = bb.tr(0);
  append(thetas, ll->findImages(P, z_s, search_area));
  P = bb.tr;
  append(thetas, ll->findImages(P, z_s, search_area));
  if (thetas.size() < 4)
    throw std::runtime_error("shear_accuracy::setTelecopeFoV did not find 4 images");
  
  search_area.ll = thetas[0];
  search_area.tr = thetas[0];
  for (int j = 1; j < thetas.size(); j++) {
    P = thetas[j];
    if (P(0) < search_area.ll(0))
      search_area.ll(0) = P(0);
    if (P(1) < search_area.ll(1))
      search_area.ll(1) = P(1);
    if (P(0) > search_area.tr(0))
      search_area.tr(0) = P(0);
    if (P(1) > search_area.tr(1))
      search_area.tr(1) = P(1);
  }
  center(0) = (search_area.tr(0) + search_area.ll(0))/2;
  center(1) = (search_area.tr(1) + search_area.ll(1))/2;
  tel.fov_x = 1.1*(search_area.tr(0) - search_area.ll(0)); // little padding
  tel.fov_y = 1.1*(search_area.tr(1) - search_area.ll(1));
}

Object im2obj_centroid(const Image<float>& im, Point<double>& theta_pixel) {
  Object obj;
  obj.resize(im.size());
  obj.grid.setSize(im.grid.getStartPosition(0), im.grid.getStartPosition(1), im.grid.getSize(0), im.grid.getSize(1));
  double sum_flux = 0;
  theta_pixel(0) = theta_pixel(1) = 0;
  for (int i=0; i < im.size(); i++) {
    obj(i) = im(i);
    if (obj(i) > 0) {
      Point<int> coords = im.grid.getCoords(i);
      sum_flux += obj(i);
      Point<double> pos(coords(0)*obj(i), coords(1)*obj(i));
      theta_pixel += pos;
    }
  }
  theta_pixel /= sum_flux;
  return obj;
}

int main(int argc, char* argv[]) {
  // parse commandline
  TCLAP::CmdLine cmd("Compare shears between full raytracing and first-order lensing only", ' ', "0.4");
  TCLAP::ValueArg<std::string> configfile("c","config","Configuration file",true,"","string", cmd);
  TCLAP::SwitchArg output("o","outfile","Whether simulated images should be written out",cmd, false);
  TCLAP::ValueArg<double> z_s("z","z_s","source redshift",true,1,"double", cmd);
  TCLAP::ValueArg<unsigned int> N("N","number","Number of different source positions",true,100,"unsigned int", cmd);
  TCLAP::ValueArg<unsigned int> seed_lens("s","seed","RNG seed for lens simulations",false,0,"unsigned int", cmd);
  cmd.parse(argc,argv);


  // for measuring computation time
  time_t t0,t1;
  t0 = time(NULL);
  std::cout << "# shear_accuracy v" << cmd.getVersion() << " (git" << STRINGIFY(GITREV) << ")" << std::endl;
  
  // read in global config file
  std::cout << "# reading global configuration from " << configfile.getValue() << std::endl;
  Property config;
  std::ifstream ifs(configfile.getValue().c_str());
  config.read(ifs);
  ifs.close();

  // set outfile: no catch, this must be set
  std::string fileroot = boost::get<std::string>(config["PROJECT"]);
  std::string outfile = boost::get<std::string>(config["OUTFILE"]);

  // connect to application DB
  SQLiteDB& db = ApplicationDB::getInstance();
  db.connect(fileroot+".db");

  // create result table
  std::string query =  "CREATE TABLE IF NOT EXISTS shear_accuracy (";
  query += "mass double NOT NULL,";
  query += "zl double NOT NULL,";
  query += "zs double NOT NULL,";
  query += "seed_lens int NOT NULL,"; // unique identifier for lens sim
  query += "seed int NOT NULL,";
  query += "Rs double NOT NULL,";
  query += "ns double NOT NULL,";
  query += "eps double NOT NULL,";
  query += "r double NOT NULL,";
  query += "Rc double NOT NULL,"; // distance from critical curve
  query += "g double NOT NULL,";
  query += "kappa double NOT NULL,";
  query += "eps_mo_t double NOT NULL,";
  query += "eps_w_t double NOT NULL,";
  query += "eps_pred_t double NOT NULL,";
  query += "eps_pred_mo_t double NOT NULL);";
  db.query(query);
  
  // get datapath
  std::string datapath = getDatapath();
  
  // set telescope and filter
  std::cout << "# setting up Telescope and Observation" << std::endl;
  Telescope tel(boost::get<std::string>(config["TELESCOPE"]),
		boost::get<std::string>(config["FILTER"]));
  Telescope tel_orig = tel;
  int exptime = boost::get<int>(config["EXPTIME"]);
  Observation obs(tel,exptime);
  obs.SUBPIXEL = boost::get<int>(config["OVERSAMPLING"]);

  // set global cosmology: default is vanilla CDM
  Cosmology& cosmo = SingleCosmology::getInstance();
  try {
    std::string cosmofile = boost::get<std::string>(config["COSMOLOGY"]);
    // FIXME: implementation missing
  } catch (std::invalid_argument) {}

  // set global RNG seed if demanded
  RNG& rng = Singleton<RNG>::getInstance();
  const gsl_rng * r = rng.getRNG();
  unsigned long seed;
  try {
    seed =  boost::get<int>(config["RNG_SEED"]);
    std::cout << "# setting RNG seed from config file to " << seed << std::endl;
  } catch (std::invalid_argument) {
    // returns long, but need int for config file: implicit mapping by overrun
    seed = abs(int(time (NULL) * getpid()));
    std::cout << "# setting RNG seed to " << seed << std::endl;
  }
  gsl_rng_set(r,seed);


  // read in lens config
  std::vector<std::string> files;
  files = boost::get<std::vector<std::string> > (config["LENSES"]);
  if (files.size() != 1)
    throw std::invalid_argument("shear_accuracy: specify exactly one lens layer in config parameter LENSES");
  std::cout << "# setting up lens from " << files[0] << std::endl;
  test_open(ifs,datapath,files[0]);
  Property lensconfig;
  lensconfig.read(ifs);

  // create lens layer
  std::string anglefile = boost::get<std::string>(lensconfig["ANGLEFILE"]);
  test_open(ifs,datapath,anglefile);
  // get mass from file
  double mass = -1;
  try {
    fitsfile* fptr = FITS::openFile(anglefile);
    FITS::readKeyword(fptr, "MVIR", mass);
    FITS::closeFile(fptr);
  } catch (std::invalid_argument) {}

  Point<double> center_lens(0,0);
  LensingLayer* ll;
  try {
    center_lens(0) = boost::get<double>(lensconfig["POS_X"]);
    center_lens(1) = boost::get<double>(lensconfig["POS_Y"]);
    ll = new LensingLayer(boost::get<double>(lensconfig["REDSHIFT"]), anglefile, &center_lens);
  } catch (std::invalid_argument) {
    ll = new LensingLayer(boost::get<double>(lensconfig["REDSHIFT"]), anglefile);
  }
  std::vector<Point<double> > critical = ll->findCriticalPoints(z_s.getValue(), 1);

  if (critical.size() == 0)
    throw std::runtime_error("shear_accuracy: no critical points found, cannot create arcs!");

  // create a layer with only one source
  // repeat N times
  SourceModelList models;
  LayerStack& ls = SingleLayerStack::getInstance();
  complex<double> I(0,1);
    
  // get caustic points
  std::vector<Point<double> > caustic;
  for (std::vector<Point<double> >::const_iterator iter = critical.begin(); iter != critical.end(); iter++)
    caustic.push_back(ll->getBeta(*iter, z_s.getValue()));

  // outputs
  fitsfile* fptr;
  if (output.isSet())
    fptr = FITS::createFile(outfile);
  // prepare statement
  sqlite3_stmt *stmt;
  query = "INSERT INTO shear_accuracy VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?);";
  db.checkRC(sqlite3_prepare_v2(db.db, query.c_str(), query.size(), &stmt, NULL));
  // prevent excessive I/O
  db.exec("BEGIN TRANSACTION", NULL);

  // COSMOS morphology catalog with Sersic models
  fitsfile* cat = FITS::openTable("/n/des/pmelchior/skylens/data_ext/sources/COSMOS/totalRAW00000.29949.fits.gz");
  int cat_len = FITS::getTableRows(cat);
  double sersic_fit[8];
  int sersic_col = FITS::getTableColumnNumber(cat, "SERSICFIT");

  // generate source model(s)
  for (int n=0; n < N.getValue(); n++) {
    data_t Rs = 0, q, ns;
    // exclude failed fits with radius = 0
    // and those to large to be plausibly lensed
    while (Rs == 0. || Rs > 2.)  {
      // pick random galaxy
      FITS::readTableValue(cat, int(floor(cat_len*gsl_rng_uniform (r))), sersic_col, *sersic_fit, 8);
      q = sersic_fit[3];
      ns = sersic_fit[2];
      Rs = sersic_fit[1] * 0.05 * sqrt(q); // correct pixel scale and ellipticity
    }
    complex<data_t> eps((1-q)/(1+q), 0);
    eps *= exp(2.*I*M_PI*gsl_rng_uniform(r));
    data_t flux = 1;
    data_t trunc_radius = 5;
    // source centroid:
    // randomly select point from fieldsize - border region
    Point<double> limit(0.5*tel_orig.fov_x, 0.5*tel_orig.fov_y);
    Point<data_t> beta (-limit(0)/2 + limit(0)*gsl_rng_uniform (r),
			-limit(1)/2 + limit(1)*gsl_rng_uniform (r));
    ShiftTransformation Z(beta);
    models.push_back(boost::shared_ptr<SourceModel>(new SersicModel(ns, Rs, flux, eps, trunc_radius, &Z)));

    // Before even raytracing: 
    // 1) check if image center is in FoV
    // will also be used later as the theoretical center of the source
    Rectangle<double> search_area;
    search_area.ll(0) = -tel_orig.fov_x/2;
    search_area.ll(1) = -tel_orig.fov_x/2;
    search_area.tr(0) = tel_orig.fov_x/2;
    search_area.tr(1) = tel_orig.fov_x/2;
    std::vector<Point<data_t> > thetas = ll->findImages(beta, z_s.getValue(), search_area);
    Point<data_t> center(0,0);
    bool inside = true, multiple_images = false;
    if (thetas.size() == 0) 
      inside = false;
    if (thetas.size() > 1)
      multiple_images = true;

    if (inside || !multiple_images) {
      // to speed up the process: adjust telescope to location and size of
      // actual image
      // estimate this by transforming the model's bounding box
      try {
	setTelescopeFoV(tel_orig, models[0]->getSupport(), ll, z_s.getValue(), tel, center);
      } catch (std::runtime_error) {
	inside = false;
      }

      // 2) check if multiple images will occur
      // in which case we'll do another galaxy
      bool multiple_images = false;
      for (std::vector<Point<double> >::const_iterator iter = caustic.begin(); iter != caustic.end(); iter++) {
	if (models[0]->contains(*iter)) {
	  multiple_images = true;
	  break;
	}
      }
    }

    if (!inside || multiple_images) {
      if (!inside)
	std::cout << "# Source not in FoV: rejected" << std::endl;
      if (multiple_images)
	std::cout << "# Source creates multiple images: rejected" << std::endl;
      n -= 1;
    }
    else { // single image
      Image<float> im;
      GalaxyLayer* gl = new GalaxyLayer(z_s.getValue(), models);
      // do the actual ray tracing
      obs.makeImage(im,&center);
      
      // convert to Object and find centroid in theta space
      Point<double> theta_pixel;
      Object obj = im2obj_centroid(im, theta_pixel);
      // measure 2nd moments around observed theta centroid
      Moments mo(obj, FlatWeightFunction(), 2, &theta_pixel);
      std::complex<double> eps_mo = shapelens::epsilon(mo);

      // now take the true center in theta frame and emulate
      // pure first-order lensing measurement
      Point<double> theta = thetas[0]; // we already know there is only one image
      theta_pixel = theta;
      im.grid.getWCS().inverse_transform(theta_pixel);
      std::complex<double> gamma = ll->getShear(theta, z_s.getValue(), true);
      std::complex<double> eps_pred = eps;
      shapelens::lensEps(gamma, eps_pred);
      
      // compute distance to nearest point on critical curve
      double Rc = tel_orig.fov_x + tel_orig.fov_y;
      for (std::vector<Point<double> >::const_iterator iter = critical.begin(); iter != critical.end(); iter++)
	Rc = std::min(Rc, distance(*iter, theta));

      // cross-check: sample model with correct ellipticity on obj grid
      // need to correct for pixel size and magnification
      Object obj_pred = obj;
      ShiftTransformation T(theta_pixel);
      data_t Rs_pixel = Rs / tel.pixsize;
      data_t flux_pixel = flux / pow2(tel.pixsize);
      data_t kappa = ll->getConvergence(theta, z_s.getValue());
      data_t mu = 1./(pow2(1 - kappa) - pow2(abs(gamma)*(1-kappa)));
      Rs_pixel *= sqrt(mu);
      flux_pixel *= mu;
      SersicModel pred(ns, Rs_pixel, flux_pixel, eps_pred, trunc_radius, &T);
      setObject(pred, obj_pred, obs.SUBPIXEL);
      Moments mo_pred(obj_pred, FlatWeightFunction(), 2, &theta_pixel);
      std::complex<double> eps_pred_mo = shapelens::epsilon(mo_pred);
      
      // 2nd cross-check: used DEIMOS to measure the weighted
      // 2nd-order moments of the actual image 
      // with a scale given by the predicted post-lensing Rs_pixel
      // at the location theta of the predicted center
      obj.centroid = theta;
      im.grid.getWCS().inverse_transform(obj.centroid);
      obj.noise_rms = 1;
      obj.noise_mean = 0;
      DEIMOSElliptical d(obj, 2, 4, Rs_pixel);
      std::complex<double> eps_w = shapelens::epsilon(d.mo);
      
      // output simulated images
      if (output.isSet()) {
	// select the bad cases only
	if (abs(eps_pred_mo - eps_pred) > 0.01) {
	  FITS::writeImage(fptr,obj);
	  FITS::writeImage(fptr,obj_pred);
	}
      }
      
      // write results to DB
      db.checkRC(sqlite3_bind_double(stmt, 1, mass));
      db.checkRC(sqlite3_bind_double(stmt, 2, ll->getRedshift()));
      db.checkRC(sqlite3_bind_double(stmt, 3, z_s.getValue()));
      db.checkRC(sqlite3_bind_int(stmt, 4, seed_lens.getValue()));
      db.checkRC(sqlite3_bind_int(stmt, 5, seed));
      db.checkRC(sqlite3_bind_double(stmt, 6, Rs));
      db.checkRC(sqlite3_bind_double(stmt, 7, ns));
      db.checkRC(sqlite3_bind_double(stmt, 8, abs(eps)));
      db.checkRC(sqlite3_bind_double(stmt, 9, sqrt(pow2(theta(0) - center_lens(0)) + pow2(theta(1) - center_lens(1)))));
      db.checkRC(sqlite3_bind_double(stmt, 10, Rc));
      db.checkRC(sqlite3_bind_double(stmt, 11, shapelens::epsTangential(gamma, theta, center_lens)));
      db.checkRC(sqlite3_bind_double(stmt, 12, kappa));
      db.checkRC(sqlite3_bind_double(stmt, 13, shapelens::epsTangential(eps_mo, theta, center_lens)));
      db.checkRC(sqlite3_bind_double(stmt, 14, shapelens::epsTangential(eps_w, theta, center_lens)));
      db.checkRC(sqlite3_bind_double(stmt, 15, shapelens::epsTangential(eps_pred, theta, center_lens)));
      db.checkRC(sqlite3_bind_double(stmt, 16, shapelens::epsTangential(eps_pred_mo, theta, center_lens)));
      if(sqlite3_step(stmt)!=SQLITE_DONE)
	std::cout << "# shear_accuracy: insertion failed: " << std::string(sqlite3_errmsg(db.db)) << std::endl;
      db.checkRC(sqlite3_reset(stmt));
      
      // delete GalaxyLayer and remove from LayerStack
      // such that its not present in the following images
      for (LayerStack::iterator iter = ls.begin(); iter != ls.end(); iter++) {
	if (iter->second->getType() == "SG") {
	  ls.erase(iter);
	  break;
	}
      }

      delete gl;

      // save snapshots along the way
      if (n%100 == 0) {
	db.exec("END TRANSACTION", NULL);
	db.exec("BEGIN TRANSACTION", NULL);
	}
    }

    // clean up models
    models.clear();
  }

  // finalize writing and close output files
  db.checkRC(sqlite3_finalize(stmt));
  db.exec("END TRANSACTION", NULL);
  if (output.isSet())
    FITS::closeFile(fptr);

  t1 = time(NULL);
  std::cout << "# computation time: " << t1-t0 << " seconds" << std::endl;
}
