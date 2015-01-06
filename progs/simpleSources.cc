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
  std::cout << "# SkyLens++: simpleSources v" << cmd.getVersion() << " (git" << STRINGIFY(GITREV) << ")" << std::endl;
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
  sourcecat.config.read(ifs);
  sourcecat.imref.pixsize = boost::get<double>(sourcecat.config["PIXSIZE"]);
  sourcecat.imref.fov = boost::get<double>(sourcecat.config["FOV"]);
  SourceCatalog::Band b;
  b.name = "tel";
  sourcecat.imref.bands.insert(b);

  // get readshifts from config file
  std::vector<double> redshifts = boost::get<std::vector<double> > (sourcecat.config["REDSHIFT"]);
  // total number of gals in FoV
  double N = fov(0)*fov(1) / 3600 * density.getValue();
  // avg. distance beween N gals in FoV
  int L = (int) floor(sqrt(N));
  std::cout << "# creating " << int(floor(N)) << " sources in FoV" << std::endl;

  // create DB table to store sersic information
  std::string dbfile = boost::get<std::string>(sourcecat.config["MODEL0_DBFILE"]);
  SQLiteDB* sdb;
  if (dbfile == "$APPLICATION_DB$")
    sdb = &db;
  else
    sdb->connect(dbfile);
  std::string dbtable =  boost::get<std::vector<std::string> >(sourcecat.config["MODEL0_DBTABLES"])[0];
  dbtable = split(dbtable,':')[1];
  std::string query = "DROP TABLE IF EXISTS " + dbtable + ";";
  sdb->query(query);
  query = "CREATE TABLE " + dbtable + " (id INTEGER PRIMARY KEY ASC, n_sersic FLOAT NOT NULL , radius FLOAT NOT NULL, ellipticity FLOAT NOT NULL);";
  sdb->query(query);
  sdb->exec("BEGIN TRANSACTION", NULL);

  // prepare statment to insert sersic infos in table
  query = "INSERT INTO " + dbtable + " VALUES (?,?,?,?);";
  sqlite3_stmt *stmt;
  sdb->checkRC(sqlite3_prepare_v2(sdb->db, query.c_str(), query.size(), &stmt, NULL));

  // set the galaxy infos
  GalaxyInfo info;
  info.model_type = 0;
  for (unsigned long i=0; i < N; i++) {
    info.object_id = i+1;
    if (regular.isSet()) {
      info.centroid(0) = (-0.5 + (0.5+(i%L))/L) * fov(0);
      info.centroid(1) = (-0.5 + (0.5+(i/L))/L) * fov(1);
    } else {
      info.centroid(0) = (-0.5 + gsl_rng_uniform(r)) * fov(0);
      info.centroid(1) = (-0.5 + gsl_rng_uniform(r)) * fov(1);
    }
    info.redshift = info.redshift_layer = redshifts[i%redshifts.size()];
    info.rotation = M_PI * gsl_rng_uniform(r);
    info.mag = mag.getValue();
    info.adus["tel"] = Conversion::photons2ADU(Conversion::flux2photons(Conversion::mag2flux(info.mag),1,tel,transmittance),tel.gain);
    info.sed = "none";
    info.sed_norm = 0;
    sourcecat.push_back(info);

    // add sersic information to table
    sdb->checkRC(sqlite3_bind_int(stmt,1,info.object_id));
    sdb->checkRC(sqlite3_bind_double(stmt,2,n_sersic.getValue()));
    sdb->checkRC(sqlite3_bind_double(stmt,3,radius.getValue()));
    sdb->checkRC(sqlite3_bind_double(stmt,4,gsl_ran_rayleigh (r,sigma_e.getValue()/M_SQRT2)));
    if(sqlite3_step(stmt)!=SQLITE_DONE)
      throw std::runtime_error("SourceCatalog: insertion failed: " + std::string(sqlite3_errmsg(sdb->db)));
    sdb->checkRC(sqlite3_reset(stmt));
  }

  // finalize statment
  sdb->checkRC(sqlite3_finalize(stmt));
  sdb->exec("END TRANSACTION", NULL);

  // save source catalog
  std::cout << "# saving sources" << std::endl;
  sourcecat.save(db);
}
