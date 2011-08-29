#include <cstring>
#include "../include/SourceCatalog.h"
#include "../include/RNG.h"
#include "../include/Layer.h"
#include "../include/Conversion.h"
#include <boost/lexical_cast.hpp>
#include <shapelens/shapelets/ShapeletObjectList.h>
#include <shapelens/utils/Property.h>

namespace skylens {
  bool SourceCatalog::Band::operator<(const Band& b) const {
    if (curve.lambdaEff() < b.curve.lambdaEff())
      return true;
    else 
      return false;
  }

  SourceCatalog::SourceCatalog() :
    std::list<GalaxyInfo> () {
    replication_ratio = 0;
  }
  
  SourceCatalog::SourceCatalog(std::string configfile) :
    std::list<GalaxyInfo> () {
    
    // read in contents of config file
    parseConfig(configfile);
    // connect to DB
    shapelens::SQLiteDB db;
    db.connect(boost::get<std::string>(config["DBFILE"]));

    // query database
    shapelens::SQLiteDB::SQLiteResult dbr = db.query(query);
    char** row;
    GalaxyInfo info;
    info.sed_norm = 0;

    // build catalog from DB entries
    while ((row = dbr.getRow()) != NULL) {
      for (int i =0; i < dbr.getFieldCount(); i++) {
	std::string name = dbr.getFieldName(i);
	if (name == "id")
	  info.object_id = boost::lexical_cast<unsigned long>(row[i]);
	else if (name == "centroid_x") {
	  if (row[i] != NULL)
	    info.centroid(0) = boost::lexical_cast<double>(row[i]);
	  else
	    info.centroid(0) = 0;
	}
	else if (name == "centroid_y") {
	  if (row[i] != NULL)
	    info.centroid(1) = boost::lexical_cast<double>(row[i]);
	  else
	    info.centroid(1) = 0;
	}
	else if (name.substr(0,3) == "mag") {
	  for (std::set<Band>::iterator iter = imref.bands.begin(); iter != imref.bands.end(); iter++) {
	    if (name == "mag_" + iter->name) {
	      if (row[i] != NULL) {
		if (dbr.getFieldName(i+1) == name+"_error") {
		  info.mags[iter->name] = 
		    std::pair<double,double> (boost::lexical_cast<double>(row[i]),
					      boost::lexical_cast<double>(row[i+1]));
		  i++; // skip next field: error of mag in current band
		} else { // next field is not mag error: set error = 0
		  info.mags[iter->name] = 
		    std::pair<double,double> (boost::lexical_cast<double>(row[i]), 0);
		}
	      }
	    }
	  }
	}
	else if (name == "redshift") {
	  if (row[i] != NULL)
	    info.redshift = boost::lexical_cast<double>(row[i]);
	  else
	    info.redshift = 0;
	}
	else if (name == "sed") {
	  if (row[i] != NULL)
	    info.sed = std::string(row[i]);
	  else
	    info.sed = "";
	}
	else if (name == "model_type") {
	  info.model_type = boost::lexical_cast<double>(row[i]);
	  need_model.insert(info.model_type);
	  // check that models are available for all required types
	  if (model_band_tables.find(info.model_type) == model_band_tables.end())
	    throw std::invalid_argument("SourceCatalog: model tables found for model_type = " + std::string(row[i]));
	}
      }
      SourceCatalog::push_back(info);
    }
    replication_ratio = 1;
  }
    
  SourceCatalog::SourceCatalog(shapelens::SQLiteDB& db, int i) :
    std::list<GalaxyInfo> () {
    
    // get config from db
    std::ostringstream tablename;
    tablename << "sources_" << i;
    std::string query = "SELECT config FROM config WHERE name='" + tablename.str() + "';";
    sqlite3_stmt *stmt;
    db.checkRC(sqlite3_prepare_v2(db.db, query.c_str(), query.size(), &stmt, NULL));
    if (sqlite3_step(stmt) == SQLITE_ROW) {
      std::stringstream is;
      is << sqlite3_column_text(stmt,0);
      config.read(is);
      // parse contents of config
      parseConfig();
    } else
      throw std::runtime_error("SourceCatalog: no config found DB for table " + tablename.str());
	
    replication_ratio = 0;

    // get contents from each source table
    GalaxyInfo info;
    query = "SELECT * FROM " + tablename.str() +";";
    db.checkRC(sqlite3_finalize(stmt));
    db.checkRC(sqlite3_prepare_v2(db.db, query.c_str(), query.size(), &stmt, NULL));
    while (sqlite3_step(stmt) == SQLITE_ROW) {
      info.object_id = sqlite3_column_int(stmt,0);
      info.redshift = sqlite3_column_double(stmt,1);
      info.sed = std::string(reinterpret_cast<const char*>(sqlite3_column_text(stmt,2)));
      info.sed_norm = sqlite3_column_double(stmt,3);
      info.mag = sqlite3_column_double(stmt,4);
      info.centroid(0) = sqlite3_column_double(stmt,5);
      info.centroid(1) = sqlite3_column_double(stmt,6);
      info.rotation = sqlite3_column_double(stmt,7);
      info.redshift_layer = sqlite3_column_double(stmt,8);
      info.model_type = sqlite3_column_int(stmt,9);
      need_model.insert(info.model_type);
      info.adus.clear();
      std::string band(reinterpret_cast<const char*>(sqlite3_column_text(stmt,10)));
      info.adus[band] = sqlite3_column_double(stmt,11);
      SourceCatalog::push_back(info);
    }
    db.checkRC(sqlite3_finalize(stmt));
  }

  void SourceCatalog::parseConfig(std::string configfile) {
    // open source config file
    std::ifstream ifs;
    std::string datapath = getDatapath();
    if (configfile.size() > 0) { // config file needs to be read out
      test_open(ifs,datapath,configfile);
      config.clear();
      config.read(ifs);
      ifs.close();
    }

    // check path to master table
    try {
      std::string dbfile = boost::get<std::string>(config["DBFILE"]);
      test_open(ifs,datapath,dbfile);
      config["DBFILE"] = dbfile;
    
      // set master table and query on it
      tablename = boost::get<std::string>(config["TABLE"]);
      query = "SELECT * FROM " + tablename ;
      where = "1";
      try {
	where = boost::get<std::string>(config["SELECTION"]);
	query += " WHERE " + where;
      } catch (std::invalid_argument) {}
    } catch (std::invalid_argument) {} 

    // store redshifts
    std::vector<double> vd = boost::get<std::vector<double> >(config["REDSHIFT"]);
    for (int i=0; i < vd.size(); i++)
      layers[vd[i]] = shapelens::SourceModelList();

    // set values in ImagingReference
    imref.pixsize = boost::get<double>(config["PIXSIZE"]);
    imref.fov = boost::get<double>(config["FOV"]);

    // define the band set for the imaging reference
    std::vector<std::string> vs, chunks;
    try {
      Band b;
      vs = boost::get<std::vector<std::string> >(config["FILTERS"]);
       // associate each filter with the table names for each model
      for (int i=0; i < vs.size(); i++) {
	// chunks[0] = identifier, chunks[1] = filename
	chunks = split(vs[i],':');
	b.name = chunks[0];
	test_open(ifs,datapath,chunks[1]);
	b.curve = astro::filter(chunks[1],"/");
	imref.bands.insert(b);
      }

      // get the sed for the sources
      vs = boost::get<std::vector<std::string> >(config["SEDS"]);
      for (int i=0; i < vs.size(); i++) {
	// chunks[0] = identifier, chunks[1] = filename
	chunks = split(vs[i],':');
	test_open(ifs,datapath,chunks[1]);
	imref.seds[chunks[0]] = astro::sed(chunks[1],"/");
      }
    } catch (std::invalid_argument) {} 

    // get the tables for all bands of all model types
    std::ostringstream os;
    for (int t = 0; t < 8; t++) { // the numbers of allowed model_type
      os.str("");
      os << "MODEL" << t << "_DBTABLES";
      try {
	vs = boost::get<std::vector<std::string> > (config[os.str()]);
	model_band_tables[t] = std::map<std::string, std::string>();
	for (int i = 0; i < vs.size(); i++) {
	  // chunks[0] = band name, chunks[1] = table name
	  chunks = split(vs[i],':');
	  model_band_tables[t][chunks[0]] = chunks[1];
	}
      } catch (std::invalid_argument) {}
    }
  }

  void SourceCatalog::adjustNumber(const shapelens::Point<double>& fov) {
    unsigned int N_ref = SourceCatalog::size();
    unsigned int N_out = (unsigned int) floor(N_ref * fov(0)*fov(1) / imref.fov);
    replication_ratio = double(N_out)/N_ref;
    RNG& rng = shapelens::Singleton<RNG>::getInstance();
    const gsl_rng* r = rng.getRNG();
    SourceCatalog::iterator iter;
    // too many objects in catalog:
    // remove objects randomly from cat
    if (N_out < N_ref) {
      int selected;
      while (SourceCatalog::size() > N_out) {
	selected = (int) floor(SourceCatalog::size()*gsl_rng_uniform(r)) - 1;
	iter = SourceCatalog::begin();
	if (selected > 0)
	  advance(iter,selected);
	SourceCatalog::erase(iter);
      }
    }
    // too few objects in catalog:
    // replicate objects from catalog (without change of shapelet_models)
    else if (N_out > N_ref) {
      int selected;
      while (SourceCatalog::size() < N_out) {
	// select only from the reference objects
	selected = (int) floor(N_ref*gsl_rng_uniform(r)) - 1;
	iter = SourceCatalog::begin();
	if (selected > 0)
	  advance(iter,selected);
	// duplicate: new entry = CatInfo[selected]
	SourceCatalog::push_back(*iter);
      }
    }
  }

  double SourceCatalog::getReplicationRatio() const {
    return replication_ratio;
  }

  void SourceCatalog::distribute(const shapelens::Point<double>& fov, bool keepPosition) {
    RNG& rng = shapelens::Singleton<RNG>::getInstance();
    const gsl_rng* r = rng.getRNG();
    for (SourceCatalog::iterator iter = SourceCatalog::begin(); iter != SourceCatalog::end(); iter++) {
      iter->redshift_layer = getRedshiftNearestLayer(iter->redshift);
      if (!keepPosition) {
	iter->centroid(0) = fov(0)*gsl_rng_uniform(r);
	iter->centroid(1) = fov(1)*gsl_rng_uniform(r);
	iter->rotation = 2*M_PI*gsl_rng_uniform(r);
	iter->rotation *= GSL_SIGN(-1 + 2*gsl_rng_uniform(r));
      }
    }
  }

  double SourceCatalog::getRedshiftNearestLayer(double z) {
    astro::cosmology& cosmo = SingleCosmology::getInstance();
    double min_dist;
    std::map<double, shapelens::SourceModelList>::iterator iter;
    for (iter = layers.begin(); iter != layers.end(); iter++) {
      if (iter == layers.begin())
	min_dist = fabs(cosmo.properDist(iter->first,z));
      else {
	double dist = fabs(cosmo.properDist(iter->first,z));
	if (dist < min_dist)
	  min_dist = dist;
	else // since redshifts are ordered by redshift, we can stop 
	  break;
      }
    }
    // since we have iterated once too much, decrement iterator
    // to get modellist with closest distance
    iter--;
    return iter->first;
  }

  void SourceCatalog::computeADUinBands(const Telescope& tel, const astro::filter& transmittance) {
    for (SourceCatalog::iterator iter = SourceCatalog::begin(); iter != SourceCatalog::end(); iter++) {
      GalaxyInfo& info = *iter;
      astro::sed s = imref.seds.find(info.sed)->second;
      s.shift(info.redshift);

      // need to compute sed normalization: 
      // use average of measured flux (from info.mags) / sed_flux in same band
      if (info.sed_norm == 0) {
	NumVector<double>norms(info.mags.size());
	unsigned int i = 0;
	double weight = 0, flux, flux_, flux_error;
	for (std::map<std::string,std::pair<double,double> >::const_iterator miter = info.mags.begin(); miter != info.mags.end(); miter++) {
	  for (std::set<Band>::iterator siter = imref.bands.begin(); siter != imref.bands.end(); siter++) {
	    // name of imref band is the one of the magnitude measurements
	    if (siter->name == miter->first) {
	      astro::sed s_ = s;
	      s_ *= siter->curve;
	      // flux in filter with filter Qe correction
	      flux_ = s_.getNorm() / siter->curve.getQe();
	      // flux from magnitude
	      flux = Conversion::mag2flux(miter->second.first);        
	      norms(i) = (flux/flux_); // ignore errors here
	      i++;
	    }
	  }
	}

	// outlier rejection
	std::set<unsigned int> usefull;
	for (unsigned int i=0; i < norms.size(); i++)
	  usefull.insert(i);
	NumVector<double> dev(norms.size());
	double max_dev;
	unsigned int max_dev_item;
	
	do {
	  // build mean and std of norms without o
	  for (std::set<unsigned int>::iterator uiter = usefull.begin();
	       uiter != usefull.end(); uiter++) {
	    unsigned int o = *uiter;
	    double mean_ = 0, tmp;
	    int j = 0;
	    for (int i=0; i< norms.size(); i++) {
	      if (i != o && usefull.find(i) != usefull.end()) {
		tmp = mean_;
		mean_ += norms(i);
		if (j>=1) {
		  dev(o) = (1-1./j)*dev(o) + (j+1)*gsl_pow_2(mean_/(j+1) - tmp/j);
		}
		j++;
	      }
	    }
	    mean_ /= j;
	    double std_ = sqrt(dev(o));
	    if (j >= 2) {// sensible variance
	      dev(o) = fabs(norms(o) - mean_)/sqrt(dev(o));
	      // find position of maximum deviation
	      if (uiter == usefull.begin()) {
		max_dev = dev(o); // initialize
		max_dev_item = o;
	      }
	      else if (dev(o) > max_dev) {
		max_dev = dev(o);
		max_dev_item = o;
	      }
	    }
	  }

	  // if maximum outlier it's larger than 10 (sigma), discard it
	  // and repeat outlier detection without it
	  if (max_dev > 10)
	    usefull.erase(usefull.find(max_dev_item));
	} while (usefull.size() > 2 && max_dev > 10);
	
	  
	// build mean from usefull
	info.sed_norm = 0;
	for (std::set<unsigned int>::iterator uiter = usefull.begin();
	     uiter != usefull.end(); uiter++) {
	  info.sed_norm += norms(*uiter);
	}
	info.sed_norm /= usefull.size();
      }

      if (info.sed_norm != 0) {
	s *= transmittance;
	double flux = info.sed_norm * s.getNorm() / transmittance.getQe();
	info.mag = Conversion::flux2mag(flux);
	double total = 0;

	// only one band: all flux for this model
	// call it "tel" band
	if (model_band_tables[info.model_type].size() == 1)
	  info.adus["tel"] = total = 1;
	else {
	  // compute the overlap of the available bands with
	  // transmittance
	  std::map<std::string, std::string>& band_tables = 
	    model_band_tables[info.model_type];
	  for (std::map<std::string, std::string>::iterator biter = band_tables.begin(); biter != band_tables.end(); biter++) {
	    astro::sed s_ = s;
	    std::set<Band>::iterator band = imref.bands.begin();
	    while (band->name != biter->first && band != imref.bands.end())
	      band++;
	    s_ *= band->curve;
	    total += s_.getNorm();
	    info.adus[band->name] = s_.getNorm();
	  }

	  // eliminate models for bands with less than 10% of total flux
	  double total_red = total;
	  for (std::map<std::string, double>::iterator aiter = info.adus.begin(); aiter != info.adus.end(); aiter++) {
	    if (aiter->second < 0.1*total) {
	      total_red -= aiter->second;
	      info.adus.erase(aiter->first);
	    }
	  }
	  total = total_red;
	}

	// for each entry in info.adus: multiply with flux and devide by total
	// then do conversion to photons and ADUS
	for (std::map<std::string, double>::iterator aiter = info.adus.begin(); aiter != info.adus.end(); aiter++) {
	  aiter->second *= flux/total;
	  aiter->second = Conversion::flux2photons(aiter->second,1,tel,transmittance);
	  aiter->second = Conversion::photons2ADU(aiter->second,tel.gain);
	  // account for change of pixels size: conservation of surface brightness
	  aiter->second /= gsl_pow_2(imref.pixsize);
	  aiter->second *= gsl_pow_2(tel.pixsize);
	}
      } 
    }
  }

  void SourceCatalog::save(shapelens::SQLiteDB& db, int i) const {
    std::ostringstream tablename;
    tablename << "sources_" << i;
    std::string query = "DROP TABLE IF EXISTS " + tablename.str() +";";
    db.query(query);
    query =  "CREATE TABLE " + tablename.str() +" (";
    query += "object_id int NOT NULL,";
    query += "redshift double NOT NULL,";
    query += "sed text NOT NULL,";
    query += "sed_norm double NOT NULL,";
    query += "mag double NOT NULL,";
    query += "centroid_x double NOT NULL,";
    query += "centroid_y double NOT NULL,";
    query += "rotation double NOT NULL,";
    query += "redshift_layer double NOT NULL,";
    query += "model_type int NOT NULL,";
    query += "band text NOT NULL,";
    query += "adu double NOT NULL)";
    db.query(query);

    sqlite3_stmt *stmt;
    query = "INSERT INTO " + tablename.str() +" VALUES (?,?,?,?,?,?,?,?,?,?,?,?);";
    db.checkRC(sqlite3_prepare_v2(db.db, query.c_str(), query.size(), &stmt, NULL));
    for (SourceCatalog::const_iterator iter = SourceCatalog::begin(); iter != SourceCatalog::end(); iter++) {
      for (std::map<std::string, double>::const_iterator aiter = iter->adus.begin(); aiter != iter->adus.end(); aiter++) {
	db.checkRC(sqlite3_bind_int(stmt,1,iter->object_id));
	db.checkRC(sqlite3_bind_double(stmt,2,iter->redshift));
	db.checkRC(sqlite3_bind_text(stmt,3,iter->sed.c_str(),iter->sed.size(),SQLITE_STATIC));
	db.checkRC(sqlite3_bind_double(stmt,4,iter->sed_norm));
	db.checkRC(sqlite3_bind_double(stmt,5,iter->mag));
	db.checkRC(sqlite3_bind_double(stmt,6,iter->centroid(0)));
	db.checkRC(sqlite3_bind_double(stmt,7,iter->centroid(1)));
	db.checkRC(sqlite3_bind_double(stmt,8,iter->rotation));
	db.checkRC(sqlite3_bind_double(stmt,9,iter->redshift_layer));
	db.checkRC(sqlite3_bind_int(stmt,10,iter->model_type));
	db.checkRC(sqlite3_bind_text(stmt,11,aiter->first.c_str(),aiter->first.size(),SQLITE_STATIC));
	db.checkRC(sqlite3_bind_double(stmt,12,aiter->second));
	if(sqlite3_step(stmt)!=SQLITE_DONE)
	  throw std::runtime_error("SourceCatalog: insertion failed: " + std::string(sqlite3_errmsg(db.db)));
	db.checkRC(sqlite3_reset(stmt));
      }
    }
    db.checkRC(sqlite3_finalize(stmt));

    // Insert config in table config
    query  = "CREATE TABLE IF NOT EXISTS config (";
    query += "name TEXT PRIMARY KEY,";
    query += "config TEXT NOT NULL);";
    db.query(query);
    query = "INSERT OR REPLACE INTO config VALUES(?,?);";
    db.checkRC(sqlite3_prepare_v2(db.db, query.c_str(), query.size(), &stmt, NULL));
    db.checkRC(sqlite3_bind_text(stmt,1,tablename.str().c_str(),tablename.str().size(),SQLITE_TRANSIENT));
    tablename.str("");
    config.write(tablename); // tablename now config!!!
    db.checkRC(sqlite3_bind_text(stmt,2,tablename.str().c_str(),tablename.str().size(),SQLITE_TRANSIENT));
    if(sqlite3_step(stmt)!=SQLITE_DONE)
      throw std::runtime_error("SourceCatalog: insertion failed: " + std::string(sqlite3_errmsg(db.db)));
    db.checkRC(sqlite3_finalize(stmt));
  }

  void SourceCatalog::setRotationMatrix(NumMatrix<double>& O, double rotation) const {
    double cos_phi = cos(fabs(rotation)), sin_phi = sin(fabs(rotation));
    // rotation matrix
    O(0,0) = cos_phi;
    O(0,1) = -sin_phi;
    O(1,0) = sin_phi;
    O(1,1) = cos_phi;
    // sign of roation indicates rotation or reflection
    if (GSL_SIGN(rotation) == -1) {
      O(0,1) *= -1;
      O(1,1) *= -1;
    }
  }



  sqlite3_stmt* setPreparedStmt(shapelens::SQLiteDB* db, std::string query) {
    sqlite3_stmt *stmt;
    sqlite3_prepare_v2(db->db, query.c_str(), query.size(), &stmt, NULL);
    return stmt;
  }

  class DBSTMT {
  public:
    shapelens::SQLiteDB* db;
    sqlite3_stmt *stmt;
    bool delete_db;
  };


  void SourceCatalog::createGalaxyLayers(double exptime) {

    // check for all needed model types
    // open the DB files
    // and prepare the fetch statements for each band
    // stmts: type -> band -> (db,stmt)
    std::map<char, std::map<std::string, DBSTMT > > stmts;
    std::ostringstream os;
    std::string datapath = getDatapath(), dbfile, query;
    std::vector<std::string> tables, chunks;
    std::ifstream ifs;
    for (std::set<char>::iterator iter = need_model.begin();
	 iter != need_model.end(); iter++) {
      os.str("");
      os << "MODEL" << int(*iter) << "_DBFILE";
      std::map<std::string, DBSTMT > band_dbs;
      DBSTMT dbstmt;
      // throws if not present....
      std::string dbfile = boost::get<std::string>(config[os.str()]);
      if (dbfile == "$APPLICATION_DB$") {
	shapelens::SQLiteDB& adb = shapelens::Singleton<shapelens::SQLiteDB>::getInstance();
	dbstmt.db = &adb;
	dbstmt.delete_db = false;
      }
      else {
	/// throws if not found
	test_open(ifs,datapath,dbfile);
	dbstmt.db = new shapelens::SQLiteDB();
	dbstmt.db->connect(dbfile);
	dbstmt.delete_db = true;
      }
      
      // prepare statements for all tables
      std::map<std::string, std::string>& band_tables = model_band_tables[*iter];
      for (std::map<std::string, std::string>::iterator biter = band_tables.begin(); biter != band_tables.end(); biter++) {
      	// define statement for different model types
	switch (*iter) {
	case 0: // sersic type
	  query = "SELECT n_sersic, radius, ellipticity FROM " + biter->second + " WHERE id = ?;";
	  break;
	case 1: // shapelet
	  query = "SELECT flags,beta,nmax,coeffs,min_x,min_y,size_x,size_y,centroid_x,centroid_y FROM " + biter->second + " WHERE id = ?;";
	  break;
	case 2: // bulge-disk
	  query = "SELECT bulge_total, bulge_ns, bulge_Re, bulge_eps, disk_Re, disk_eps FROM " + biter->second + " WHERE id = ?;";
	}
	dbstmt.stmt = setPreparedStmt(dbstmt.db,query);
	band_dbs[biter->first] = dbstmt;
      }
      stmts[*iter] = band_dbs;
    }

    // iterator through SourceCatalog
    // query dbs for each element (in each band)
    // create the appropriate model from it then an abstract SourceModel
    NumMatrix<double> O(2,2);
    for (SourceCatalog::const_iterator iter = SourceCatalog::begin(); iter != SourceCatalog::end(); iter++) {
      const GalaxyInfo& info = *iter;

      // get correct layer from info.redshift_layer
      // or create it if it doesn't exist
      std::map<double, shapelens::SourceModelList>::iterator liter = layers.find(info.redshift_layer);
      if (liter == layers.end())
	liter = layers.insert(layers.begin(),std::pair<double,shapelens::SourceModelList>(info.redshift_layer,shapelens::SourceModelList()));
	
      // set rotation/parity flip matrix
      setRotationMatrix(O, info.rotation);
      // account for original pixel size
      O *= imref.pixsize;
      // apply linear transformation from O
      shapelens::LinearTransformation A(O);
      // and shift the centroid
      shapelens::ShiftTransformation Z(info.centroid);
      A *= Z;
      
      // get the pair (db,stmt) from model_type and band
      for (std::map<std::string, double>::const_iterator aiter = info.adus.begin(); aiter != info.adus.end(); aiter++) {
	DBSTMT& dbstmt = stmts[info.model_type][aiter->first];
	double flux = aiter->second * exptime;
	double missing_flux = 0;
	// bind stmt to object id
	dbstmt.db->checkRC(sqlite3_bind_int(dbstmt.stmt,1,info.object_id));
	if (sqlite3_step(dbstmt.stmt) == SQLITE_ROW) {
	  // get data out of stmt and populate a model

	  if (info.model_type == 0) { // Sersic model
	    double n_sersic = sqlite3_column_double(dbstmt.stmt,0);
	    double radius = sqlite3_column_double(dbstmt.stmt,1);
	    double e = sqlite3_column_double(dbstmt.stmt,2);
	    double b_a = 1 - e;
	    double epsilon = e/(1+b_a);
	    std::complex<double> eps(epsilon,0); // random orientation via A
	    // add missing_flux from invalid models
	    if (missing_flux > 0) {
	      flux += missing_flux;
	      missing_flux = 0;
	    }
	    liter->second.push_back(boost::shared_ptr<shapelens::SourceModel>(new shapelens::SersicModel(n_sersic, radius, flux , eps, 5, &A, info.object_id)));
	  }

	  else if (info.model_type == 1) { // Shapelet models
	    std::bitset<16> flags(sqlite3_column_int(dbstmt.stmt,0));
	    if (!flags.test(15)) { // only use valid models
	      double beta = sqlite3_column_double(dbstmt.stmt,1);
	      int nmax = sqlite3_column_int(dbstmt.stmt,2);
	      shapelens::CoefficientVector<shapelens::data_t> coeffs(nmax);
	      memcpy(coeffs.c_array(), reinterpret_cast<const shapelens::data_t*>(sqlite3_column_blob(dbstmt.stmt,3)), sqlite3_column_bytes(dbstmt.stmt,3));
	      shapelens::Grid grid(sqlite3_column_int(dbstmt.stmt,4),
				   sqlite3_column_int(dbstmt.stmt,5),
				   sqlite3_column_int(dbstmt.stmt,6),
				   sqlite3_column_int(dbstmt.stmt,7));
	      shapelens::Point<shapelens::data_t> centroid(sqlite3_column_double(dbstmt.stmt,8),
							   sqlite3_column_double(dbstmt.stmt,9));
	      // FIXME: if replication_ratio > 1, it would be usefull
	      // to store sobjs because they will be reused.
	      boost::shared_ptr<shapelens::ShapeletObject> sobj(new shapelens::ShapeletObject(coeffs,beta,centroid,grid));
	      sobj->setID(info.object_id);
	      // add missing_flux from invalid models
	      if (missing_flux > 0) {
		flux += missing_flux;
		missing_flux = 0;
	      }
	      liter->second.push_back(boost::shared_ptr<shapelens::SourceModel>(new shapelens::ShapeletModel(sobj, flux ,&A)));
	    } else {
	      missing_flux += flux;
	    }
	  }

	  else if (info.model_type == 2) { // bulge-disk
	    // add missing_flux from invalid models
	    if (missing_flux > 0) {
	      flux += missing_flux;
	      missing_flux = 0;
	    }
	    // bulge-total luminosity ratio
	    double b_t = fabs(sqlite3_column_double(dbstmt.stmt,0));
	    // bulge component
	    double n_sersic = sqlite3_column_double(dbstmt.stmt,1);
	    double radius = sqlite3_column_double(dbstmt.stmt,2);
	    double e = sqlite3_column_double(dbstmt.stmt,3);
	    double b_a = 1 - e;
	    double epsilon = e/(1+b_a);
	    std::complex<double> eps(epsilon,0); // random orientation via A
	    liter->second.push_back(boost::shared_ptr<shapelens::SourceModel>(new shapelens::SersicModel(n_sersic, radius, b_t*flux , eps, 5, &A, info.object_id)));
	    // disk component available?
	    if (sqlite3_column_type(dbstmt.stmt,4) != SQLITE_NULL) {
	      n_sersic = 1;
	      radius = sqlite3_column_double(dbstmt.stmt,4);
	      e = sqlite3_column_double(dbstmt.stmt,5);
	      b_a = 1 - e;
	      epsilon = e/(1+b_a);
	      real(eps) = epsilon; 
	      liter->second.push_back(boost::shared_ptr<shapelens::SourceModel>(new shapelens::SersicModel(n_sersic, radius, (1-b_t)*flux , eps, 5, &A, info.object_id)));
	    }

	  } else {
	    os.str("");
	    os << "SourceCatalog: no information available for model_type" << int(info.model_type) << "!";
	    throw std::runtime_error(os.str());
	  }

	} // end: sqlite3_step == SQLITE_ROW

	else { // no data for model in DB
	  os.str("");
	  os << "SourceCatalog: no entry found for model_type" << int(info.model_type) << " and object " << info.object_id << "!";
	  throw std::runtime_error(os.str());
	}
	// reset statement for usage with new id
	dbstmt.db->checkRC(sqlite3_reset(dbstmt.stmt));
      }
     }

    // finalize stmts and close DBs
    for (std::map<char, std::map<std::string, DBSTMT > >::iterator iter = stmts.begin(); iter != stmts.end(); iter++) {
      for (std::map<std::string, DBSTMT >::iterator diter = iter->second.begin(); diter != iter->second.end(); diter++) {
	DBSTMT& dbstmt = diter->second;
	sqlite3_finalize(dbstmt.stmt);
	if (diter == --(iter->second.end()))
	  if (dbstmt.delete_db)
	    delete dbstmt.db;
      }
    }

    // create GalaxyLayer for each SourceModelList in layers
    for (std::map<double, shapelens::SourceModelList>::iterator liter = layers.begin(); liter != layers.end(); liter++)
      new GalaxyLayer(liter->first,liter->second);
  }


} // end namespace
