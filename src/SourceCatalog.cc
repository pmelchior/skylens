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
	  info.object_id = boost::lexical_cast<unsigned long>(row[0]);
	else if (name.substr(0,3) == "mag") {
	  for (std::set<Band>::iterator iter = imref.bands.begin(); iter != imref.bands.end(); iter++) {
	    if (name == "mag_" + iter->name) {
	      if (row[i] != NULL) {
		info.mags[iter->name] = 
		  std::pair<double,double> (boost::lexical_cast<double>(row[i]),
					    boost::lexical_cast<double>(row[i+1]));
		i++; // skip next field: error of mag in current band
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
	else if (name == "radius") {
	  if (row[i] != NULL)
	    info.radius = boost::lexical_cast<double>(row[i]);
	  else
	    info.radius = 0;
	}
	else if (name == "ellipticity") {
	  if (row[i] != NULL)
	    info.ellipticity = boost::lexical_cast<double>(row[i]);
	  else
	    info.ellipticity = 0;
	}
	else if (name == "n_sersic") {
	  if (row[i] != NULL)
	    info.n_sersic = boost::lexical_cast<double>(row[i]);
	  else
	    info.n_sersic = 0;
	}
	else if (name == "model_type") {
	  info.model_type = boost::lexical_cast<double>(row[i]);
	  if (info.model_type > 0) {
	    if (imref.bands.begin()->dbdetails.find(info.model_type) == imref.bands.begin()->dbdetails.end()) {
	      std::ostringstream out;
	      out << "SourceCatalog: No information for model_type " << int(info.model_type) << "!";
	      throw std::invalid_argument(out.str());
	    }
	  }
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
    } else
      throw std::runtime_error("SourceCatalog: no config found DB for table " + tablename.str());
	
    // parse contents of config file
    parseConfig();

    // get additional infos from db
    query = "SELECT band_overlap,replication_ratio FROM source_info WHERE name = '" + tablename.str() + "';";
    db.checkRC(sqlite3_finalize(stmt));
    db.checkRC(sqlite3_prepare_v2(db.db, query.c_str(), query.size(), &stmt, NULL));
    if (sqlite3_step(stmt) == SQLITE_ROW) {
      std::string overlap(reinterpret_cast<const char*>(sqlite3_column_text(stmt,0)));
      std::vector<std::string> chunks = split(overlap,',');
      for (int i=0; i < chunks.size(); i++) {
	std::vector<std::string> bchunks = split(chunks[i],':');
	for (std::set<Band>::iterator biter = imref.bands.begin(); biter != imref.bands.end(); biter++)
	  if (biter->name == bchunks[0])
	    const_cast<Band&>(*biter).overlap = boost::lexical_cast<double>(bchunks[1]);
      }
      replication_ratio = sqlite3_column_double(stmt,1);
    } else
      throw std::runtime_error("SourceCatalog: no infos found DB for table " + tablename.str());

    // get contents from each source table
    GalaxyInfo info;
    query = "SELECT * FROM " + tablename.str() +";";
    db.checkRC(sqlite3_finalize(stmt));
    db.checkRC(sqlite3_prepare_v2(db.db, query.c_str(), query.size(), &stmt, NULL));
    while (sqlite3_step(stmt) == SQLITE_ROW) {
      info.object_id = sqlite3_column_int(stmt,0);
      info.redshift = sqlite3_column_double(stmt,1);
      info.sed = std::string(reinterpret_cast<const char*>(sqlite3_column_text(stmt,2)));
      info.radius = sqlite3_column_double(stmt,3);
      info.ellipticity = sqlite3_column_double(stmt,4);
      info.n_sersic = sqlite3_column_double(stmt,5);
      info.sed_norm = sqlite3_column_double(stmt,6);
      info.mag = sqlite3_column_double(stmt,7);
      info.centroid(0) = sqlite3_column_double(stmt,8);
      info.centroid(1) = sqlite3_column_double(stmt,9);
      info.rotation = sqlite3_column_double(stmt,10);
      info.redshift_layer = sqlite3_column_double(stmt,11);
      info.model_type = sqlite3_column_int(stmt,12);
      info.adus.clear();
      std::string band(reinterpret_cast<const char*>(sqlite3_column_text(stmt,13)));
      info.adus[band] = sqlite3_column_double(stmt,14);
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

    // store redshifts
    std::vector<double> vd = boost::get<std::vector<double> >(config["REDSHIFT"]);
    for (int i=0; i < vd.size(); i++)
      layers[vd[i]] = shapelens::SourceModelList();

    // set values in ImagingReference
    imref.pixsize = boost::get<double>(config["PIXSIZE"]);
    imref.fov = boost::get<double>(config["FOV"]);
    // define the band set
    Band b;
    std::vector<std::string> vs = boost::get<std::vector<std::string> >(config["FILTERS"]);
    // get table names for all model types
    std::map<char, std::vector<std::string> > modeltables;
    std::ostringstream os;
    for (char m=1; m < 8; m++) {
      os.str("");
      os << "MODEL" << int(m) << "_DBTABLES";
      try {
	std::vector<std::string> filtertables = boost::get<std::vector<std::string> >(config[os.str()]);
	if (filtertables.size() != vs.size())
	  throw std::runtime_error("SourceCatalog: " + os.str() + " has wrong number of entries!");
	modeltables[m] = filtertables;
      } catch (std::invalid_argument) {}
    }

    // associate each filter with the table names for each model
    std::vector<std::string> chunks;
    for (int i=0; i < vs.size(); i++) {
      // chunks[0] = identifier, chunks[1] = filename
      chunks = split(vs[i],':');
      b.name = chunks[0];
      test_open(ifs,datapath,chunks[1]);
      b.curve = filter(chunks[1],"/");
      for (char m=0; m < 8; m++) {
	std::map<char,std::vector<std::string> >::iterator iter = 
	  modeltables.find(m);
	if (iter != modeltables.end()) {
	  b.dbdetails[m] = iter->second[i];
	}
      }
      imref.bands.insert(b);
    }

    // get the sed for the sources
    vs = boost::get<std::vector<std::string> >(config["SEDS"]);
    for (int i=0; i < vs.size(); i++) {
      // chunks[0] = identifier, chunks[1] = filename
      chunks = split(vs[i],':');
      test_open(ifs,datapath,chunks[1]);
      imref.seds[chunks[0]] = sed(chunks[1],"/");
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

  void SourceCatalog::distribute(const shapelens::Point<double>& fov) {
    RNG& rng = shapelens::Singleton<RNG>::getInstance();
    const gsl_rng* r = rng.getRNG();
    for (SourceCatalog::iterator iter = SourceCatalog::begin(); iter != SourceCatalog::end(); iter++) {
      iter->centroid(0) = fov(0)*gsl_rng_uniform(r);
      iter->centroid(1) = fov(1)*gsl_rng_uniform(r);
      iter->rotation = 2*M_PI*gsl_rng_uniform(r);
      iter->rotation *= GSL_SIGN(-1 + 2*gsl_rng_uniform(r));
      iter->redshift_layer = getRedshiftNearestLayer(iter->redshift);
    }
    
  }

  void SourceCatalog::selectOverlapBands(const filter& transmittance, double fraction) {
    // compute overlap with telescopes total transmittance
    double sum = 0;
    std::set<Band>::iterator iter;
    for (iter = imref.bands.begin(); iter != imref.bands.end(); iter++) {
      filter f = iter->curve;
      f *= transmittance;
      // trick to modify element of set:
      // iterator does formally not allow write access to it
      const_cast<Band&>(*iter).overlap = f.getQe();
      sum += f.getQe();
    }

    // to be included band must fraction of overlap
    double limit = sum * fraction;
    sum = 0;
    std::cout << "get overlapping bands" << std::endl;
    iter = imref.bands.begin();
    for (iter = imref.bands.begin(); iter != imref.bands.end(); iter++) {
      std::cout << iter->name << "\t" << iter->overlap << "\t" << limit << std::endl;
      if (iter->overlap < limit) 
	const_cast<Band&>(*iter).overlap = 0;
      else
	sum += iter->overlap;
    }

    for (iter = imref.bands.begin(); iter != imref.bands.end(); iter++)
      const_cast<Band&>(*iter).overlap /= sum;
  }

  
  double SourceCatalog::getRedshiftNearestLayer(double z) {
    cosmology& cosmo = SingleCosmology::getInstance();
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

  void SourceCatalog::computeADUinBands(const Telescope& tel, const filter& transmittance) {
    for (SourceCatalog::iterator iter = SourceCatalog::begin(); iter != SourceCatalog::end(); iter++) {
      GalaxyInfo& info = *iter;
      std::cout << info.object_id << "\t" << info.sed << std::endl;
      sed s = imref.seds.find(info.sed)->second;
      s.shift(info.redshift);

      // need to compute sed normalization: 
      // use average of measured flux (from info.mags) / sed_flux in same band
      if (info.sed_norm == 0) {
	double weight = 0, flux, flux_, flux_error;
	for (std::map<std::string,std::pair<double,double> >::const_iterator iter = info.mags.begin(); iter != info.mags.end(); iter++) {
	  for (std::set<Band>::iterator siter = imref.bands.begin(); siter != imref.bands.end(); siter++) {
	    // name of imref band is the one of the magnitude measurements
	    if (siter->name == iter->first && siter-> name != "B") {
	      sed s_ = s;
	      s_ *= siter->curve;
	      // flux in filter with filter Qe correction
	      flux_ = s_.getNorm() / siter->curve.getQe(); 
	      flux = Conversion::mag2flux(iter->second.first);        // flux from magnitude
	      // FIXME: magnitude errors
	      flux_error = 1;//Conversion::mag2flux(iter->second.second); // weight from magnitude error
	      weight += 1./flux_error;
	      info.sed_norm += (flux/flux_)/flux_error;
	      // FIXME: need outlier rejection: blue band has catastrophic outliers!!!
	      std::cout << info.object_id << "\t" << iter->second.first << "\t" << flux << "\t" << flux_ << "\t" << flux/flux_ << std::endl;
	    }
	  }
	}
	info.sed_norm /= weight;
      }

      if (info.sed_norm != 0) {
	s *= transmittance;
	double flux = info.sed_norm * s.getNorm() / transmittance.getQe();
	info.mag = Conversion::flux2mag(flux);
	std::cout << "\t" << info.mag << "\t" << flux << std::endl;
	double total = 0;
	if (info.model_type > 0) {
	  // find overlap of the filtered s with imref's bands to get their relative weights
	  for (std::set<Band>::iterator siter = imref.bands.begin(); siter != imref.bands.end(); siter++) {
	    // only consider bands with non-vanishing overlap to transmittance
	    if (siter->overlap > 0) {
	      sed s_ = s;
	      s_ *= siter->curve;
	      total += s_.getNorm();
	      info.adus[siter->name] = s_.getNorm();
	    }
	  }
	} else { // only one band
	  info.adus["tel"] = total = 1;
	}

	// for each entry in info.adus: multiply with flux and devide by total
	// then do conversion to photons and ADUS
	for (std::map<std::string, double>::iterator iter = info.adus.begin(); iter != info.adus.end(); iter++) {
	  iter->second *= flux/total;
	  iter->second = Conversion::flux2photons(iter->second,1,tel,transmittance);
	  iter->second = Conversion::photons2ADU(iter->second,tel.gain);
	  // account for change of pixels size: conservation of surface brightness
	  iter->second /= gsl_pow_2(imref.pixsize);
	  iter->second *= gsl_pow_2(tel.pixsize);
	}
	std::cout << std::endl;
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
    query += "radius double NOT NULL,";
    query += "ellipticity double NOT NULL,";
    query += "n_sersic double NOT NULL,";
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
    query = "INSERT INTO " + tablename.str() +" VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?);";
    db.checkRC(sqlite3_prepare_v2(db.db, query.c_str(), query.size(), &stmt, NULL));
    for (SourceCatalog::const_iterator iter = SourceCatalog::begin(); iter != SourceCatalog::end(); iter++) {
      for (std::map<std::string, double>::const_iterator aiter = iter->adus.begin(); aiter != iter->adus.end(); aiter++) {
	db.checkRC(sqlite3_bind_int(stmt,1,iter->object_id));
	db.checkRC(sqlite3_bind_double(stmt,2,iter->redshift));
	db.checkRC(sqlite3_bind_text(stmt,3,iter->sed.c_str(),iter->sed.size(),SQLITE_STATIC));
	db.checkRC(sqlite3_bind_double(stmt,4,iter->radius));
	db.checkRC(sqlite3_bind_double(stmt,5,iter->ellipticity));
	db.checkRC(sqlite3_bind_double(stmt,6,iter->n_sersic));
	db.checkRC(sqlite3_bind_double(stmt,7,iter->sed_norm));
	db.checkRC(sqlite3_bind_double(stmt,8,iter->mag));
	db.checkRC(sqlite3_bind_double(stmt,9,iter->centroid(0)));
	db.checkRC(sqlite3_bind_double(stmt,10,iter->centroid(1)));
	db.checkRC(sqlite3_bind_double(stmt,11,iter->rotation));
	db.checkRC(sqlite3_bind_double(stmt,12,iter->redshift_layer));
	db.checkRC(sqlite3_bind_int(stmt,13,iter->model_type));
	db.checkRC(sqlite3_bind_text(stmt,14,aiter->first.c_str(),aiter->first.size(),SQLITE_STATIC));
	db.checkRC(sqlite3_bind_double(stmt,15,aiter->second));
	if(sqlite3_step(stmt)!=SQLITE_DONE)
	  throw std::runtime_error("SourceCatalog: insertion failed: " + std::string(sqlite3_errmsg(db.db)));
	db.checkRC(sqlite3_reset(stmt));
      }
    }
    db.checkRC(sqlite3_finalize(stmt));

    // Insert additional infos in table source_info
    query  = "CREATE TABLE IF NOT EXISTS source_info (";
    query += "name TEXT PRIMARY KEY,";
    query += "band_overlap TEXT NOT NULL,";
    query += "replication_ratio DOUBLE NOT NULL);";
    db.query(query);
    query = "INSERT OR REPLACE INTO source_info VALUES(?,?,?);";
    db.checkRC(sqlite3_prepare_v2(db.db, query.c_str(), query.size(), &stmt, NULL));
    db.checkRC(sqlite3_bind_text(stmt,1,tablename.str().c_str(),tablename.str().size(),SQLITE_TRANSIENT));
    std::ostringstream line;
    for (std::set<Band>::iterator biter = imref.bands.begin(); biter != imref.bands.end(); biter++) {
      line << biter->name << ":" << biter->overlap;
      if (biter != --(imref.bands.end()))
	line << ",";
    }
    db.checkRC(sqlite3_bind_text(stmt,2,line.str().c_str(),line.str().size(),SQLITE_TRANSIENT));
    db.checkRC(sqlite3_bind_double(stmt,3,replication_ratio));
    if(sqlite3_step(stmt)!=SQLITE_DONE)
      throw std::runtime_error("SourceCatalog: insertion failed: " + std::string(sqlite3_errmsg(db.db)));
    db.checkRC(sqlite3_finalize(stmt));

    // Insert config in table config
    query  = "CREATE TABLE IF NOT EXISTS config (";
    query += "name TEXT PRIMARY KEY,";
    query += "config TEXT NOT NULL);";
    db.query(query);
    query = "INSERT OR REPLACE INTO config VALUES(?,?);";
    db.checkRC(sqlite3_prepare_v2(db.db, query.c_str(), query.size(), &stmt, NULL));
    db.checkRC(sqlite3_bind_text(stmt,1,tablename.str().c_str(),tablename.str().size(),SQLITE_TRANSIENT));
    line.str("");
    config.write(line);
    db.checkRC(sqlite3_bind_text(stmt,2,line.str().c_str(),line.str().size(),SQLITE_TRANSIENT));
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


  void SourceCatalog::createGalaxyLayers(double exptime) {
    // model_type = 0: use catalog infos to create SersicModels
    // model_type = 1: query ShapeletDB with a join query 
    //  (select subset of the catalog query with model_type = 1)
    //  use a sorted lookup table to relate entries in catalog with entries in ShapeletObjectList

    // since smodels is temporary, all unneeded shapelet models
    // will be destroyed after smodel goes out of scope
    std::map<std::string, ShapeletObjectCat> smodels = getShapeletModels();

    // FIXME: test if shapelet models are required
    //if (imref.bands.begin()->dbdetails.find(char(1)) != imref.bands.begin()->dbdetails.end())
    //  smodels = getShapeletModels();

    NumMatrix<double> O(2,2);
    for (SourceCatalog::const_iterator iter = SourceCatalog::begin(); iter != SourceCatalog::end(); iter++) {
      const GalaxyInfo& info = *iter;

      // set rotation/parity flip matrix
      setRotationMatrix(O, info.rotation);
      // account for original pixel size
      O *= imref.pixsize;
      // apply linear transformation from O
      shapelens::LinearTransformation A(O);
      // and shift the centroid
      shapelens::ShiftTransformation Z(info.centroid);
      A *= Z;
      
      double missing_flux = 0;
      // model switch
      if (info.model_type == 1) {
	// create model for every band in info.adus
	for (std::map<std::string, double>::const_iterator aiter = info.adus.begin(); aiter != info.adus.end(); aiter++) {
	  // selects models for given band
	  ShapeletObjectCat& sl = smodels[aiter->first];
	  // weigh model with its relative flux
	  double flux = aiter->second * exptime;
	  // push it onto correct layer
	  ShapeletObjectCat::iterator siter = sl.find(info.object_id);
	  if (siter != sl.end()) {
	    // if it has a valid model:
	    if (!(siter->second->getFlags().test(15)))
	      layers[info.redshift_layer].push_back(boost::shared_ptr<shapelens::SourceModel>(new shapelens::ShapeletModel(siter->second, flux, &A)));
	    else // if not: remember flux and create Sersic model with it
	      missing_flux += flux;
	  } else {
	    std::ostringstream error;
	    error << "SourceCatalog: no shapelet model for object" << info.object_id; 
	    throw std::runtime_error(error.str());
	  }
	}
      }
      else if (info.model_type == 0 || missing_flux) {
	double e = info.ellipticity;
	double b_a = 1 - e;
	double epsilon = e/(1+b_a);
	std::complex<double> eps(epsilon,0);                // random orientation via A
	double flux;
	if (missing_flux) // make Sersic model for invalid ones
	  flux = missing_flux;
	else // there is only a single band for model_type = 0
	  flux = info.adus.begin()->second * exptime;
	layers[info.redshift_layer].push_back(boost::shared_ptr<shapelens::SourceModel>(new shapelens::SersicModel(info.n_sersic, info.radius, flux , eps, &A, info.object_id)));
      }

      else // no implementation for model_type > 1 !
	throw std::invalid_argument("SourceCatalog: creation of models with model_type > 1 not implemented");
    }

    // create GalaxyLayer for each SourceModelList in layers
    for (std::map<double, shapelens::SourceModelList>::iterator iter = layers.begin(); iter != layers.end(); iter++)
      new GalaxyLayer(iter->first,iter->second);

  }


  std::map<std::string, ShapeletObjectCat> SourceCatalog::getShapeletModels() {
    std::map<std::string, ShapeletObjectCat> models;

    // get DBFILE for model_type = 1
    std::string dbfile = boost::get<std::string>(config["MODEL1_DBFILE"]);
    shapelens::SQLiteDB db;
    std::ifstream ifs;
    std::string datapath = getDatapath();
    test_open(ifs,datapath,dbfile);
    db.connect(dbfile);
    
    for (std::set<Band>::iterator biter = imref.bands.begin(); biter != imref.bands.end(); biter++) {
      if (biter->overlap > 0) {
	ShapeletObjectCat scat;
	//select correct table from bandname
	std::map<char, std::string>::const_iterator iter =  biter->dbdetails.find(1);
	if (iter != biter->dbdetails.end()) {
	  shapelens::ShapeletObjectList sl(db,iter->second);
	  // insert entries of sl into models
	  for (shapelens::ShapeletObjectList::iterator iter = sl.begin();
	       iter != sl.end(); iter++) {
	    scat[(*iter)->getObjectID()] = *iter;
	  }
	  models[biter->name] = scat;
	} else {
	  throw std::runtime_error("SourceCatalog: table for model_type = 1 and band " + biter->name + " missing!");
	}
      }
    }
    return models;
  }


} // end namespace
