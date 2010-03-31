#include "../include/SourceCatalog.h"
#include "../include/RNG.h"
#include "../include/Layer.h"
#include "../include/Conversion.h"
#include <shapelens/utils/MySQLDB.h>
#include <shapelens/utils/Property.h>
#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>
#include <math.h>

namespace skylens {
  bool SourceCatalog::Band::operator<(const Band& b) const {
    if (curve.lambdaEff() < b.curve.lambdaEff())
      return true;
    else 
      return false;
  }
  
  SourceCatalog::SourceCatalog(std::string configfile) : 
    std::map<unsigned long, GalaxyInfo> () {
    
    // read in contents of config file
    readConfig(configfile);
    // connect to DB
    shapelens::MySQLDB db;
    db.connect("SKYLENSDBCONF");

    // query database
    shapelens::DBResult dbr = db.query(query);
    MYSQL_ROW row;
    MYSQL_FIELD* fields = dbr.getFields();
    GalaxyInfo info;

    // build catalog from DB entries
    while (row = dbr.getRow()) {
      for (int i =0; i < dbr.getFieldCount(); i++) {
	std::string name = fields[i].name;
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
	else if (name == "model_type")
	  info.model_type = boost::lexical_cast<double>(row[i]);
      }
      std::map<unsigned long, GalaxyInfo>::insert(std::pair<unsigned long, GalaxyInfo>(info.object_id,info));
    }
  }

  SourceCatalog::SourceCatalog(std::string configfile, std::string catfile) :
    std::map<unsigned long, GalaxyInfo> () {
    
    // read in contents of config file
    readConfig(configfile);

    // read in contents of catfile and populate SourceCat from it
    GalaxyInfo info;
    std::string line;
    typedef boost::tokenizer<boost::char_separator<char> > Tok;
    std::vector<std::string> chunks;
    boost::char_separator<char> fieldsep("\t");
    boost::char_separator<char> vectorsep(", ");
    std::ifstream ifs (catfile.c_str());
    while(getline(ifs, line)) {
      if (line[0] != '#' && line != "") {
	Tok tok(line, fieldsep);
	Tok::iterator tok_iter = tok.begin();
	info.object_id = boost::lexical_cast<unsigned long>(*tok_iter);
	tok_iter++;
	info.redshift = boost::lexical_cast<double>(*tok_iter);
	tok_iter++;
	info.sed = *tok_iter;
	tok_iter++;
	info.radius = boost::lexical_cast<double>(*tok_iter);
	tok_iter++;
	info.ellipticity = boost::lexical_cast<double>(*tok_iter);
	tok_iter++;
	info.n_sersic = boost::lexical_cast<double>(*tok_iter);
	tok_iter++;
	Tok maptok(*tok_iter, vectorsep);
	info.mags.clear();
	for(Tok::iterator mtok_iter = maptok.begin(); mtok_iter != maptok.end(); ++mtok_iter) {	      
	  // chunks[0]: band name, chunks[1]: mag, chunks[2]: mag error
	  chunks = split(*mtok_iter,':');
	  info.mags[chunks[0]] = 
	    std::pair<double, double>(boost::lexical_cast<double>(chunks[1]),
				      boost::lexical_cast<double>(chunks[2]));
	}
	tok_iter++;
	info.model_type =  boost::lexical_cast<int>(*tok_iter);
	tok_iter++;
	    
	Tok maptok2(*tok_iter, vectorsep);
	info.adus.clear();
	for(Tok::iterator mtok_iter = maptok2.begin(); mtok_iter != maptok2.end(); ++mtok_iter) {	      
	  // chunks[0]: band name, chunks[1]: adu
	  chunks = split(*mtok_iter,':');
	  info.adus[chunks[0]] = boost::lexical_cast<double>(chunks[1]);
	}
	tok_iter++;
	info.mag = boost::lexical_cast<double>(*tok_iter);
	tok_iter++;
	info.centroid(0) = boost::lexical_cast<double>(*tok_iter);
	tok_iter++;
	info.centroid(1) = boost::lexical_cast<double>(*tok_iter);
	tok_iter++;
	info.rotation = boost::lexical_cast<double>(*tok_iter);
	tok_iter++;
	info.redshift_layer = boost::lexical_cast<double>(*tok_iter);
	// insert into catalog
	SourceCatalog::operator[](SourceCatalog::size() + 1) = info;
      }
    }
  }

  void SourceCatalog::readConfig(std::string configfile) { 
    // open source config file
    std::ifstream ifs;
    std::string datapath = getDatapath();
    test_open(ifs,datapath,configfile);
    shapelens::Property config;
    config.read(ifs);
    ifs.close();

    // set master table and query on it
    tablename = boost::get<std::string>(config["DBTABLE"]);
    query = "SELECT * FROM " + tablename ;
    try {
      std::string where = boost::get<std::string>(config["SELECTION"]);
      query += " WHERE " + where;
    } catch (std::invalid_argument) {}

    // store redshifts
    std::vector<double> vd = boost::get<std::vector<double> >(config["REDSHIFT"]);
    for (int i=0; i < vd.size(); i++)
      redshifts.insert(vd[i]);

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

  void SourceCatalog::adjustNumber(const Telescope& tel) {
    unsigned int N_ref = std::map<unsigned long, GalaxyInfo>::size();
    unsigned int N_out = (unsigned int) floor(N_ref * tel.fov_x*tel.fov_y / imref.fov);
    RNG& rng = shapelens::Singleton<RNG>::getInstance();
    const gsl_rng* r = rng.getRNG();
    SourceCatalog::iterator iter;
    // too many objects in catalog:
    // remove objects randomly from cat
    if (N_out < N_ref) {
      int selected;
      while (std::map<unsigned long, GalaxyInfo>::size() > N_out) {
	selected = (int) floor(std::map<unsigned long, GalaxyInfo>::size()*gsl_rng_uniform(r)) - 1;
	iter = std::map<unsigned long, GalaxyInfo>::begin();
	if (selected > 0)
	  advance(iter,selected);
	std::map<unsigned long, GalaxyInfo>::erase(iter);
      }
    }
    // too few objects in catalog:
    // replicate objects from catalog (without change of shapelet_models)
    else if (N_out > N_ref) {
      int selected;
      // get maximum id in cat
      iter = std::map<unsigned long, GalaxyInfo>::end();
      iter--;
      unsigned int max = iter->first;
      while (std::map<unsigned long, GalaxyInfo>::size() < N_out) {
	// select only from the reference objects
	selected = (int) floor(N_ref*gsl_rng_uniform(r)) - 1;
	max++;
	iter = std::map<unsigned long, GalaxyInfo>::begin();
	if (selected > 0)
	  advance(iter,selected);
	// duplicate: cat[max] = CatInfo[selected]
	std::map<unsigned long, GalaxyInfo>::insert(std::pair<unsigned long, GalaxyInfo>(max,iter->second));
      }
    }
  }

  void SourceCatalog::distribute(const Telescope& tel) {
    RNG& rng = shapelens::Singleton<RNG>::getInstance();
    const gsl_rng* r = rng.getRNG();
    for (SourceCatalog::iterator iter = SourceCatalog::begin(); iter != SourceCatalog::end(); iter++) {
      iter->second.centroid(0) = tel.fov_x*gsl_rng_uniform(r);
      iter->second.centroid(1) = tel.fov_y*gsl_rng_uniform(r);
      iter->second.rotation = 2*M_PI*gsl_rng_uniform(r);
      iter->second.rotation *= GSL_SIGN(-1 + 2*gsl_rng_uniform(r));
      iter->second.redshift_layer = getRedshiftNearestLayer(iter->second.redshift);
    }
    
  }

  void SourceCatalog::selectOverlapBands(const Telescope& tel, double fraction) {
    // compute overlap with telescopes total transmittance
    double sum = 0;
    std::set<Band>::iterator iter;
    for (iter = imref.bands.begin(); iter != imref.bands.end(); iter++) {
      filter f = iter->curve;
      f *= tel.total;
      // trick to modify element of set:
      // iterator does formally not allow write access to it
      const_cast<Band&>(*iter).overlap = f.getQe();
      sum += f.getQe();
    }

    // to be included band must fraction of overlap
    double limit = sum * fraction;
    sum = 0;
    iter = imref.bands.begin();
    for (iter = imref.bands.begin(); iter != imref.bands.end(); iter++) {
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
    std::set<double>::iterator iter;
    for (iter = redshifts.begin(); iter != redshifts.end(); iter++) {
      if (iter == redshifts.begin())
	min_dist = fabs(cosmo.properDist(*iter,z));
      else {
	double dist = fabs(cosmo.properDist(*iter,z));
	if (dist < min_dist)
	  min_dist = dist;
	else // since redshifts are ordered by redshift, we can stop 
	  break;
      }
    }
    // since we have iterated once too much, decrement iterator
    // to get modellist with closest distance
    iter--;
    return *iter;
  }

  void SourceCatalog::computeADUinBands(const Telescope& tel) {
    for (SourceCatalog::iterator iter = SourceCatalog::begin(); iter != SourceCatalog::end(); iter++)
      computeADU(iter->second,tel);
  }

  void SourceCatalog::computeADU(GalaxyInfo& info, const Telescope& tel) {
    sed s = imref.seds.find(info.sed)->second;
    s.shift(info.redshift);
    // need sed normalization: 
    // use average of measured flux (from info.mags) / sed_flux in same band
    double norm = 0, weight = 0, flux, flux_, flux_error;
    for (std::map<std::string,std::pair<double,double> >::const_iterator iter = info.mags.begin(); iter != info.mags.end(); iter++) {
      sed s_ = s;
      for (std::set<Band>::iterator siter = imref.bands.begin(); siter != imref.bands.end(); siter++) {
	// name of imref band is the one of the magnitude measurements
	if (siter->name == iter->first) {
	  s_ *= siter->curve;
	  // flux in filter with filter Qe correction
	  flux_ = s_.getNorm() / siter->curve.getQe(); 
	  flux = Conversion::mag2flux(iter->second.first);        // flux from magnitude
	  flux_error = Conversion::mag2flux(iter->second.second); // weight from magnitude error
	  weight += 1./flux_error;
	  norm += (flux/flux_)/flux_error;
	}
      }
    }
    if (norm != 0 && weight != 0) { // some magnitudes were available
      norm /= weight;
      s *= tel.total;
      flux = norm * s.getNorm() / tel.total.getQe();
      info.mag = Conversion::flux2mag(flux);
      norm = 0; //reused below
      if (info.model_type > 0) {
	// find overlap of the filtered s with imref's bands to get their relative weights
	for (std::set<Band>::iterator siter = imref.bands.begin(); siter != imref.bands.end(); siter++) {
	  // only consider bands with non-vanishing overlap to tel.total
	  if (siter->overlap > 0) {
	    sed s_ = s;
	    s_ *= siter->curve;
	    norm += s_.getNorm();
	    info.adus[siter->name] = s_.getNorm();
	  }
	}
      } else { // only one band
	info.adus["tel"] = norm = 1;
      }

      // for each entry in info.adus: multiply with flux and devide by sum
      // then do conversion to photons and ADUS
      for (std::map<std::string, double>::iterator iter = info.adus.begin(); iter != info.adus.end(); iter++) {
	iter->second *= flux/norm;
	// FIXME: is tel.total the corect transmission curve here???
	iter->second = Conversion::flux2photons(iter->second,1,tel,tel.total);
	iter->second = Conversion::photons2ADU(iter->second,tel.gain);
	// account for size of tel's pixels
	iter->second *= gsl_pow_2(tel.pixsize);
      }
    }
  }

  void SourceCatalog::save(std::string filename) const {
    std::ofstream ofs (filename.c_str());
    ofs << "# SourceCatalog file" << std::endl;
    ofs << "# SQL: " << query << std::endl;
    for (SourceCatalog::const_iterator iter = SourceCatalog::begin(); iter != SourceCatalog::end(); iter++) {
      ofs << iter->second.object_id << "\t" << iter->second.redshift << "\t";
      ofs << iter->second.sed << "\t" << iter->second.radius << "\t";
      ofs << iter->second.ellipticity << "\t" << iter->second.n_sersic << "\t";
      for (std::map<std::string, std::pair<double,double> >::const_iterator miter = iter->second.mags.begin(); miter != iter->second.mags.end(); miter++) {
	ofs << miter->first << ":" << miter->second.first << ":" << miter->second.second;
	if (miter != --(iter->second.mags.end()))
	  ofs << ", ";
      } 
      ofs << "\t" << iter->second.model_type << "\t";
      for (std::map<std::string, double>::const_iterator miter = iter->second.adus.begin(); miter != iter->second.adus.end(); miter++) {
	ofs << miter->first << ":" << miter->second;
	if (miter != --(iter->second.adus.end()))
	  ofs << ", ";
      }
      ofs << "\t" << iter->second.mag << "\t";
      ofs << iter->second.centroid(0) << "\t" << iter->second.centroid(1) << "\t";
      ofs << iter->second.rotation << "\t" << iter->second.redshift_layer << std::endl;
    }
    ofs.close();
  }

  void SourceCatalog::createGalaxyLayers() const {
    // FIXME: implementation
    // model_type = 0: use catalog infos to create SersicModels
    // model_type = 1: query ShapeletDB with a join query 
    //  (select subset of the catalog query with model_type = 1)
    //  use a sorted lookup table to relate entries in catalog with entries in ShapeletObjectList
  }

} // end namespace
