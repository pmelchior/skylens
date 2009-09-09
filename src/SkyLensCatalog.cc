#include "../include/SkyLensCatalog.h"
#include "../include/RNG.h"
#include <shapelens/utils/MySQLDB.h>
#include <boost/lexical_cast.hpp>
#include <math.h>

namespace skylens {

  SkyLensCatalog::SkyLensCatalog(const std::map<std::string, std::string>& where) : 
    std::map<unsigned long, GalaxyInfo> (), dbname("galaxies"), tablename("skylens"), where(where) {

    shapelens::MySQLDB db;
    db.connect("SHAPELENSDBCONF");
    db.selectDatabase(dbname);
    std::string query = "SELECT * FROM `" + tablename + "` WHERE ";
    for (std::map<std::string, std::string>::const_iterator iter = where.begin(); iter!= where.end(); iter++) {
      query += "(`" + iter->first + "` " + iter->second + ")";
      if (iter != --(where.end()))
	query += " AND ";
    }
    shapelens::DBResult dbr = db.query(query); // only select shapelet models
    MYSQL_ROW row;
    GalaxyInfo info;
    while (row = dbr.getRow()) {
      unsigned long id = boost::lexical_cast<unsigned long>(row[0]);
      if (row[1] != NULL)
	info.mags["B"] = boost::lexical_cast<double>(row[1]);
      else
	info.mags["B"] = 0;
      if (row[2] != NULL)
	info.mags["V"] = boost::lexical_cast<double>(row[2]);
      else
	info.mags["V"] = 0;
      if (row[3] != NULL)
	info.mags["i"] = boost::lexical_cast<double>(row[3]);
      else
	info.mags["i"] = 0;
      if (row[4] != NULL)
	info.mags["z"] = boost::lexical_cast<double>(row[4]);
      else
	info.mags["z"] = 0;
      info.redshift = boost::lexical_cast<double>(row[5]);
      info.sed = boost::lexical_cast<double>(row[6]);
      if (row[7] != NULL)
	info.radius = boost::lexical_cast<double>(row[7]);
      else
	info.radius = 0;
      if (row[8] != NULL)
	info.ellipticity = boost::lexical_cast<double>(row[8]);
      else
	info.ellipticity = 0;
      if (row[9] != NULL)
	info.n_sersic = boost::lexical_cast<double>(row[9]);
      else
	info.n_sersic = 0;
      info.model_type = boost::lexical_cast<double>(row[10]);
      info.object_id = id;
      std::map<unsigned long, GalaxyInfo>::insert(std::pair<unsigned long, GalaxyInfo>(id,info));
    }
  }

  void SkyLensCatalog::adjustGalaxyNumber(double FoV_ref, double FoV_out) {
    unsigned int N_ref = std::map<unsigned long, GalaxyInfo>::size();
    unsigned int N_out = (unsigned int) floor(N_ref * FoV_out / FoV_ref);
    RNG& rng = shapelens::Singleton<RNG>::getInstance();
    const gsl_rng* r = rng.getRNG();
    SkyLensCatalog::iterator iter;
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

  shapelens::Catalog SkyLensCatalog::getCatalog(const shapelens::CoordinateTransformation<double>& wc2pc) {
    shapelens::Catalog cat;
    shapelens::CatObject co;
    shapelens::Point<double> P;
    shapelens::Rectangle<double> bb;
    co.FLAGS = 0;
    for (std::map<unsigned long, GalaxyInfo>::iterator iter = std::map<unsigned long, GalaxyInfo>::begin(); iter != std::map<unsigned long, GalaxyInfo>::end(); iter++) {
      bb = iter->second.bb;
      bb.apply(wc2pc);
      co.XMIN = bb.ll(0);
      co.YMIN = bb.ll(1);
      co.XMAX = bb.tr(0);
      co.YMAX = bb.tr(1);
      P = iter->second.centroid;
      wc2pc.transform(P);
      co.XCENTROID = P(0);
      co.YCENTROID = P(1);
      co.FLUX = iter->second.flux;
      co.CLASSIFIER = iter->second.model_type;
      co.PARENT = iter->second.object_id;
      cat[iter->first] = co;
    }
    return cat;
  }

} // end namespace
