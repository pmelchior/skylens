#include <skylens/Layer.h>
#include <skylens/Telescope.h>
#include <skylens/Observation.h>
#include <skylens/Conventions.h>
#include <skylens/Conversion.h>
#include <skylens/RNG.h>
#include <iostream>
#include <map>
#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>
#include <shapelens/ShapeLens.h>
#include <shapelens/utils/MySQLDB.h>
#include <time.h>
#include <tclap/CmdLine.h>

using namespace skylens;
using namespace shapelens;
using boost::format;

struct GalaxyInfo {
  std::map<std::string, data_t> mags;
  data_t redshift;
  data_t sed;
  data_t radius;
  data_t ellipticity;
  data_t n_sersic;
  unsigned int model_type;
  
  friend std::ostream& operator<<(std::ostream& os, const GalaxyInfo& info) {
    for (std::map<std::string, data_t>::const_iterator it = info.mags.begin(); it != info.mags.end(); it++)
      os << it->first << "\t" << it->second << "\t";
    os << info.redshift << "\t" << info.sed << "\t";
    os << info.radius << "\t" << info.ellipticity << "\t";
    os << info.n_sersic << "\t" << info.model_type << std::endl;
    return os;
  }
};
typedef std::map<unsigned long, GalaxyInfo> SkyCat;

struct Band {
  int lambda_c;
  filter curve;
};

struct ImagingReference {
  data_t pixsize;
  data_t fov;
  std::string table;
  std::map<std::string, Band> bands;
  std::map<data_t, sed> seds;
};

data_t computeADU(const Telescope& tel, const filter& transmission, data_t time, const GalaxyInfo& info, const ImagingReference& ref) {
  sed s = ref.seds.find(info.sed)->second;
  s.shift(info.redshift);
  // need sed normalization: 
  // use average of measured flux (from info.mags) / sed_flux in same band
  // ATTENTION: this assumes equal weight for all bands, ignores errors
  data_t norm = 0, flux, flux_;
  int count = 0;
  for (std::map<std::string,data_t>::const_iterator iter = info.mags.begin(); iter != info.mags.end(); iter++) {
    if (iter->second != 0 && iter->second < 30) {
      count++;
      sed s_ = s;
      s_ *= ref.bands.find(iter->first)->second.curve;
      flux_ = s_.getNorm() / ref.bands.find(iter->first)->second.curve.getQe();
      flux = Conversion::mag2flux(iter->second);
      norm += flux/flux_;
    }
  }
  if (norm != 0) {
    norm /= count;
    s *= tel.total;
    flux_ = norm * s.getNorm() / tel.total.getQe();
    double photons_ = Conversion::flux2photons(flux_,time,tel,transmission);
    double ADU_ = Conversion::photons2ADU(photons_,tel.gain);
    ADU_ *= gsl_pow_2(tel.pixsize);
    return ADU_;
  } else
    return 0;
}

data_t computeADU(const Telescope& tel, const filter& transmission, data_t time, const GalaxyInfo& info, const std::map<std::string, data_t>& overlap_bands) {
  data_t mag_ = 0, weights = 0, mag;
  for (std::map<std::string,data_t>::const_iterator iter = overlap_bands.begin(); iter != overlap_bands.end(); iter++) {
    mag = info.mags.find(iter->first)->second;
    if (mag > 0) {
      weights += iter->second;
      mag_ += iter->second * mag;
    }
  }
  if (weights > 0) {
    mag_ /= weights;
    double flux_ = Conversion::mag2flux(mag_); 
    double photons_ = Conversion::flux2photons(flux_,time,tel,transmission);
    double ADU_ = Conversion::photons2ADU(photons_,tel.gain);
    ADU_ *= gsl_pow_2(tel.pixsize);
    return ADU_;
  } else
    return 0;
}

std::map<std::string, data_t> getOverlapBands(const Telescope& tel, const ImagingReference& ref) {
  std::map<std::string,data_t> overlap;
  data_t sum = 0;
  for (std::map<std::string, Band>::const_iterator iter = ref.bands.begin(); iter != ref.bands.end(); iter++) {
    filter f = iter->second.curve;
    f *= tel.total;
    overlap[iter->first] = f.getQe();
    sum += f.getQe();
  }
  // to be included band must have > 10% of overlap flux
  data_t limit = sum/10;
  sum = 0;
  std::map<std::string, data_t>::iterator iter;
  for (iter = overlap.begin(); iter != overlap.end(); iter++) {
    if (iter->second < limit)
      overlap.erase(iter);
    else
      sum += iter->second;
  }
  for (iter = overlap.begin(); iter != overlap.end(); iter++)
    iter->second /= sum;

  return overlap;
}

void adjustGalaxyNumbers(double FoV_ref, double FoV_out, SkyCat& skycat, std::map<std::string, ShapeletObjectList>& shapelet_models,std::map<unsigned long, unsigned long>& model_index) {
  unsigned int N_ref = skycat.size();
  unsigned int N_out = (unsigned int) floor(N_ref * FoV_out / FoV_ref);
  RNG& rng = Singleton<RNG>::getInstance();
  const gsl_rng* r = rng.getRNG();
  SkyCat::iterator iter;
  // too many objects in catalog:
  // remove objects randomly from cat
  if (N_out < N_ref) {
    int selected;
    while (skycat.size() > N_out) {
      selected = (int) floor(skycat.size()*gsl_rng_uniform(r)) - 1;
      iter = skycat.begin();
      if (selected > 0)
	advance(iter,selected);
      skycat.erase(iter);
    }
  }
  // too few objects in catalog:
  // replicate objects from catalog (without change of shapelet_models)
  else if (N_out > N_ref) {
    int selected;
    // get maximum id in cat
    iter = skycat.end();
    iter--;
    unsigned int max = iter->first;
    while (skycat.size() < N_out) {
      // select only from the reference objects
      selected = (int) floor(N_ref*gsl_rng_uniform(r)) - 1;
      max++;
      iter = skycat.begin();
      if (selected > 0)
	advance(iter,selected);
      // duplicate: cat[max] = CatInfo[selected]
      skycat[max] = iter->second;
      // add new cat entry to model lookup:
      // cat entry corresponds to existing entry of shapelet_models
      if (iter->second.model_type == 1)
	model_index[max] = model_index[iter->first];
    }
  }
}

// resize to account for different pixel sizes
void adjustShapeletSize(ShapeletObject& sobj, data_t pix_origin) {
  data_t beta_ = sobj.getBeta()*pix_origin;
  sobj.setBeta(beta_);
  Grid grid = sobj.getGrid();
  Point<data_t> centroid = sobj.getCentroid();
  ScalarTransformation<double> S(pix_origin);
  grid.setWCS(S);
  S.transform(centroid);
  sobj.setGrid(grid);
  sobj.setCentroid(centroid);
}

std::map<std::string, ShapeletObjectList> getShapeletModels(const std::map<std::string, data_t>& overlap, const ImagingReference& ref, std::string catalogdb) {
  std::map<std::string, ShapeletObjectList> models;
  ShapeletObjectDB sdb;
  sdb.useHistory(false);
  sdb.useCovariance(false);
  for (std::map<std::string, data_t>::const_iterator iter = overlap.begin(); iter != overlap.end(); iter++) {
    // select correct table from bandname
    std::string table = str(format(ref.table) %(iter->first));
    sdb.selectTable(table);
    // get shapelet models from table
    models[iter->first] = sdb.load("`" + catalogdb + "`.`skylens`.`model_type` = 1","`" + catalogdb + "`.`skylens` ON (`" + catalogdb + "`.`skylens`.`id` = `" + table + "`.`id`)");
    // adjust size of each model to account for change in pixel scale
    ShapeletObjectList& thislist = models[iter->first];
    for (ShapeletObjectList::iterator siter = thislist.begin(); siter != thislist.end(); siter++)
      adjustShapeletSize(*(*siter), ref.pixsize);
  }
  return models;
}

int sign(data_t x) {
  if (x>=0)
    return 1;
  else
    return -1;
}

void setRandomOrthogonalMatrix(NumMatrix<data_t>& M) {
  if (M.getRows() != 2 || M.getColumns() != 2)
    M.resize(2,2);
  RNG& rng = Singleton<RNG>::getInstance();
  const gsl_rng * r = rng.getRNG();
  
  M(0,0) = -1 + 2*gsl_rng_uniform(r);         // in [-1,1)
  M(1,0) = sqrt(1-(M(0,0)*M(0,0)));           // normal column
  M(1,0)*= sign(-1 + 2*gsl_rng_uniform(r));   // random sign
  data_t s = sign(-1 + 2*gsl_rng_uniform(r)); // either rotation or reflection
  M(0,1) = s*M(1,0);
  M(1,1) = -s*M(0,0);
}

int main(int argc, char* argv[]) {
  // parse commandline
  TCLAP::CmdLine cmd("Run SkyLens++ simulator", ' ', "0.2");
  TCLAP::ValueArg<std::string> telname("t","telescope","Name of telescope",true,"","string", cmd);
  TCLAP::ValueArg<std::string> bandname("b","band","Name of filter band",true,"","string", cmd);
  TCLAP::ValueArg<double> expt("e","exposure_time","Exposure time (in seconds)",true,0,"double", cmd);
  TCLAP::ValueArg<std::string> skyspectrum("s","sky_spectrum","SED of the sky",true,"","string");
  TCLAP::ValueArg<double> sky_mag("S","sky_mag","Exposure time (in seconds)",true,0,"double");
  cmd.xorAdd(skyspectrum,sky_mag);
  TCLAP::ValueArg<double> airmass("m","air_mass","Airmass",false,1,"double",cmd);
  TCLAP::ValueArg<std::string> atm("a","atmosphere_spectrum","Atmospheric absorption spectrum",false,"","string");
  TCLAP::ValueArg<double> absorption("A","air_absorption","Average atmospheric absorption",false,0,"double");
  cmd.xorAdd(atm,absorption);
  TCLAP::SwitchArg noise("n","noise","Add noise to image", cmd, false);
  cmd.parse(argc,argv);
  
  // for measuring computation time
  time_t t0,t1;
  t0 = time(NULL);

  Telescope tel(telname.getValue(),bandname.getValue());
  double exptime = expt.getValue();
  Observation obs(tel,exptime);

  if (skyspectrum.isSet())
    obs.createSkyFluxLayer(sed(skyspectrum.getValue(),datapath));
  else
    obs.createSkyFluxLayer(sky_mag.getValue());
  if (atm.isSet())
    obs.computeTransmittance(filter(atm.getValue(),datapath),airmass.getValue());
  else if (absorption.isSet())
    obs.computeTransmittance(absorption.getValue(),airmass.getValue());
  if (noise.isSet())
    obs.setNoise();
  
  const filter& transmittance = obs.getTotalTransmittance();

  // set up reference information
  ImagingReference hudf_ref;
  hudf_ref.pixsize = 0.03;
  hudf_ref.fov = 39600; // in arcsec^2
  hudf_ref.table = "hudf_acs_%s";
  Band b;
  b.lambda_c = 435;
  b.curve = filter("HST/ACS_WFC/filter_F435W.fits",datapath);
  hudf_ref.bands["B"] = b;
  b.lambda_c = 606;
  b.curve = filter("HST/ACS_WFC/filter_F606W.fits",datapath);
  hudf_ref.bands["V"] = b;
  b.lambda_c = 775;
  b.curve = filter("HST/ACS_WFC/filter_F775W.fits",datapath);
  hudf_ref.bands["i"] = b;
  b.lambda_c = 850;
  b.curve = filter("HST/ACS_WFC/filter_F850LP.fits",datapath);
  hudf_ref.bands["z"] = b;
  // hudf_ref.seds[1] = sed("sed/bpz_1.sed.fits",datapath);
//   hudf_ref.seds[1.333] = sed("sed/bpz_1.333.sed.fits",datapath);
//   hudf_ref.seds[1.667] = sed("sed/bpz_1.667.sed.fits",datapath);
//   hudf_ref.seds[2] = sed("sed/bpz_2.sed.fits",datapath);
//   hudf_ref.seds[2.333] = sed("sed/bpz_2.333.sed.fits",datapath);
//   hudf_ref.seds[2.667] = sed("sed/bpz_2.667.sed.fits",datapath);
//   hudf_ref.seds[3] = sed("sed/bpz_3.sed.fits",datapath);
//   hudf_ref.seds[3.333] = sed("sed/bpz_3.333.sed.fits",datapath);
//   hudf_ref.seds[3.667] = sed("sed/bpz_3.667.sed.fits",datapath);
//   hudf_ref.seds[4] = sed("sed/bpz_4.sed.fits",datapath);
//   hudf_ref.seds[4.333] = sed("sed/bpz_4.333.sed.fits",datapath);
//   hudf_ref.seds[4.667] = sed("sed/bpz_4.667.sed.fits",datapath);
//   hudf_ref.seds[5] = sed("sed/bpz_5.sed.fits",datapath);
//   hudf_ref.seds[5.333] = sed("sed/bpz_5.333.sed.fits",datapath);
//   hudf_ref.seds[5.667] = sed("sed/bpz_5.667.sed.fits",datapath);
//   hudf_ref.seds[6] = sed("sed/bpz_6.sed.fits",datapath);
//   hudf_ref.seds[6.333] = sed("sed/bpz_6.333.sed.fits",datapath);
//   hudf_ref.seds[6.667] = sed("sed/bpz_6.667.sed.fits",datapath);
//   hudf_ref.seds[7] = sed("sed/bpz_7.sed.fits",datapath);
//   hudf_ref.seds[7.333] = sed("sed/bpz_7.333.sed.fits",datapath);
//   hudf_ref.seds[7.667] = sed("sed/bpz_7.667.sed.fits",datapath);
//   hudf_ref.seds[8] = sed("sed/bpz_8.sed.fits",datapath);

  // make wavelength-sorted list of all reference bands
  // this is for fast band selection in case of missing bands
  std::map<int, std::string> sorted;
  for (std::map<std::string,Band>::iterator biter = hudf_ref.bands.begin(); biter != hudf_ref.bands.end(); biter++)
    sorted[biter->second.lambda_c] = biter->first;

  // find reference bands with overlap to tel
  std::map<std::string, data_t> overlap_bands = getOverlapBands(tel,hudf_ref);

  // set up galaxies
  // 0) get information from DB
  MySQLDB db;
  db.connect("SHAPELENSDBCONF");
  std::string catalogdb = "galaxies";
  db.selectDatabase(catalogdb);
  DBResult dbr = db.query("SELECT * FROM `skylens` WHERE `i_mag` IS NOT NULL AND `model_type` = 1"); // only select shapelet models
  MYSQL_ROW row;
  GalaxyInfo info;
  SkyCat skycat;
  SkyCat::iterator iter;
  while (row = dbr.getRow()) {
    unsigned long id = boost::lexical_cast<unsigned long>(row[0]);
    if (row[1] != NULL)
      info.mags["B"] = boost::lexical_cast<data_t>(row[1]);
    else
      info.mags["B"] = 0;
    if (row[2] != NULL)
      info.mags["V"] = boost::lexical_cast<data_t>(row[2]);
    else
      info.mags["V"] = 0;
    if (row[3] != NULL)
      info.mags["i"] = boost::lexical_cast<data_t>(row[3]);
    else
      info.mags["i"] = 0;
    if (row[4] != NULL)
      info.mags["z"] = boost::lexical_cast<data_t>(row[4]);
    else
      info.mags["z"] = 0;
    info.redshift = boost::lexical_cast<data_t>(row[5]);
    info.sed = boost::lexical_cast<data_t>(row[6]);
    if (row[7] != NULL)
      info.radius = boost::lexical_cast<data_t>(row[7]);
    else
      info.radius = 0;
    if (row[8] != NULL)
      info.ellipticity = boost::lexical_cast<data_t>(row[8]);
    else
      info.ellipticity = 0;
    if (row[9] != NULL)
      info.n_sersic = boost::lexical_cast<data_t>(row[9]);
    else
      info.n_sersic = 0;
    info.model_type = boost::lexical_cast<data_t>(row[10]);
    
    skycat[id] = info;
  }
  
  // get shapelet models for those objects
  std::map<std::string, ShapeletObjectList> shapelet_models = getShapeletModels(overlap_bands, hudf_ref, catalogdb);
  // create searchable map: object ID -> vector index
  std::map<unsigned long, unsigned long> model_index;
  ShapeletObjectList& sl = shapelet_models.begin()->second;
  for (unsigned long i=0; i < sl.size(); i++)
    model_index[sl[i]->getObjectID()] = i;

  // account for change of FoV with respect to reference
  adjustGalaxyNumbers(hudf_ref.fov, tel.fov_x*tel.fov_y, skycat, shapelet_models, model_index);

  // set up RNG
  RNG& rng = Singleton<RNG>::getInstance();
  const gsl_rng * r = rng.getRNG();

  // populate galaxy layer with objects from skycat
  SourceModelList galaxies;
  std::complex<data_t> I(0,1);
  data_t flux;
  Point<data_t> centroid, zero(0,0);
  NumMatrix<data_t> O(2,2);
  for (iter = skycat.begin(); iter != skycat.end(); iter++) {
    GalaxyInfo& info = iter->second;
    // compute flux
    //  too slow to be usefull, yet
    //data_t flux = computeADU(tel, t, info, hudf_ref);
    data_t flux = computeADU(tel, transmittance, exptime, info, overlap_bands);
    // set random centroid ..
    centroid(0) = tel.fov_x*gsl_rng_uniform(r);
    centroid(1) = tel.fov_y*gsl_rng_uniform(r);
    // ... and rotation/flip matrix
    setRandomOrthogonalMatrix(O);
    // form affine transformation with O and centroid
    // assuming all SourceModels live at (0,0)
    AffineTransformation<data_t> A(O,zero,centroid);
    if (info.model_type == 1) { // Shapelet models
      unsigned long index = model_index[iter->first];
      std::map<std::string,data_t> flux_ = overlap_bands;
      std::map<std::string, data_t>::iterator fiter;
      // check whether all overlapping bands have valid models
      bool missing = false;
      for (std::map<std::string, ShapeletObjectList>::iterator siter = shapelet_models.begin(); siter != shapelet_models.end(); siter++) {
	if (siter->second[index]->getFlags().test(15)) {
	  missing = true;
	  flux_[siter->first] *= -1; // make it negative, but remember it
	}
      }

      // some bands are missing: 
      // shift their flux to next redder available band
      if (missing) {
	data_t flux__ = 0;
	for (std::map<int, std::string>::iterator siter = sorted.begin(); siter != sorted.end(); siter++) {
	  fiter = flux_.find(siter->second);
	  if (fiter != flux_.end()) {
	    if (fiter->second < 0) {
	      flux__ -= fiter->second;
	      flux_.erase(fiter);
	    } else {
	      fiter->second += flux__;
	      flux__ = 0;
	    }
	  }
	}
	// last band had no model: need to add flux to second-to-last
	// available band
	if (flux__ != 0 && flux_.size() > 0) {
	  fiter = flux_.end();
	  fiter--;
	  fiter->second += flux__;
	}
      }

      // add model for all available bands to galaxy list
      for (fiter = flux_.begin(); fiter!= flux_.end(); fiter++) {
	ShapeletObjectList& sl = shapelet_models[fiter->first];
	// weight model with its relative flux
	galaxies.push_back(boost::shared_ptr<SourceModel>(new ShapeletModel(*sl[index],flux * (fiter->second),&A)));
      }
    }
    // else if (info.model_type == 0) { // Sersic model
//       data_t n = iter->second.n_sersic;
//       data_t Re = iter->second.radius;
//       data_t e = iter->second.ellipticity;
//       Re *= hudf_ref.pixsize;   // correct for size change
//       data_t b_a = 1 - e;
//       data_t epsilon = e/(1+b_a);
//       std::complex<data_t> eps(epsilon,0);// random orientation via A
//       galaxies.push_back(boost::shared_ptr<SourceModel>(new SersicModel(n,Re,flux,eps,&A,iter->first)));
    //}
  }

  // define names for outputs
  std::ostringstream filename;
  std::string tname = tel.name;
  std::string::size_type loc = tname.find( "/",0);
  while( loc != std::string::npos ) {
    tname.replace(loc,1,"_");
    loc = tname.find( "/",loc);
  }
  filename << tname << "_" << exptime << "s_" << tel.bandname << ".cat";

  // prepare catalog of source
  // - remove duplicates (from multible band models)
  // - compute flux to mags
  // - adjust positions from arcsec to pixels
  Catalog cat = galaxies.getCatalog();
  Catalog::iterator citer_;
  data_t zeropoint = Conversion::zeroPoint(tel,transmittance,exptime);
  data_t flux_ = 0;
  for (Catalog::iterator citer = cat.begin(); citer != cat.end(); citer++) {
    CatObject& ci = citer->second;
    flux_ = ci.FLUX;
    citer_ = citer;
    citer_++;
    if (citer_!=cat.end()) {
      if (citer_->second.PARENT == ci.PARENT) {
	flux_ += citer_->second.FLUX;
	cat.erase(citer_);
    }
      }
    ci.FLUX = Conversion::flux2mag(flux_/gsl_pow_2(tel.pixsize)*pow(10.,-0.4*(zeropoint+48.6)));
  }
  ScalarTransformation<data_t> S(tel.pixsize);
  cat.apply(*S.getInverse());
  cat.save(filename.str());

  new GalaxyLayer(0.75,galaxies);
  //new ConvolutionLayer(L*tel.pixsize,tel.pixsize,tel.psf);


  Image<double> im;
  obs.makeImage(im);

  filename.str("");
  filename << tname << "_" << exptime << "s_" << tel.bandname << ".fits";
  fitsfile* fptr = IO::createFITSFile(filename.str());
  IO::writeFITSImage(fptr,im);
  IO::updateFITSKeywordString(fptr,"TELESCOPE",tel.name);
  IO::updateFITSKeywordString(fptr,"BANDNAME",tel.bandname);
  IO::updateFITSKeyword(fptr,"EXPTIME",exptime,"exposure time");
  IO::updateFITSKeyword(fptr,"MAGZPT",Conversion::zeroPoint(tel,transmittance,1),"magnitude zeropoint");
  IO::updateFITSKeyword(fptr,"GAIN",tel.gain,"CCD gain");
  IO::updateFITSKeyword(fptr,"RNG_SEED",gsl_rng_default_seed);
  // add WCS parameters here...
  IO::closeFITSFile(fptr);

  t1 = time(NULL);
  std::cout << "Computation time: " << t1-t0 << " seconds" << std::endl;
}
