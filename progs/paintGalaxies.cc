#include <shapelens/ShapeLens.h>
#include <skylens/SkyLens.h>
#include <iostream>
#include <map>
#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>
#include <libastro/cosmology.h>
#include <time.h>
#include <tclap/CmdLine.h>

using namespace skylens;
using namespace shapelens;
using boost::format;

typedef std::map<unsigned long, boost::shared_ptr<ShapeletObject> > ShapeletObjectCat;

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

std::map<std::string, ShapeletObjectCat> getShapeletModels(const std::map<std::string, data_t>& overlap, const ImagingReference& ref, const SkyLensCatalog& skycat) {
  std::map<std::string, ShapeletObjectCat> models;
  ShapeletObjectDB sdb;
  sdb.useHistory(false);
  sdb.useCovariance(false);
  for (std::map<std::string, data_t>::const_iterator iter = overlap.begin(); iter != overlap.end(); iter++) {
    // create entry in models
    models[iter->first] = ShapeletObjectCat();
    // select correct table from bandname
    std::string table = str(format(ref.table) %(iter->first));
    sdb.selectTable(table);
    // get shapelet models from table
    std::string query;
    for (std::map<std::string, std::string>::const_iterator witer = skycat.where.begin(); witer!= skycat.where.end(); witer++) {
      query += "(`" + skycat.dbname + "`.`" + skycat.tablename + "`.`" + witer->first + "` " + witer->second + ")";
      if (witer != --(skycat.where.end()))
	query += " AND ";
    }
    // only load models of object for which we trust model
    query += " AND (`" + skycat.dbname + "`.`" + skycat.tablename + "`.`model_type` = 1)";
    ShapeletObjectList sl = sdb.load(query,"`" + skycat.dbname + "`.`" + skycat.tablename + "` ON (`" + skycat.dbname + "`.`" + skycat.tablename +"`.`id` = `" + table + "`.`id`)");
    // insert each entry of sl into models and
    for (ShapeletObjectList::iterator siter = sl.begin(); siter != sl.end(); siter++)
      models[iter->first][(*siter)->getObjectID()] = *siter;
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

std::map<double, SourceModelList> createGalaxyLists(std::string s) {
  std::map<double, SourceModelList> gls;
  double z;
  std::string::size_type front = s.find(",",0);
  if (front != std::string::npos) {
    std::string::size_type back = 0;
    while (front != std::string::npos) {
      z = boost::lexical_cast<double>(s.substr(back,front-back));
      gls[z] = SourceModelList();
      back = front+1;
      front = s.find(",",back);
    }
    front = s.size();
    z = boost::lexical_cast<double>(s.substr(back,front-back));
    gls[z] = SourceModelList();
  } else {
    z = boost::lexical_cast<double>(s);
    gls[z] = SourceModelList();
  }
  return gls;
}

SourceModelList& findNearestLayer(std::map<double, SourceModelList>& gals, cosmology& cosmo, double redshift) {
  double min_dist;
  std::map<double, SourceModelList>::iterator iter;
  for (iter = gals.begin(); iter != gals.end(); iter++) {
    if (iter == gals.begin())
      min_dist = fabs(cosmo.properDist(iter->first,redshift));
    else {
      double dist = fabs(cosmo.properDist(iter->first,redshift));
      if (dist < min_dist)
	min_dist = dist;
      else // since gals are ordered by redshift, we can stop 
	break;
    }
  }
  // since we have iterated once too much, decrement iterator
  // to get modellist with closest distance
  iter--;
  return iter->second;
}

int main(int argc, char* argv[]) {
  // parse commandline
  TCLAP::CmdLine cmd("Run SkyLens++ simulator", ' ', "0.3");
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
  TCLAP::ValueArg<std::string> gallayerlist("g","galaxy_layers","CSV list of redshifts for the galaxy layers",true,"1","string",cmd);
  cmd.parse(argc,argv);
  
  // for measuring computation time
  time_t t0,t1;
  t0 = time(NULL);

  Telescope tel(telname.getValue(),bandname.getValue());
  double exptime = expt.getValue();
  Observation obs(tel,exptime);
  // vanilla LCDM, change it here to have effect on 
  // all cosmological calculations.
  cosmology& cosmo = SingleCosmology::getInstance();
  
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

  // get reference catalog from skylens DB
  // select according to following properties
  std::map<std::string, std::string> where;
  where["i_mag"] = "IS NOT NULL";
  where["model_type"] = "=1";
  SkyLensCatalog skycat(where);
  // account for change of FoV with respect to reference
  skycat.adjustGalaxyNumber(hudf_ref.fov, tel.fov_x*tel.fov_y);
  
  // get shapelet models for those objects
  std::map<std::string, ShapeletObjectCat> shapelet_models = getShapeletModels(overlap_bands, hudf_ref, skycat);

  // access global RNG
  RNG& rng = Singleton<RNG>::getInstance();
  const gsl_rng * r = rng.getRNG();

  // populate galaxy layer with objects from skycat
  std::map<double, SourceModelList> gals = createGalaxyLists(gallayerlist.getValue());

  std::complex<data_t> I(0,1);
  Point<data_t> zero(0,0);
  NumMatrix<data_t> O(2,2);
  for (SkyLensCatalog::iterator iter = skycat.begin(); iter != skycat.end(); iter++) {
    GalaxyInfo& info = iter->second;
    // determine galaxy layer on which to place this source
    SourceModelList& galaxies = findNearestLayer(gals,cosmo,info.redshift);
    // compute flux
    //  too slow to be usefull, yet
    //data_t flux = computeADU(tel, t, info, hudf_ref);
    info.flux = computeADU(tel, transmittance, exptime, info, overlap_bands);
    // set random centroid ..
    info.centroid(0) = tel.fov_x*gsl_rng_uniform(r);
    info.centroid(1) = tel.fov_y*gsl_rng_uniform(r);
    // ... and rotation/flip matrix
    setRandomOrthogonalMatrix(O);
    // account for original pixe size
    O *= hudf_ref.pixsize;
    // conservation of surface brightness:
    info.flux /= hudf_ref.pixsize*hudf_ref.pixsize; 
    // form affine transformation with O and centroid
    // assuming all SourceModels live at (0,0)
    AffineTransformation<data_t> A(O,zero,info.centroid);
    // store bounding box of model
    if (info.model_type == 1) { // Shapelet models
      std::map<std::string,data_t> flux_ = overlap_bands;
      std::map<std::string, data_t>::iterator fiter;
      // check whether all overlapping bands have valid models
      bool missing = false;
      for (std::map<std::string, ShapeletObjectCat>::iterator siter = shapelet_models.begin(); siter != shapelet_models.end(); siter++) {
	if (siter->second[info.object_id]->getFlags().test(15)) {
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
	ShapeletObjectCat& sl = shapelet_models[fiter->first];
	// weight model with its relative flux
	galaxies.push_back(boost::shared_ptr<SourceModel>(new ShapeletModel(*sl[info.object_id],info.flux * (fiter->second),&A)));
	// get bounding box of this sourcemodel for skycat
	// as they are the same for each of them: once is enough
	if (fiter == flux_.begin())
	  info.bb = galaxies.back()->getSupport();
      }
    }
    else if (info.model_type == 0) { // Sersic model
      data_t e = info.ellipticity;
      data_t b_a = 1 - e;
      data_t epsilon = e/(1+b_a);
      std::complex<data_t> eps(epsilon,0);// random orientation via A
      galaxies.push_back(boost::shared_ptr<SourceModel>(new SersicModel(info.n_sersic, info.radius, info.flux,eps,&A,info.object_id)));
      // get bounding box of this sourcemodel for skycat
      info.bb = galaxies.back()->getSupport();
    }
  }
  // for each SourceModelList: create a GalaxyLayer at appropriate redshift
  for (std::map<double, SourceModelList>::iterator sliter = gals.begin(); sliter != gals.end(); sliter++)
    new GalaxyLayer(sliter->first,sliter->second);

  
  new LensingLayer(0.2975,"data/deflector/alpha_vectors.fits");
  //new ConvolutionLayer(L*tel.pixsize,tel.pixsize,tel.psf);

  // define names for outputs
  std::ostringstream filename;
  std::string tname = tel.name;
  std::string::size_type loc = tname.find( "/",0);
  while( loc != std::string::npos ) {
    tname.replace(loc,1,"_");
    loc = tname.find( "/",loc);
  }
  
  //Image<double> im(1000,1000); 
  //im.grid.setWCS(ScalarTransformation<double>(tel.pixsize)); 
  //obs.makeImage(im,false); 
  
  Image<double> im;
  obs.makeImage(im);

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

  // get catalog of all sources used
  filename.str("");
  filename << tname << "_" << exptime << "s_" << tel.bandname << ".cat";
  Catalog cat = skycat.getCatalog(im.grid.getWCS().getWC2PC());
  cat.save(filename.str());


  t1 = time(NULL);
  std::cout << "Computation time: " << t1-t0 << " seconds" << std::endl;
}
