#include <shapelens/ShapeLens.h>
#include <skylens/SkyLens.h>

using namespace skylens;
using namespace shapelens;

complex<data_t> epsilon(const Moment2& m2) {
  complex<data_t> e(m2(0,0) - m2(1,1),2*m2(0,1));
  e/= (complex<data_t>(m2(0,0) + m2(1,1)) + 2.*sqrt(complex<data_t>(m2(0,0)*m2(1,1) - gsl_pow_2(m2(0,1)))));
  return e;
}

int main() {
  // some definitions
  Telescope tel;
  tel.pixsize = 0.2;
  tel.fov_x = tel.fov_y = 1000;
  double exptime = 1000;
  Observation obs(tel,exptime);
  double n = 40; // number density of gals per arcmin^2
  double N = tel.fov_x*tel.fov_y / 3600 * n; // total number of gals in FoV
  int L = (int) floor(sqrt(N)); // avg. distance beween N gals in FoV

  // access global RNG
  RNG& rng = Singleton<RNG>::getInstance();
  const gsl_rng * r = rng.getRNG();

  // set up galaxies
  SourceModelList gals;
  Point<data_t> centroid;
  double eps_intr = 0., flux = 1*tel.pixsize, n_sersic = 1.5, radius = 0.35;
  complex<data_t> I(0,1);
  for (int i=0; i < N; i++) {
    centroid(0) = (0.5+(i%L))/L * tel.fov_x;
    centroid(1) = (0.5+(i/L))/L * tel.fov_y;
    ShiftTransformation<data_t> T(centroid);
    std::complex<data_t> eps(eps_intr,0);
    eps *= exp(I*2.*M_PI*gsl_rng_uniform(r));
    gals.push_back(boost::shared_ptr<SourceModel>(new SersicModel(n_sersic, radius, flux,eps,&T,i+1)));
  }
  // place them on layer
  new GalaxyLayer(0.6,gals);

  // define the lens
  LensingLayer* ll = new LensingLayer(0.2975,"data/deflector/g_cluster.fits");

  // make an image: use Frame as we want to find objects later
  // use 1000 pixels centered at the cluster center
  Frame f;
  int F = 1000;
  f.resize(F*F);
  f.grid.setSize(0,0,F,F);
  Point<data_t> cluster_center = ll->getCenter(), frame_center(F/2,F/2);
  NumMatrix<data_t> S(2,2);
  S(0,0) = S(1,1) = tel.pixsize;
  f.grid.setWCS(AffineTransformation<data_t>(S,frame_center,cluster_center));
  obs.makeImage(f,false);

  // and store it
  fitsfile* fptr = IO::createFITSFile("lensingTest.fits");
  IO::writeFITSImage(fptr,f);
  IO::closeFITSFile(fptr);

  // detect objects
  ShapeLensConfig::ADD_BORDER = 0.1;
  ShapeLensConfig::DETECT_THRESHOLD = 0.1;
  ShapeLensConfig::MIN_THRESHOLD = 1e-6;
  f.setNoiseMeanRMS(0,0);
  f.findObjects();
  
  // open file for storing lensing parameters
  std::ofstream ofs("lensingTest.lenscat");
    
  // run through all objects
  const Catalog& cat = f.getCatalog();
  Catalog::const_iterator iter;
  for(iter = cat.begin(); iter != cat.end(); iter++) {
    // for clearity:
    unsigned long id = (*iter).first;
    // select all correctly segmented objects
    if ((*iter).second.FLAGS == 0) {
      // "cut out" the object from whole frame and put it into Object obj
      Object obj;
      f.fillObject(obj,iter);
      Moment2 m2(obj);
      complex<data_t> eps = epsilon(m2);
      ofs << id << " " << obj.centroid(0) << " " << obj.centroid(1) << " " << real(eps) << " " << imag(eps) << std::endl;
    }
  }
  ofs.close();
  
}
