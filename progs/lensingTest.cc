#include <shapelens/ShapeLens.h>
#include <skylens/SkyLens.h>
#include <boost/lexical_cast.hpp>

using namespace skylens;
using namespace shapelens;

complex<data_t> ellipticity(const Quadrupole& Q) {
  complex<data_t> I(0,1);
  complex<data_t> Q11(Q(0,0),0),Q22(Q(1,1),0),Q12(Q(0,1),0);
  return (Q11 - Q22 + data_t(2)*I*Q12)/(Q11+Q22 + 2.*sqrt(Q11*Q22-Q12*Q12));
}

int main() {
  // some definitions
  Telescope tel;
  tel.pixsize = 0.14871;
  tel.fov_x = tel.fov_y = 300;
  double exptime = 1000;
  Observation obs(tel,exptime);
  int n = 100;
  int N = n*n;

  // set up galaxies
  SourceModelList gals;
  Point<data_t> centroid;
  double epsilon = 0, flux = 1*tel.pixsize, n_sersic = 0.5, radius = 0.3;
  for (int i=0; i < N; i++) {
    centroid(0) = (0.5+(i%n))/n * tel.fov_x;
    centroid(1) = (0.5+(i/n))/n * tel.fov_y;
    ShiftTransformation<data_t> T(centroid);
    std::complex<data_t> eps(epsilon,0);
    gals.push_back(boost::shared_ptr<SourceModel>(new SersicModel(n_sersic, radius, flux,eps,&T,i+1)));
  }
  // place them on layer
  new GalaxyLayer(0.6,gals);

  // define the lens
  new LensingLayer(0.2975,"data/deflector/alpha_vectors.fits");

  // make an image: use Frame as we want to find objects later
  Frame f;
  obs.makeImage(f);

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
      obj.computeMoments();
      // use obj.Q, obj.O, obj.H from now on
      complex<data_t> epsilon = ellipticity(obj.Q);
      ofs << id << " " << obj.centroid(0) << " " << obj.centroid(1) << " " << real(epsilon) << " " << imag(epsilon) << std::endl;
    }
  }
  ofs.close();
  
}