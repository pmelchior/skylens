#include <skylens/Layer.h>
#include <gsl/gsl_randist.h>
#include <iostream>
#include <shapelens/ShapeLens.h>

using namespace skylens;
using namespace shapelens;

#define N 1000
#define NSTARS 10
#define L 1024
#define PAD 12
#define MIN_N 0.25
#define MAX_N 4
#define MIN_RE 2
#define MAX_RE 10
#define EPS_STD 0.2

int main() {
  PSF psf("data/SUBARU/psf.sif");
  DitherLayer ld(0.2,0.5);
 
  std::list<Polygon<double> > masks;
  std::list<Point<double> > points;
  points.push_back(Point<double>(0.,0.));
  points.push_back(Point<double>(0.,double(L)/4));
  points.push_back(Point<double>(double(L)/4,double(L)/4));
  points.push_back(Point<double>(double(L)/4,0.));
  masks.push_back(Polygon<double>(points));
  //MaskLayer ml(masks);
  //SkyFluxLayer sky(2e2);

  const gsl_rng_type * T;
  gsl_rng * r;
  T = gsl_rng_mt19937;
  r = gsl_rng_alloc (T);
  SourceModelList galaxies;
  // *** SERSIC GALAXIES ***
  for (int i=0; i < N; i++) {
    double n = MIN_N + (MAX_N - MIN_N) * gsl_rng_uniform(r);
    double Re = MIN_RE + (MAX_RE - MIN_RE) * gsl_rng_uniform(r);
    Point<double> centroid(L*gsl_rng_uniform(r),L*gsl_rng_uniform(r));
    complex<double> eps(gsl_ran_gaussian (r,EPS_STD),gsl_ran_gaussian (r,EPS_STD));
    galaxies.push_back(boost::shared_ptr<SourceModel>(new SersicModel(n,Re,1e4,eps,centroid))); 
  }
  GalaxyLayer lg1(0.75,galaxies);

  // SourceModelList stars;
//   for (int i=0; i < NSTARS; i++) {
//     Point<double> centroid(L*gsl_rng_uniform(r),L*gsl_rng_uniform(r));
//     stars.push_back(boost::shared_ptr<SourceModel>(new ShapeletModel(psf.getShape(),5e4,centroid)));
//   }
//   StarLayer ls1(stars);

//   ConvolutionLayer lc(L,1,psf,0);


  LayerStack& ls = SingleLayerStack::getInstance();
  
  Image<double> im(L,L);
  fitsfile* fptr = IO::createFITSFile("testLayer.fits");
  Layer* front = ls.begin()->second;
  for (unsigned long i=0; i < im.size(); i++)
    im(i) = front->getFlux(im.grid(i));
  IO::writeFITSImage(fptr,im);
  IO::closeFITSFile(fptr);

  std::cout << ls;

  return 0;
}
