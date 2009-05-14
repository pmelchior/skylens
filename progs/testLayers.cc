#include <skylens/Layer.h>
#include <gsl/gsl_randist.h>
#include <iostream>
#include <skydb/ObjectList.h>
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
  std::list<Point2D<double> > points;
  points.push_back(Point2D<double>(0.,0.));
  points.push_back(Point2D<double>(0.,double(L)/4));
  points.push_back(Point2D<double>(double(L)/4,double(L)/4));
  points.push_back(Point2D<double>(double(L)/4,0.));
  masks.push_back(Polygon<double>(points));
  //MaskLayer ml(masks);
  //SkyFluxLayer sky(2e2);

  const gsl_rng_type * T;
  gsl_rng * r;
  T = gsl_rng_mt19937;
  r = gsl_rng_alloc (T);
  SourceModelList galaxies;
  // *** SERSIC GALAXIES ***
  // for (int i=0; i < N; i++) {
//     double n = MIN_N + (MAX_N - MIN_N) * gsl_rng_uniform(r);
//     double Re = MIN_RE + (MAX_RE - MIN_RE) * gsl_rng_uniform(r);
//     Point2D<double> centroid(L*gsl_rng_uniform(r),L*gsl_rng_uniform(r));
//     complex<double> eps(gsl_ran_gaussian (r,EPS_STD),gsl_ran_gaussian (r,EPS_STD));
//     galaxies.push_back(boost::shared_ptr<SourceModel>(new SersicModel(n,Re,1e4,eps,centroid))); 
//   }
  
  // *** SHAPELET GALAXIES ***
  std::map<std::string,std::string> where;
  std::ostringstream num;
  num << N;
  where["hudf_acs"] = "1 ORDER BY RAND() LIMIT " + num.str();
  skydb::ObjectList ol(where);
  for (int i=0; i < ol.size(); i++) {
    const skydb::ImagingFiles* imaging = ol[i]->getImagingFiles();
    if (imaging != NULL) {
      std::map<std::string, std::string>::const_iterator iter = imaging->model.find("b");      
      if (iter != imaging->model.end()) {
	std::string siffile_b = iter->second;
	iter = imaging->model.find("z");
	std::string siffile_z = iter->second;
	ShapeletObject* b = new ShapeletObject(siffile_b);
	ShapeletObject* z = new ShapeletObject(siffile_z);
	Point2D<double> centroid(L*gsl_rng_uniform(r),L*gsl_rng_uniform(r));
	galaxies.push_back(boost::shared_ptr<SourceModel>(new ShapeletModel(*b,b->getShapeletFlux(),centroid)));
	galaxies.push_back(boost::shared_ptr<SourceModel>(new ShapeletModel(*z,z->getShapeletFlux(),centroid)));
      }
    }
  } 
  GalaxyLayer lg1(0.75,galaxies);

  // SourceModelList stars;
//   for (int i=0; i < NSTARS; i++) {
//     Point2D<double> centroid(L*gsl_rng_uniform(r),L*gsl_rng_uniform(r));
//     stars.push_back(boost::shared_ptr<SourceModel>(new ShapeletModel(psf.getShape(),5e4,centroid)));
//   }
//   StarLayer ls1(stars);

//   ConvolutionLayer lc(L,1,psf,0);


  LayerStack& ls = SingleLayerStack::getInstance();
  
  Image<double> im(L,L);
  fitsfile* fptr = IO::createFITSFile("testLayer.fits");
  Layer* front = ls.begin()->second;
  for (int i=0; i < im.getSize(0); i++)
    for (int j=0; j < im.getSize(1); j++)
      im(i,j) = front->getFlux(i+0.5,j+0.5); // centered pixellation
  IO::writeFITSImage(fptr,im);
  IO::closeFITSFile(fptr);

  for (LayerStack::iterator iter = ls.begin(); iter != ls.end(); iter++) {
    std::cout << iter->second->getRedshift() << "\t" << iter->second->getType() << std::endl;
 }
  return 0;
}
