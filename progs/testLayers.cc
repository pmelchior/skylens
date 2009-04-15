#include <Layer.h>
#include <iostream>



using namespace skylens;
using namespace shapelens;

int main() {
  LensingLayer ll1(0.5,"");
  LensingLayer ll2(1,"");
  GalaxyLayer lg1(0.75);
  GalaxyLayer lg2(2);
  PSF psf("data/SUBARU/psf.fits");
  StarLayer lS(1e-3,psf);
  
  
  std::list<Polygon<double> > masks;
  std::list<Point2D<double> > points;
  points.push_back(Point2D<double>(0,0));
  points.push_back(Point2D<double>(0,10));
  points.push_back(Point2D<double>(10,10));
  masks.push_back(Polygon<double>(points));
  points.clear();
  points.push_back(Point2D<double>(10,10));
  points.push_back(Point2D<double>(10,12));
  points.push_back(Point2D<double>(12,12));
  points.push_back(Point2D<double>(12,10));
  masks.push_back(Polygon<double>(points));

  MaskLayer ml(-1,masks);
  DitherLayer ld(-2,0.1,0.2);

  LayerStack& ls = SingleLayerStack::getInstance();
  for (LayerStack::iterator iter = ls.begin(); iter != ls.end(); iter++) {
    std::cout << iter->second->getRedshift() << "\t" << iter->second->getType() << std::endl;
  }
  
  std::cout << "total flux at (0,0) = " << ls.begin()->second->getFlux(0,0) << std::endl;

}
