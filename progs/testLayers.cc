#include <Layer.h>
#include <iostream>

using namespace skylens;
using namespace shapelens;

int main() {
  // LensingLayer ll1(0.5,"");
//   LensingLayer ll2(1,"");
//   GalaxyLayer lg1(0.75);
//   GalaxyLayer lg2(2);
//   PSF psf("data/SUBARU/psf.fits");
//   StarLayer lS(1e-3,psf);
//   DitherLayer ld(-2,0.1,0.2);
  
  std::list<Polygon<double> > masks;
  std::list<Point2D<double> > points;
  points.push_back(Point2D<double>(0,0));
  points.push_back(Point2D<double>(2,10));
  points.push_back(Point2D<double>(10,10));
  masks.push_back(Polygon<double>(points));
  points.clear();
  points.push_back(Point2D<double>(10,10));
  points.push_back(Point2D<double>(10,12));
  points.push_back(Point2D<double>(12,12));
  points.push_back(Point2D<double>(12,10));
  masks.push_back(Polygon<double>(points));
  std::ofstream ofs("mask.txt");
  for (std::list<Polygon<double> >::iterator iter = masks.begin(); iter != masks.end(); iter++)
    ofs << *iter;
  ofs.close();

  MaskLayer ml(-1,1,"mask.txt");
  ConstFluxLayer cl(1,1e3);

  LayerStack& ls = SingleLayerStack::getInstance();
  
  Image<double> im(20,20);
  Layer* front = ls.begin()->second;
  for (int i=0; i < im.getSize(0); i++)
    for (int j=0; j < im.getSize(1); j++)
      im(i,j) = front->getFlux(i+0.5,j+0.5); // centered pixellation
  im.save("testLayer.fits");

  for (LayerStack::iterator iter = ls.begin(); iter != ls.end(); iter++) {
    std::cout << iter->second->getRedshift() << "\t" << iter->second->getType() << std::endl;
  }
}
