#include <Observation.h>
#include <Layer.h>
#include <utils/IO.h>
#include <Conventions.h>
using namespace skylens;
using namespace shapelens;

int main() {
  SUBARU subaru("I");
  filter sky("sky/zodiacal.fits",datapath);
  Observation obs(subaru,sky,1000);
  LayerStack& ls = SingleLayerStack::getInstance();
  for (LayerStack::iterator iter = ls.begin(); iter != ls.end(); iter++) {
    std::cout << iter->second->getRedshift() << "\t" << iter->second->getType();
    if (iter->second->getType() == "SS")
      std::cout << "\tsky brightness = " << iter->second->getFlux(0,0);
    std::cout << std::endl;
  }
  
  Image<double> im(100,100);
  obs.makeImage(im,false);
  fitsfile* fptr = IO::createFITSFile("testObservation.fits");
  IO::writeFITSImage(fptr,im);
  IO::closeFITSFile(fptr);
  return 0;
}
