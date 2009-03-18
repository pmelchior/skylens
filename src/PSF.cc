#include <PSF.h>
#include <shapelets/ShapeletObject.h>

using namespace skylens;

PSF::PSF() {
}

PSF::PSF(std::string filename) {
  // check for file extension
  // 1) .fits -> Image
  // 2) .sif  -> ShapeletObject
  size_t pos = filename.find_last_of(".");
  std::string extension = filename.substr(pos);
  if (extension == ".fits")
    psf = shapelens::Image<double>(filename);
  else if (extension == ".sif") {
    shapelens::ShapeletObject sobj(filename);
    psf = sobj.getModel();
  }
}

const shapelens::Object& PSF::getShape() const {
  return psf;
}
