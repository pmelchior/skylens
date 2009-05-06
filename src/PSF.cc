#include "../include/PSF.h"

using namespace skylens;

PSF::PSF() {
}

PSF::PSF(std::string filename) {
  psf = shapelens::ShapeletObject(filename);
}

const shapelens::ShapeletObject& PSF::getShape() const {
  return psf;
}
