#include <Telescope.h>
#include <Conventions.h>
#include <fstream>
#include <stdexcept>
#include <iostream>

using namespace skylens;

double Telescope::getDiameter() const {
  return d;
}

double Telescope::getPixelScale() const {
  return px;
}

double Telescope::getGain() const {
  return gain;
}

double Telescope::getReadOutNoise() const {
  return ron;
}

std::string Telescope::getName() const {
  return name;
}

std::string Telescope::getFilterName() const {
  return band;
}

const filter& Telescope::getFilter() const {
  return total;
}

const PSF& Telescope::getPSF() const {
  return psf;
}

void Telescope::readConfig(std::string path) {
  std::ifstream ifs((path+"/telescope.conf").c_str());
  if (ifs.bad())
    std::cerr << "Telescope: telescope.conf missing" << std::endl;

  std::string key;
  double value,qe_ccd,qe_mirror,qe_optics;
  qe_ccd = qe_mirror = qe_optics = 0; // OK for comparison
  while(ifs >> key >> value) {
    if (key == "#DIAMETER")
      d = value;
    else if (key == "#FLAT-ACC")
      flat_acc = value;
    else if (key == "#GAIN")
      gain = value;
    else if (key == "#PIXSIZE")
      px = value;
    else if (key == "#RON")
      ron = value;
    else if (key == "#CCD")
      qe_ccd = value;
    else if (key == "#MIRROR")
      qe_mirror = value;
    else if (key == "#OPTICS")
      qe_optics = value;
   }
  // either open spectral shape files
  if (qe_ccd == 0)
    total*=filter("ccd.fits",path);
  // or multiply with constant qe
  else
    total*=qe_ccd;
  // same for mirror ...
  if (qe_mirror == 0)
    total*=filter("mirror.fits",path);
  else
    total*=qe_mirror;
  /// ... and optics
  if (qe_optics == 0)
    total*=filter("optics.fits",path);
  else
    total*=qe_optics;
}

SUBARU::SUBARU(std::string b) : Telescope() {
  name = "SUBARU";
  band = b;
  std::string path = datapath + "/" + name;

  // open all spectral curves files
  try {
    total = filter("filter_"+band+".fits",path);
    readConfig(path);
    psf = PSF(path + "/psf.fits");
  } catch (std::exception & e) {
    throw std::invalid_argument("SUBARU: configuration or data files missing!\n");//+std::string(e.what));
  }
}

HST::HST(std::string b, std::string i) {
  name = "HST";
  band = b;
  std::string path = datapath + "/" + name + "/" + i;

  try {
    total = filter("filter_"+band+".fits",path);
    readConfig(path);
    psf = PSF(path + "/psf.fits");
  } catch (std::exception & e) {
    std::cerr << "HST_" << i <<": configuration or data files missing!" << std::endl;
  }
}
