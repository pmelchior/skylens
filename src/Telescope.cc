#include "../include/Telescope.h"
#include "../include/Conventions.h"
#include "../include/Layer.h"
#include <fstream>
#include <stdexcept>
#include <iostream>

using namespace skylens;

void Telescope::readConfig(std::string path) {
  std::ifstream ifs((path+"/telescope.conf").c_str());
  if (ifs.bad())
    std::cerr << "Telescope: telescope.conf missing" << std::endl;

  std::string key;
  double value,qe_ccd,qe_mirror,qe_optics;
  qe_ccd = qe_mirror = qe_optics = 0; // OK for comparison
  while(ifs >> key >> value) {
    if (key == "#DIAMETER")
      diameter = 100*value; // convert to centimeter
    else if (key == "#FLAT-ACC")
      flat_acc = value;
    else if (key == "#GAIN")
      gain = value;
    else if (key == "#PIXSIZE")
      pixsize = value;
    else if (key == "#FOV_X")
      fov_x = value;
    else if (key == "#FOV_Y")
      fov_y = value;
    else if (key == "#RON")
      ron = value;
    else if (key == "#CCD")
      qe_ccd = value;
    else if (key == "#MIRROR")
      qe_mirror = value;
    else if (key == "#OPTICS")
      qe_optics = value;
  }
  ifs.close();

  total = band;
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

  // check for mask.txt
  ifs.open((path+"/mask.txt").c_str());
  // if present: create MaskLayer from it
  if (ifs.good()) {
    new MaskLayer(path+"/mask.txt");
  }
  ifs.close();
}

SUBARU::SUBARU(std::string b) : Telescope() {
  name = "SUBARU";
  bandname = b;
  std::string path = datapath + "/" + name;

  // open all spectral curves files
  try {
    band = filter("filter_"+bandname+".fits",path);
    readConfig(path);
    psf = PSF(path + "/psf.sif");
  } catch (std::exception & e) {
    throw std::invalid_argument("SUBARU: configuration or data files missing!\n");//+std::string(e.what));
  }
}

HST::HST(std::string b, std::string i) {
  name = "HST";
  bandname = b;
  std::string path = datapath + "/" + name + "/" + i;

  try {
    band = filter("filter_"+bandname+".fits",path);
    readConfig(path);
    psf = PSF(path + "/psf.sif");
  } catch (std::exception & e) {
    std::cerr << "HST_" << i <<": configuration or data files missing!" << std::endl;
  }
}
