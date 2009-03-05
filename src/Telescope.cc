#include <Telescope.h>
#include <Conventions.h>
#include <fstream>
#include <CCfits/CCfits>

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
  return filter_name;
}

const filter& Telescope::getFilter() const {
  return total;
}

void Telescope::readConfig(std::string path) {
  std::ifstream ifs((path+"/telescope.conf").c_str());
  if (ifs.bad())
    std::cerr << "Telescope: telescope.conf missing" << std::endl;

  std::string key;
  double value;
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
  }
}

SUBARU::SUBARU(std::string band) {
  std::string path = datapath + "/SUBARU";
  Telescope::readConfig(path);
  Telescope::filter_name = band;
  // open all spectral curves files
  try {
    filter bandf("filter_"+band+".fits",path);
    filter ccd("ccd.fits",path);
    // mirror and optics are just set to 90% each, so
    // we multiply ccd with 0.81
    ccd*=0.81;
    // total = product spectral curve
    total = bandf*ccd;
  } catch (CCfits::FitsException& e) {
    std::cerr << "SUBARU: configuration files missing!" << std::endl;
    std::cerr << "        check for filter_"+band+".fits and ccd.fits" << std::endl;
  }
}
