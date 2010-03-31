#include "../include/Telescope.h"
#include "../include/Layer.h"
#include "../include/Helpers.h"
#include <shapelens/utils/Property.h>
#include <fstream>
#include <stdexcept>
#include <iostream>
#include <fenv.h>

namespace skylens {

  
Telescope::Telescope() : name("Dummy") , bandname("Dummy") {
  diameter = flat_acc = pixsize = fov_x = fov_y = gain = ron = 0;
}

Telescope::Telescope(std::string configfile, std::string filterfile) {
  std::string datapath = std::string(getenv("SKYLENSDATAPATH"));
  // test if filter file is valid
  std::ifstream ifs;
  test_open(ifs,datapath,filterfile);
  total = filter(filterfile,"/");
  // read config file
  readConfig(datapath,configfile);
  // set name of telescope and band
  name = split(split(configfile,'/').back(),'.').front();
  bandname = split(split(filterfile,'/').back(),'.').front();
}

void Telescope::readConfig(std::string datapath, std::string configfile) {
  std::ifstream ifs;
  test_open(ifs,datapath,configfile);

  shapelens::Property config;
  config.read(ifs);
  ifs.close();

  diameter = boost::get<shapelens::data_t>(config["DIAMETER"]);
  flat_acc = boost::get<shapelens::data_t>(config["FLAT-ACC"]);
  gain     = boost::get<shapelens::data_t>(config["GAIN"]);
  pixsize  = boost::get<shapelens::data_t>(config["PIXSIZE"]);
  fov_x    = boost::get<shapelens::data_t>(config["FOV_X"]);
  fov_y    = boost::get<shapelens::data_t>(config["FOV_Y"]);
  ron      = boost::get<shapelens::data_t>(config["RON"]);
  
  try { // either numbers: flat spectrum
    double ccd = boost::get<shapelens::data_t>(config["CCD"]);
    total *= ccd;
  } catch (boost::bad_get) { // or string: filter file
    std::string ccd = boost::get<std::string>(config["CCD"]);
    test_open(ifs,datapath,ccd);
    total *= filter(ccd,"/");
  }
  try {
    double mirror = boost::get<shapelens::data_t>(config["MIRROR"]);
    total *= mirror;
  } catch (boost::bad_get) {
    std::string mirror = boost::get<std::string>(config["MIRROR"]);
    test_open(ifs,datapath,mirror);
    total *= filter(mirror,"/");
  }
  try {
    double optics = boost::get<shapelens::data_t>(config["OPTICS"]);
    total *= optics;
  } catch (boost::bad_get) {
    std::string optics = boost::get<std::string>(config["OPTICS"]);
    test_open(ifs,datapath,optics);
    total *= filter(optics,"/");
  }

  // FIXME: PSF stuff not implemented
  // psf = PSF(path + "/psf.sif");

  // check for mask
  try {
    std::string mask = boost::get<std::string>(config["MASK"]);
    test_open(ifs,datapath,mask);
    new MaskLayer(mask);
  } catch (std::invalid_argument) {}
}

} // end namespace
