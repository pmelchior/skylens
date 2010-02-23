#include "../include/Layer.h"
#include <shapelens/utils/Interpolation.h>
#include <shapelens/utils/IO.h>
#include <libastro/constants.h>

using namespace skylens;

LensingLayer::LensingLayer(double z_, std::string angle_file) :
  // automatically creates a single instance of LayerStack
  ls(SingleLayerStack::getInstance()),
  li(SingleLensingInformation::getInstance()),
  cosmo(SingleCosmology::getInstance())
{
  Layer::z = z_;
  Layer::transparent = false;
  me = ls.insert(std::pair<double,Layer*>(z,this));

  // check if this is the first lensing layer
  li.z_first_lens = std::min(li.z_first_lens, z);

  // open file with two real-valued images of first and second component
  fitsfile* fptr = shapelens::IO::openFITSFile(angle_file);
  // read in lens parameters
  double sidelength, z_lens, z_source, h, omega, lambda;
  shapelens::IO::readFITSKeyword(fptr,"SIDEL",sidelength);
  shapelens::IO::readFITSKeyword(fptr,"ZLENS",z_lens);
  shapelens::IO::readFITSKeyword(fptr,"ZSOURCE",z_source);
  shapelens::IO::readFITSKeyword(fptr,"OMEGA",omega);
  shapelens::IO::readFITSKeyword(fptr,"LAMBDA",lambda);
  shapelens::IO::readFITSKeyword(fptr,"H",h);
  // compute rescaling factor of deflection angle
  // in the cosmology from the FITS file
  cosmology cosmo_l(omega,lambda,h);
  const constants& consts = cosmo_l.getConstants();
  double Dl, Dls, Ds, c_H0;
  // D in units [c/H0] = [cm] -> [Mpc/h]
  c_H0 = consts.get("c")/consts.get("H0")*h/consts.get("Mpc");
  Dl = cosmo_l.angularDist(0,z_lens)*c_H0;
  Ds = cosmo_l.angularDist(0,z_source)*c_H0;
  Dls = cosmo_l.angularDist(z_lens,z_source)*c_H0;
  // take out lensing efficiency factor (which is reinserted in getFlux())
  // and convert to arcsec
  scale0 = Ds/Dls * 180/M_PI * 3600;
  
  // Massimo's convention for deflection angles contains 
  // rescaling factor of sidelength/Dl
  // Do we need to apply it: Yes if keyword LENSRESC is set
  bool lensresc;
  try {
    shapelens::IO::readFITSKeyword(fptr,"LENSRESC",lensresc);
    scale0 *= sidelength/Dl;
  } catch (std::exception) {}
  
  // read in complex deflection angle field
  shapelens::IO::readFITSImage(fptr,a);

  // compute angular rescaling factor: 
  // arcsec -> pixel position in angle map
  // if the lens needs to be moved, we must set it here!
  double theta0 = (sidelength/Dl)*(180/M_PI)*3600/a.grid.getSize(0);
  std::cout << "SIDEL [arcsec]" << theta0*a.grid.getSize(0) << std::endl;
  a.grid.setWCS(shapelens::ScalarTransformation<double>(theta0));

  shapelens::IO::closeFITSFile(fptr);
}

// sum all fluxes until the next transformation layer is found
double LensingLayer::getFlux(const shapelens::Point<double>& P) const {
  double flux = 0;
  std::string type;
  LayerStack::iterator iter = me;
  iter++; // next layer

  // first lensing layer must selectively switch on source layers
  if (li.z_first_lens == z) {
    // first time only: 
    // find source layers and their distances
    if (li.Ds.size() == 0) {
      // compute c/H0: units of D
      const constants& consts = SingleCosmology::getInstance().getConstants();
      li.c_H0 = consts.get("c")/consts.get("H0")/consts.get("Mpc")*consts.get("h100");
      for (iter; iter != ls.end(); iter++) {
	type = iter->second->getType();
	if (type[0] == 'S') {
	  li.Ds[iter->first] = cosmo.angularDist(0,iter->first)*li.c_H0;
	}
      }
    }

    // raytrace thru lensed source layers: one at a time
    for (li.current_source = li.Ds.begin(); li.current_source != li.Ds.end(); li.current_source++) {
      // switch off all but one source layer: current_source
      for (std::map<double, double>::iterator source = li.Ds.begin(); source != li.Ds.end(); source++) {
	iter = ls.find(source->first);
	if (iter->first == li.current_source->first)
	  iter->second->transparent = false;
	else
	  iter->second->transparent = true;
      }

      // apply lens equation: beta = theta - Dls/Ds*alpha(theta)
      complex<float> p(P(0),P(1));
      if (!transparent)
	p -= a.interpolate(P) * scale0* float(cosmo.angularDist(z, li.current_source->first)*li.c_H0 / li.current_source->second);
      iter = me;
      iter++;
      for (iter; iter != ls.end(); iter++) { 
	type = iter->second->getType();
	flux += iter->second->getFlux(shapelens::Point<double>(real(p),imag(p)));
	// since current_source is the last and only shining
	// source layer, we can stop here
	if (type[0] == 'T' || iter->first == li.current_source->first)
	  break;
      }
    }
  }
  
  // any farther lensing layer acts normally
  else {
    // apply lens equation:
    // beta = theta - Dls/Ds*alpha(theta)
    complex<float> p(P(0),P(1));
    if (!transparent)
      p -= a.interpolate(P) * scale0 * float(cosmo.angularDist(z, li.current_source->first)*li.c_H0 / li.current_source->second); 
    for (iter; iter != ls.end(); iter++) {
      std::string type = iter->second->getType();
      flux += iter->second->getFlux(shapelens::Point<double>(real(p),imag(p)));
      // since current_source is the last and only shining
      // source layer, we can stop here
      if (type[0] == 'T' || iter->first == li.current_source->first)
	break;
    }
  }
 
  return flux;
}

std::string LensingLayer::getType() const {
  return "TL";
}

shapelens::Point<double> LensingLayer::getCenter() const {
  shapelens::Point<int> pc(a.grid.getSize(0)/2,a.grid.getSize(1)/2);
  return a.grid(a.grid.getPixel(pc));
}
