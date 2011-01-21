#include "../include/Layer.h"
#include <shapelens/utils/Interpolation.h>
#include <shapelens/utils/IO.h>
#include <libastro/constants.h>

using namespace skylens;
using std::complex;

LensingLayer::LensingLayer(double z_, std::string angle_file, const shapelens::Point<double>& center) :
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
  astro::cosmology cosmo_l(omega,lambda,h);
  const astro::constants& consts = cosmo_l.getConstants();
  double Dl, Dls, Ds, c_H0;
  // D in units [c/H0] = [cm] -> [Mpc/h]
  c_H0 = consts.get_lightspeed()/consts.get_Hubble()*h/consts.get_Megaparsec();
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
    if (lensresc)
      scale0 *= sidelength/Dl;
  } catch (std::exception) {}
  
  // read in complex deflection angle field
  shapelens::IO::readFITSImage(fptr,a);

  // compute angular rescaling factor: 
  // arcsec -> pixel position in angle map
  // if the lens needs to be moved, we must set it here!
  theta0 = (sidelength/Dl)*(180/M_PI)*3600/a.grid.getSize(0);
  shapelens::ScalarTransformation S(theta0);

  if (center(0) != 0 || center(1) != 0) {
    shapelens::Point<double> center_image(0.5*a.grid.getSize(0),0.5*a.grid.getSize(1));
    shapelens::ShiftTransformation Z(-center_image);
    shapelens::ShiftTransformation ZF(center);
    S *= ZF;
    Z *= S;
    a.grid.setWCS(Z);
  } else
    a.grid.setWCS(S);

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
    // find source layers and compute their distances to redshift 0 and z_l
    if (li.Ds.size() == 0) {
      li.Dls.insert(std::pair<double, std::map<double,double> >(z,std::map<double,double>()));
      // compute c/H0: units of D
      const astro::constants& consts = SingleCosmology::getInstance().getConstants();
      li.c_H0 = consts.get_lightspeed()/consts.get_Hubble()*consts.get_h100()/consts.get_Megaparsec();
      for (iter; iter != ls.end(); iter++) {
	type = iter->second->getType();
	if (type[0] == 'S') {
	  li.Ds[iter->first] = cosmo.angularDist(0,iter->first)*li.c_H0;
	  li.Dls[z][iter->first] = cosmo.angularDist(z,iter->first)*li.c_H0;
	}
      }
    }

    // raytrace thru lensed source layers: one at a time
    for (li.current_source = li.Ds.begin(); li.current_source != li.Ds.end(); li.current_source++) {
      if (li.Ds.size() > 1) {
	// switch off all but one source layer: current_source
	for (std::map<double, double>::iterator source = li.Ds.begin(); source != li.Ds.end(); source++) {
	  iter = ls.find(source->first);
	  if (iter->first == li.current_source->first)
	    iter->second->transparent = false;
	  else
	    iter->second->transparent = true;
	}
      }

      // apply lens equation: beta = theta - Dls/Ds*alpha(theta)
      complex<float> p(P(0),P(1));
      if (!transparent)
	p -= a.interpolate(P) * scale0* float(li.Dls[z][li.current_source->first] / li.current_source->second);
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
    // calculate Dls for this layer if not not in li yet
    if (li.Dls.find(z) == li.Dls.end()) {
      li.Dls.insert(std::pair<double, std::map<double,double> >(z,std::map<double,double>()));
      for (iter; iter != ls.end(); iter++) {
	type = iter->second->getType();
	if (type[0] == 'S') {
	  li.Dls[z][iter->first] = cosmo.angularDist(z,iter->first)*li.c_H0;
	}
      }
    }
    // apply lens equation:
    // beta = theta - Dls/Ds*alpha(theta)
    complex<float> p(P(0),P(1));
    if (!transparent)
      p -= a.interpolate(P) * scale0 * float(li.Dls[z][li.current_source->first] / li.current_source->second); 
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

std::map<shapelens::Point<double>, shapelens::Point<double> > LensingLayer::findCriticalPoints(double zs) {
  const astro::constants& consts = cosmo.getConstants();
  double c_H0 = consts.get_lightspeed()/consts.get_Hubble()*consts.get_h100()/consts.get_Megaparsec();
  double D_ls = cosmo.angularDist(z,zs)*c_H0;
  double D_s = cosmo.angularDist(0,zs)*c_H0;
  std::map<shapelens::Point<double>, shapelens::Point<double> > cpoints;
  for (long i = 2; i < a.grid.getSize(0) - 2; i++) {
    for (long j = 2; j < a.grid.getSize(1) - 2; j++) {
      double phixx = (-real(a(i+2,j)) + 8.0*real(a(i+1,j)) 
		      - 8.0*real(a(i-1,j)) +real(a(i-2,j)))/(12*theta0);
      double phiyy = (-imag(a(i,j+2)) + 8.0*imag(a(i,j+1))
		      - 8.0*imag(a(i,j-1)) + imag(a(i,j-2)))/(12*theta0);
      double phixy = (-real(a(i,j+2)) + 8.0*real(a(i,j+1))
		      - 8.0*real(a(i,j-1)) + real(a(i,j-2)))/(12*theta0);
      double phiyx = (-imag(a(i+2,j)) + 8.0*imag(a(i+1,j))
		      - 8.0*imag(a(i-1,j)) + imag(a(i-2,j)))/(12*theta0);
      double kappa = scale0 * float(D_ls / D_s) * 0.5*(phixx + phiyy);
      double gamma1 = scale0 * float(D_ls / D_s)* 0.5*(phixx - phiyy);
      double gamma2 = scale0 * float(D_ls / D_s)* phixy;
      double gamma = sqrt(gamma1*gamma1+gamma2*gamma2);
      double jacdet = (1.0-kappa)*(1.0-kappa) - gamma*gamma;
      if (fabs(jacdet) < 1e-2) {
	shapelens::Point<double> critical(a.grid(a.grid.getPixel(shapelens::Point<int>(i,j))));
	complex<float> alpha = a(i,j) * scale0 * float(D_ls / D_s);
	shapelens::Point<double> caustic(critical(0)-real(alpha), critical(1) - imag(alpha));
	cpoints.insert(std::pair<shapelens::Point<double>, shapelens::Point<double> >(critical, caustic));
      }
    }
  }
  return cpoints;
}
