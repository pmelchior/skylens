#include "../include/Layer.h"
#include <shapelens/FITS.h>

#ifdef HAS_OpenMP
#include <omp.h>
#endif

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

#pragma omp critical
  {
  // check if this is the first lensing layer
  li.z_first_lens = std::min(li.z_first_lens, z);
  // compute c/H0: units of D
  li.c_H0 = cosmo.getc()/cosmo.getH0()*cosmo.h100/cosmo.getMpc();  
  }

  // open file with two real-valued images of first and second component
  fitsfile* fptr = shapelens::FITS::openFile(angle_file);
  // read in lens parameters
  double sidelength, z_lens, z_source, h, omega, lambda;
  shapelens::FITS::readKeyword(fptr,"SIDEL",sidelength);
  shapelens::FITS::readKeyword(fptr,"ZLENS",z_lens);
  shapelens::FITS::readKeyword(fptr,"ZSOURCE",z_source);
  shapelens::FITS::readKeyword(fptr,"OMEGA",omega);
  shapelens::FITS::readKeyword(fptr,"LAMBDA",lambda);
  shapelens::FITS::readKeyword(fptr,"H",h);
  // compute rescaling factor of deflection angle
  // in the cosmology from the FITS file
  Cosmology cosmo_l(omega,lambda,h);
  double Dl, Dls, Ds, c_H0;
  // D in units [c/H0] = [cm] -> [Mpc/h]
  c_H0 = cosmo_l.getc()/cosmo_l.getH0()*h/cosmo_l.getMpc();
  Dl = cosmo_l.Dang(z_lens)*c_H0;
  Ds = cosmo_l.Dang(z_source)*c_H0;
  Dls = cosmo_l.Dang(z_source,z_lens)*c_H0;
  // take out lensing efficiency factor (which is reinserted in getFlux())
  // and convert to arcsec
  scale0 = Ds/Dls * 180/M_PI * 3600;
  
  // Massimo's convention for deflection angles contains 
  // rescaling factor of sidelength/Dl
  // Do we need to apply it: Yes if keyword LENSRESC is set
  bool lensresc;
  try {
    shapelens::FITS::readKeyword(fptr,"LENSRESC",lensresc);
    if (lensresc)
      scale0 *= sidelength/Dl;
  } catch (std::exception) {}
  
  // read in complex deflection angle field
  Image<float> a1;
  shapelens::FITS::readImage(fptr, a1);
  a.resize(a1.size());
  a.grid = a1.grid;
  a.realPart() = a1;
  shapelens::FITS::moveToExtension(fptr, 2);
  shapelens::FITS::readImage(fptr, a1);
  a.imagPart() = a1;

  // compute angular rescaling factor: 
  // arcsec -> pixel position in angle map
  // if the lens needs to be moved, we must set it here!
  theta0 = (sidelength/Dl)*(180/M_PI)*3600/a.grid.getSize(0);
  shapelens::ScalarTransformation S(theta0);

  if (center(0) != 0 || center(1) != 0) {
    shapelens::Point<double> center_image(-0.5*a.grid.getSize(0),-0.5*a.grid.getSize(1));
    shapelens::ShiftTransformation Z(center_image);
    shapelens::ShiftTransformation ZF(center);
    S *= ZF;
    Z *= S;
    a.grid.setWCS(Z);
  } else
    a.grid.setWCS(S);

  shapelens::FITS::closeFile(fptr);
  
}

// iter points to current lensing layer
void LensingLayer::setDistances(const LayerStack::iterator& iter_) const {
  LayerStack::iterator iter = iter_;
  iter++;
  std::string type;
  li.Dls.insert(std::pair<double, std::map<double,double> >(z, std::map<double,double>()));
  for (iter; iter != ls.end(); iter++) {
    type = iter->second->getType();
    if (type[0] == 'S') {
      li.Ds[iter->first] = cosmo.Dang(iter->first)*li.c_H0;
      li.Dls[iter_->first][iter->first] = cosmo.Dang(iter->first, iter_->first)*li.c_H0;
    }
    else if (type == "TL") {
      // recursively call setDistances
      setDistances(iter);
    }
  }
}

// sum all fluxes until the next transformation layer is found
double LensingLayer::getFlux(const shapelens::Point<double>& P, double* z_) const {
  double flux = 0;
  if (z_ != NULL)
    if (*z_ < z)
      return flux;

  LayerStack::iterator iter = me;
  iter++;
  std::string type;

  // first lensing layer must selectively switch on source layers
  if (li.z_first_lens == z) {
    if (z_ != NULL)
      throw std::invalid_argument("LensingLayer: getFlux() source redshift set before first layer!");

    // first time only: 
    // find source layers and compute their distances to redshift 0 and z_l
#pragma omp critical
    if (li.Ds.size() == 0) {
      setDistances(me);
    }

    // raytrace thru lensed source layers: one at a time
    for (std::map<double, double>::iterator source_iter = li.Ds.begin(); source_iter != li.Ds.end(); source_iter++) {
      // active source layer 
      double zs = source_iter->first;

      // apply lens equation: beta = theta - Dls/Ds*alpha(theta)
      complex<float> p(P(0),P(1));
      if (!transparent)
	p -= a.interpolate(P) * scale0* float(li.Dls[z][zs] / source_iter->second);
      for (iter; iter != ls.end(); iter++) { 
	type = iter->second->getType();
	flux += iter->second->getFlux(shapelens::Point<double>(real(p),imag(p)), &zs);
	// since source_iter is the last and only shining
	// source layer, we can stop here
	if (type[0] == 'T' || iter->first == zs)
	  break;
      }
    }
  }
  
  // any farther lensing layer just computes the deflection angle for the
  // specified source redshify
  else {
    // apply lens equation:
    // beta = theta - Dls/Ds*alpha(theta)
    complex<float> p(P(0),P(1));
    if (!transparent)
      p -= a.interpolate(P) * scale0 * float(li.Dls[z][*z_] / li.Ds[*z_]);

    for (iter; iter != ls.end(); iter++) {
      std::string type = iter->second->getType();
      flux += iter->second->getFlux(shapelens::Point<double>(real(p),imag(p)), z_);
      // since given source is the last and only shining
      // source layer, we can stop here
      if (type[0] == 'T' || iter->first == *z_)
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
  double c_H0 = cosmo.getc()/cosmo.getH0()*cosmo.h100/cosmo.getMpc(); 
  double D_ls = cosmo.Dang(zs,z)*c_H0;
  double D_s = cosmo.Dang(zs)*c_H0;
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
