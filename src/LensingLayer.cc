#include "../include/Layer.h"
#include <shapelens/utils/Interpolation.h>
#include <shapelens/utils/IO.h>

using namespace skylens;

LensingLayer::LensingLayer(double z, std::string angle_file) :
  // automatically creates a single instance of LayerStack
  ls(SingleLayerStack::getInstance()),
  li(SingleLensingInformation::getInstance()),
  cosmo(SingleCosmology::getInstance())
{
  Layer::z = z;
  Layer::transparent = false;
  me = ls.insert(std::pair<double,Layer*>(z,this));

  // check if this is the first lensing layer
  li.z_first_lens = std::min(li.z_first_lens, z);

  // open file with two real-valued images of first and second component
  fitsfile* fptr = shapelens::IO::openFITSFile(angle_file);
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
  double scale, Dl, Dls, Ds;
  // FIXME: D in units [c/H0]!!!
  Dl = cosmo_l.angularDist(0,z_lens);
  Ds = cosmo_l.angularDist(0,z_source);
  Dls = cosmo_l.angularDist(z_lens,z_source);
  // FIXME: is xi0 = sidelength or xi0 = sidelength*h ???
  scale = sidelength*Ds/(Dl*Dls);
  // std::cout << "scale = " <<scale << std::endl;
  
  // file contains two real-valued images of first and second component
  shapelens::Image<float> component;
  shapelens::IO::readFITSImage(fptr,component);
  a.resize(component.size());
  a.grid = component.grid;

  // compute angular rescaling factor: 
  // arcsec -> pixel position in angle map
  double theta0 = sidelength/(Dl*a.grid.getSize(0)); // FIXME: sidelength
  a.grid.setWCS(shapelens::ScalarTransformation<double>(theta0));
  // std::cout << "theta0 = " << theta0 << std::endl;

  // copy 1st component
  for(unsigned long i=0; i < component.size(); i++)
    a(i) = scale*component(i);

  // move on to 2nd component
  int status = 0, hdutype;
  fits_movrel_hdu(fptr, 1, &hdutype, &status);
  shapelens::IO::readFITSImage(fptr,component);
  for(unsigned long i=0; i < component.size(); i++)
    a(i) += complex<float>(0,scale*component(i));
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
      for (iter; iter != ls.end(); iter++) {
	type = iter->second->getType();
	if (type[0] == 'S') {
	  li.Ds[iter->first] = cosmo.angularDist(0,iter->first);
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
      //std::cout << P << " -> " << p << " - ";
      if (!transparent)
	p -= a.interpolate(P) * (cosmo.angularDist(z, li.current_source->first) / li.current_source->second);
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
      p -= a.interpolate(P) * (cosmo.angularDist(z, li.current_source->first) / li.current_source->second); 
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
