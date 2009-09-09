#include "../include/Layer.h"
#include <shapelens/utils/Interpolation.h>
#include <shapelens/utils/IO.h>

using namespace skylens;

LensingLayer::LensingLayer(double z, std::string angle_file) :
  // automatically creates a single instance of LayerStack
  ls(SingleLayerStack::getInstance()),
  li(SingleLensingInformation::getInstance())
{
  Layer::z = z;
  Layer::transparent = false;
  me = ls.insert(std::pair<double,Layer*>(z,this));

  // check if this is the first lensing layer
  li.z_first_lens = std::min(li.z_first_lens, z);

  // open file with two real-valued images of first and second component
  fitsfile* fptr = shapelens::IO::openFITSFile(angle_file);
  double sidelength, z_lens, z_source;
  shapelens::IO::readFITSKeyword(fptr,"SIDEL",sidelength);
  shapelens::IO::readFITSKeyword(fptr,"ZLENS",z_lens);
  shapelens::IO::readFITSKeyword(fptr,"ZSOURCE",z_source);
  // FIXME: need to compute lensing parameters from these values
  shapelens::Image<float> component;
  shapelens::IO::readFITSImage(fptr,component);
  a.resize(component.size());
  for(unsigned long i=0; i < component.size(); i++)
    a(i) = component(i);
  // move on to 2nd component
  int status = 0, hdutype;
  fits_movrel_hdu(fptr, 1, &hdutype, &status);
  shapelens::IO::readFITSImage(fptr,component);
  for(unsigned long i=0; i < component.size(); i++)
    a(i) += complex<float>(0,component(i));
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
    LayerStack::iterator source;
    // first time only: define first lensed source layer 
    if (li.z_first_lensed_source == 0) {
      for (source = iter; source != ls.end(); source++) {
	type = source->second->getType();
	if (type[0] == 'S') {
	  li.z_first_lensed_source = li.z_current_source = source->first;
	  break; // layer are sorted: we're done
	}
      }
    }

    // raytrace thru lensed source layers: one at a time
    for (source = ls.find(li.z_first_lensed_source); source != ls.end(); source++) {
      type = source->second->getType();
      // iterate only over source layers
      if (type[0] == 'S') {
	// switch off all but one source layer: current_source
	for (LayerStack::iterator siter = ls.find(li.z_first_lensed_source); siter != ls.end(); siter++) {
	  type = siter->second->getType();
	  if (type[0] == 'S') {
	    if (siter->first == li.z_current_source)
	      siter->second->transparent = false;
	    else
	      siter->second->transparent = true;
	  }
	}
	// apply lens equation
	complex<float> p(P(0),P(1));
	if (!transparent)
	  p -= a.interpolate(P); // FIXME: need factor from current_source
	iter = me;
	iter++;
	for (iter; iter != ls.end(); iter++) { 
	  type = iter->second->getType();
	  flux += iter->second->getFlux(shapelens::Point<double>(real(p),imag(p)));
	  // since current_source is the last and only shining
	  // source layer, we can stop here
	  if (type[0] == 'T' || iter->first == li.z_current_source)
	    break;
	}

	// find the next current_source
	iter = ls.find(li.z_current_source);
	iter++;
	for (iter; iter != ls.end(); iter++) {
	  type = iter->second->getType();
	  if (type[0] == 'S')
	    li.z_current_source = iter->first;
	}
    }
    }
  }
  
  // any farther lensing layer acts normally
  else {
    // apply lens equation
    complex<float> p(P(0),P(1));
    if (!transparent)
      p -= a.interpolate(P); // FIXME: need factor from current_source
    for (iter; iter != ls.end(); iter++) {
      std::string type = iter->second->getType();
      flux += iter->second->getFlux(shapelens::Point<double>(real(p),imag(p)));
      // since current_source is the last and only shining
      // source layer, we can stop here
      if (type[0] == 'T' || iter->first == li.z_current_source)
	break;
    }
  }
 
  return flux;
}

std::string LensingLayer::getType() const {
  return "TL";
}
