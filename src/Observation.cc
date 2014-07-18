#include "../include/Observation.h"
#include "../include/Conversion.h"
#include "../include/RNG.h"
#include <shapelens/MathHelper.h>

namespace skylens {
  using shapelens::pow_int;

  Observation::Observation (const Telescope& tel_, double exptime) : tel(tel_), time(exptime), ron(0), flat_field(0), hasNoise(false), SUBPIXEL(1) {
    // construct NullLayer to connect all Layers behind
    new NullLayer();
  }

  void Observation::makeImage(shapelens::Image<float>& im, const shapelens::Point<double>* center) const {
    int npix_x = tel.fov_x/tel.pixsize, npix_y = tel.fov_y/tel.pixsize;
    if (im.getSize(0) != npix_x || im.getSize(1) != npix_y)
      im.resize(npix_x*npix_y);
    im.grid.setSize(0,0,npix_x,npix_y);
    shapelens::ScalarTransformation S(tel.pixsize);
    if (center != NULL) {
      shapelens::Point<double> center_image(-0.5*npix_x, -0.5*npix_y);
      shapelens::ShiftTransformation Z(center_image);
      shapelens::ShiftTransformation ZF(*center);
      S *= ZF;
      Z *= S;
      im.grid.setWCS(Z);
    } else   
      im.grid.setWCS(S);

    Layer* front = SingleLayerStack::getInstance().begin()->second;
    shapelens::Point<double> P, P_;
    RNG& rng = Singleton<RNG>::getInstance();
    double offset = 1./SUBPIXEL;  // regular subpixel shift
    const gsl_rng* r = rng.getRNG();
    for (unsigned long i=0; i < im.size(); i++) {
      P = im.grid(i);  // assumes proper WCS of im
      im(i) = 0;       // resets prior content of im !!!
      // regular subpixel sampling
      for (int n1 = 0; n1 < SUBPIXEL; n1++) {
	P_(0) = P(0) + ((0.5+n1)*offset - 0.5)*tel.pixsize;
	for (int n2 = 0; n2 < SUBPIXEL; n2++) {
	  P_(1) = P(1) + ((0.5+n2)*offset - 0.5)*tel.pixsize;
	  im(i) += front->getFlux(P_);
	}
      }
      // adding noise if demanded
      if (hasNoise)
	addNoise(r,im(i));
    }
  }

  // create sky background layer
  void Observation::createSkyFluxLayer(const SED& sky) {
    // since sky is in flux/arcsec^2, we need pixelsize
    double photons_pixel = Conversion::emission2photons(sky,time,tel,tel.total)*pow_int(tel.pixsize, 2);
    new SkyFluxLayer(Conversion::photons2ADU(photons_pixel,tel.gain));
  }

  // create sky background layer
  void Observation::createSkyFluxLayer(double sky_mag) {
    double sky_flux = Conversion::mag2flux(sky_mag);
    double sky_photons = Conversion::flux2photons(sky_flux,time,tel,tel.total);
    double sky_ADU = Conversion::photons2ADU(sky_photons,tel.gain)*pow_int(tel.pixsize, 2);
    new SkyFluxLayer(sky_ADU);
  }

  void Observation::computeTransmittance(const Filter& atmosphere, double airmass) {
    // compute total filter curve, including air mass extinction
    total_air = tel.total;
    Filter air_transmittance(atmosphere);
    for (Filter::iterator iter=air_transmittance.begin(); iter != air_transmittance.end(); iter++)
      iter->second = pow(10., -0.4*airmass*iter->second);
    total_air *= air_transmittance;
  }

  void Observation::computeTransmittance(double extinction, double airmass) {
    // compute total filter curve, including air mass extinction
    total_air = tel.total;
    total_air *= pow(10., -0.4*airmass*extinction);
  }
    

  void Observation::setNoise(int nexp) {
    hasNoise = true;
    // set up noise quantities
    // eq. (31) in Meneghetti et al. (2008)
    ron = nexp*pow_int(tel.ron/tel.gain, 2);
    double f = 0; // should contain residuals after flat-field subtraction
    // where could we get this from???
    flat_field = f + pow_int(tel.flat_acc/nexp, 2);
  }

  void Observation::addNoise(const gsl_rng* r, float& flux) const {
    flux += gsl_ran_gaussian_ziggurat (r,sqrt(ron + fabs(flux) + flat_field*flux*flux));
  }

  const Filter& Observation::getTotalTransmittance() const {
    return total_air;
  }

} // end namespace
