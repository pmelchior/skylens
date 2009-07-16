#include "../include/Observation.h"
#include "../include/Conversion.h"
#include <gsl/gsl_randist.h>

using namespace skylens;

Observation::Observation (const Telescope& tel, double exptime) : tel(tel), time(exptime), ron(0), flat_field(0) {
  // construct NullLayer to connect all Layers behind
  new NullLayer();
}

Observation::~Observation() {
  gsl_rng_free (r);
}

void Observation::makeImage(shapelens::Image<double>& im, bool adjust) {
  if (adjust) {
    int npix_x = tel.fov_x/tel.pixsize, npix_y = tel.fov_y/tel.pixsize;
    if (im.getSize(0) != npix_x || im.getSize(0) != npix_y)
      im.resize(npix_x*npix_y);
    im.grid.setSize(0,0,npix_x,npix_y);
    im.grid.apply(shapelens::ScalarTransformation<double>(tel.pixsize));
  }
  Layer* front = SingleLayerStack::getInstance().begin()->second;
  shapelens::Point<double> P;
  for (unsigned long i=0; i < im.size(); i++) {
    P = im.grid(i); // assumes proper WCS of im
    im(i) = front->getFlux(P);
    addNoise(im(i));
  }
}

// create sky background layer
void Observation::createSkyFluxLayer(const sed& sky) {
  // since sky is in flux/arcsec^2, we need pixelsize
  double photons_pixel = Conversion::emission2photons(sky,time,tel,tel.total)*gsl_pow_2(tel.pixsize);
  new SkyFluxLayer(Conversion::photons2ADU(photons_pixel,tel.gain));
}

// create sky background layer
void Observation::createSkyFluxLayer(double sky_mag) {
  double sky_flux = Conversion::mag2flux(sky_mag);
  double sky_photons = Conversion::flux2photons(sky_flux,time,tel,tel.total);
  double sky_ADU = Conversion::photons2ADU(sky_photons,tel.gain)*gsl_pow_2(tel.pixsize);
  new SkyFluxLayer(sky_ADU);
}

void Observation::computeTransmittance(const filter& atmosphere, double airmass) {
  // compute total filter curve, including air mass extinction
  total_air = tel.total;
  filter air_transmittance = atmosphere;
  arr1d<float>& curve =  air_transmittance.getCurve();
  for (unsigned int i=0; i < curve.size(); i++)
    curve[i] = pow(10.,-0.4*(airmass)*curve[i]);
  total_air *= air_transmittance;
}

void Observation::computeTransmittance(double extinction, double airmass) {
  // compute total filter curve, including air mass extinction
  total_air = tel.total;
  total_air *= pow(10.,-0.4*(airmass)*extinction);
}
    

void Observation::setNoise(int nexp) {
  // set up RNG
  const gsl_rng_type * T;
  gsl_rng_env_setup(); // read env variables to set seed and RNG type
  T = gsl_rng_default;
  r = gsl_rng_alloc(T);

  // set up noise quantities
  // eq. (31) in Meneghetti et al. (2008)
  ron = nexp*gsl_pow_2(tel.ron/tel.gain);
  double f = 0; // should contain residuals after flat-field subtraction
                // where could we get this from???
  flat_field = f + gsl_pow_2(tel.flat_acc/nexp);
}

void Observation::addNoise(double& flux) {
  flux += gsl_ran_gaussian_ziggurat (r,sqrt(ron + fabs(flux) + flat_field*flux*flux));
}

const filter& Observation::getTotalTransmittance() const {
  return total_air;
}
