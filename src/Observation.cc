#include "../include/Observation.h"
#include "../include/Conversion.h"
#include <gsl/gsl_randist.h>

using namespace skylens;

unsigned int Observation::SUBPIX = 1;

Observation::Observation (const Telescope& tel, const sed& sky, double time, int n_exposures) :
  tel(tel), time(time), nexp(n_exposures)
{
  // construct NullLayer to connect all Leyers behind
  new NullLayer();
  
  // create up sky background layer
  // since sky is in flux/arcsec^2, we need pixelsize
  double photons_pixel = Conversion::emission2photons(sky,time,tel)*gsl_pow_2(tel.pixsize);
  new SkyFluxLayer(Conversion::photons2ADU(photons_pixel,tel.gain));
  
  // compute noise characteristic
  setNoise();
}

Observation::Observation (const Telescope& tel, double sky_mag, double time, int n_exposures) :
  tel(tel), time(time), nexp(n_exposures)
{
  // construct NullLayer to connect all Leyers behind
  new NullLayer();

  // create up sky background layer
  // since sky is in flux/arcsec^2, we need pixelsize
  double sky_flux = Conversion::mag2flux(sky_mag);
  double sky_photons = Conversion::flux2photons(sky_flux,time,tel);
  double sky_ADU = Conversion::photons2ADU(sky_photons,tel.gain)*gsl_pow_2(tel.pixsize);
  new SkyFluxLayer(sky_ADU);
  
  // compute noise characteristic
  setNoise();
}

Observation::~Observation() {
  gsl_rng_free (r);
}

void Observation::makeImage(shapelens::Image<double>& im, bool adjust) {
  if (adjust) {
    int npix_x = tel.fov_x/tel.pixsize, npix_y = tel.fov_y/tel.pixsize;
    if (im.getSize(0) != npix_x || im.getSize(0) != npix_y)
      im.resize(npix_x*npix_y);
    im.grid = shapelens::Grid(0,0,npix_x,npix_y);
  }
  Layer* front = SingleLayerStack::getInstance().begin()->second;
  shapelens::Point2D<int> P;
  double x,y;
  double subpix_dist = 1./SUBPIX;
  for (unsigned long i=0; i < im.size(); i++) {
    P = im.grid.getCoords(i);
    im(i) = 0;
    for (int s1 = 0; s1 < SUBPIX; s1++) {
      x = tel.pixsize*(P(0)+(s1+0.5)*subpix_dist); // center of subpixel
      for (int s2 = 0; s2 < SUBPIX; s2++) {
	y = tel.pixsize*(P(1)+(s2+0.5)*subpix_dist);
	//std::cout << x << "\t" << y << "\t" << front->getFlux(x,y) << std::endl;
	im(i) += front->getFlux(x,y);
      }
    }
    im(i) /= SUBPIX*SUBPIX;
    addNoise(im(i));
  }
}

void Observation::setNoise() {
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
