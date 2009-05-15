#include "../include/Observation.h"
#include "../include/Conversion.h"

using namespace skylens;

Observation::Observation (const Telescope& tel, const sed& sky, double time, int n_exposures) :
  tel(tel), time(time), sky(sky), nexp(n_exposures)
{
  // create up sky background layer
  // since sky is in flux/arcsec^2, we need pixelsize
  double photons_pixel = Conversion::emission2photons(sky,time,tel)*gsl_pow_2(tel.pixsize);
  new SkyFluxLayer(Conversion::photons2ADU(photons_pixel,tel.gain));

  // alternative: sky magnitude is given ...
  // double sky_mag = ???;
//   double sky_flux = Conversion::mag2flux(sky_mag);
//   double sky_photons = Conversion::flux2photons(sky_flux,time,tel);
//   double sky_ADU = Conversion::photons2ADU(sky_photons,tel.gain)*gsl_pow_2(tel.pixsize);
//   new SkyFluxLayer(sky_ADU);

  // set up noise layer
  // eq. (31) in Meneghetti et al. (2008)
  double ron = nexp*gsl_pow_2(tel.ron/tel.gain);
  double f = 0; // should contain residuals after flat-field subtraction
                // where could we get this from???
  double flat_field = f + gsl_pow_2(tel.flat_acc/nexp);
  new NoiseLayer(ron,flat_field);
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
  for (unsigned long i=0; i < im.size(); i++) {
    P = im.grid.getCoords(i);
    x = tel.pixsize*(P(0)+0.5); // centered pixellation
    y = tel.pixsize*(P(1)+0.5);
    im(i) = front->getFlux(x,y);
  }
}
