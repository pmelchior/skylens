#include "../include/Observation.h"
#include "../include/Conversion.h"

using namespace skylens;

Observation::Observation (const Telescope& tel, const sed& sky, double time, int n_exposures) :
  tel(tel), time(time), sky(sky), nexp(n_exposures)
{
  // create up sky background layer
  new SkyFluxLayer(Conversion::photons2ADU(Conversion::emission2photons(sky,time,tel),tel.gain));

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
