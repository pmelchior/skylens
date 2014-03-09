#include "../include/Conversion.h"
#include <shapelens/MathHelper.h>

using namespace skylens;
using shapelens::pow_int;

double Conversion::mag2flux(double mag) {
  return pow(10.,-0.4*(mag + 48.6));
}

double Conversion::flux2mag(double flux) {
  return -2.5*log10(flux) - 48.6;
}

double Conversion::flux2photons(double flux, double time, const Telescope& tel, const astro::filter& total) {
  // telescope diameter in meters, needs to be in cm
  return M_PI*pow_int(tel.diameter*100, 2)*time*flux*total.getQe()/(4*6.62606885e-27);
}

double Conversion::photons2flux(double photons, double time, const Telescope& tel, const astro::filter& total) {
  return (photons*4*6.62606885e-27)/(M_PI*pow_int(tel.diameter*100, 2)*time*total.getQe());
}

double Conversion::emission2photons(const astro::sed& emission, double time, const Telescope& tel, const astro::filter& total) {
  astro::sed em_filtered = emission;
  em_filtered *= total;
  return M_PI*pow_int(tel.diameter*100, 2)*time*em_filtered.getNorm()/(4*6.62606885e-27);
}

double Conversion::photons2ADU(double photons, double gain) {
  return photons/gain;
}

double Conversion::ADU2photons(double ADU, double gain) {
  return ADU*gain;
}

double Conversion::ADU2flux(double ADU, double zeropoint) {
  return ADU*mag2flux(zeropoint);
}

double Conversion::flux2ADU(double flux, double zeropoint) {
  return flux/mag2flux(zeropoint);
}

double Conversion::zeroPoint(const Telescope& tel, const astro::filter& total, double time) {
  return 2.5*log10(0.43035e-3*pow_int(tel.diameter*100, 2)*time*total.getQe()/tel.gain) + 25;
}
