#include "../include/Conversion.h"
#include <shapelens/MathHelper.h>

using namespace skylens;
using shapelens::pow2;

double Conversion::mag2flux(double mag) {
  return pow(10.,-0.4*(mag + 48.6));
}

double Conversion::flux2mag(double flux) {
  return -2.5*log10(flux) - 48.6;
}

double Conversion::flux2photons(double flux, double time, const Telescope& tel, const Filter& total) {
  // telescope diameter in meters, needs to be in cm
  return M_PI*pow2(tel.diameter*100)*time*flux*total.computeNorm()/(4*6.62606885e-27);
}

double Conversion::photons2flux(double photons, double time, const Telescope& tel, const Filter& total) {
  return (photons*4*6.62606885e-27)/(M_PI*pow2(tel.diameter*100)*time*total.computeNorm());
}

double Conversion::emission2photons(const SED& emission, double time, const Telescope& tel, const Filter& total) {
  SED em_filtered(emission);
  em_filtered *= total;
  return M_PI*pow2(tel.diameter*100)*time*em_filtered.computeNorm()/(4*6.62606885e-27);
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

double Conversion::zeroPoint(const Telescope& tel, const Filter& total, double time) {
  return 2.5*log10(0.43035e-3*pow2(tel.diameter*100)*time*total.computeNorm()/tel.gain) + 25;
}
