#include "../include/Conversion.h"

using namespace skylens;

double Conversion::mag2flux(double mag) {
  return pow(10.,-0.4*(mag + 48.6));
}

double Conversion::flux2mag(double flux) {
  return -2.5*log10(flux) - 48.6;
}

double Conversion::flux2photons(double flux, double time, const Telescope& tel, const filter& total) {
  // telescope diameter in meters, needs to be in cm
  return M_PI*gsl_pow_2(tel.diameter*100)*time*flux*total.getQe()/(4*6.62606885e-27);
}

double Conversion::photons2flux(double photons, double time, const Telescope& tel, const filter& total) {
  return (photons*4*6.62606885e-27)/(M_PI*gsl_pow_2(tel.diameter*100)*time*total.getQe());
}

double Conversion::emission2photons(const sed& emission, double time, const Telescope& tel, const filter& total) {
  sed em_filtered = emission;
  em_filtered *= total;
  return M_PI*gsl_pow_2(tel.diameter*100)*time*em_filtered.getNorm()/(4*6.62606885e-27);
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

double Conversion::zeroPoint(const Telescope& tel, const filter& total, double time) {
  return 2.5*log10(0.43035e-3*gsl_pow_2(tel.diameter*100)*time*total.getQe()/tel.gain) + 25;
}
