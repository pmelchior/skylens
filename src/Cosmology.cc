#include "../include/Cosmology.h"
#include "../include/Integrator.h"
#include <shapelens/MathHelper.h>

namespace skylens {

  Cosmology::Cosmology(double Omega_m_, double Omega_l_, double h100_) : Omega_m(Omega_m_ ), Omega_l(Omega_l_), Omega_c(1. - Omega_m - Omega_l), h100(h100_) {}

  double Cosmology::a(double z) {
    return 1./(1+z);
  }

  double Cosmology::E(double a) {
    return sqrt(Omega_m*1./(a*a*a) + Omega_c*1./(a*a) + Omega_l);
  }
  double Cosmology::Dcom(double z, double z_ref) {
    Integrator<Cosmology, &Cosmology::angKernel> integrate (*this);
  return integrate(1 + z_ref, 1 + z);
  }
  double Cosmology::Dang(double z, double z_ref) {
    Integrator<Cosmology, &Cosmology::angKernel> integrate (*this);
    double rk = sqrt(fabs(Omega_c));
    double d = integrate (z_ref + 1, z + 1);
    if (rk*d>0.01) {
      if (Omega_c > 0) 
	d=sinh(rk*d)/rk;
      if (Omega_c < 0) 
	d=sin(rk*d)/rk;
    }
    return d/(1+z);
  }
  double Cosmology::Dprop(double z, double z_ref) {
    Integrator<Cosmology, &Cosmology::propKernel> integrate (*this);
    return integrate(1 + z_ref, 1 + z);
  }
  double Cosmology::Dlum(double z, double z_ref) {
    return shapelens::pow2((1 + z)/(1 + z_ref))*Dang(z, z_ref);
  }
  double Cosmology::Dlens(double z, double z_lens) {
    return Dang(z_lens)*Dang(z, z_lens)/Dang(z);
  }
  double Cosmology::angKernel(double x) {
    return 1./E(1./x);
  }
  double Cosmology::propKernel(double x) {
    return DpropDz(x-1);
  }
  double Cosmology::DpropDz(double z) {
    double a_ = a(z);
    return a_/E(a_);
  }

  double Cosmology::getc() { return 2.99792458e+10; } // in cm/s
  double Cosmology::getH0() { return 100*1e5/getMpc()*h100;} // in s
  double Cosmology::getMpc() { return 3.0856775806e24;} // in cm
  
} // end namespace
