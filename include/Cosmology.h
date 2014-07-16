#ifndef SKYLENS_COSMOLOGY_H
#define SKYLENS_COSMOLOGY_H

namespace skylens {

  /// Minimal cosmology class.
  /// Basic cosmological calculations such as angular diameter distances.
  /// 
  /// Derived from Matthias Bartelmann's libastro 
  class Cosmology {
  public:
    /// Default constructor.
    /// Vanilla LambdaCDM universe: radiation is assumed to be zero and 
    /// Dark Energy constant with w = -1, but curvature is arbitrary.
    Cosmology(double Omega_m=0.3, double Omega_l=0.7, double h100=0.7);
    double a(double z);
    double E(double a);
    double Dcom(double z, double z_ref=0);
    double Dang(double z, double z_ref=0);
    double Dprop(double z, double z_ref=0);
    double Dlum(double z, double z_ref=0);
    double Dlens(double z, double z_lens);
    double getc();
    double getH0();
    double getMpc();
    double Omega_m, Omega_l, Omega_c, h100;
  protected:
    double angKernel(double x);
    double propKernel(double x);
    double DpropDz(double z);
  };

} // end namespace

#endif
