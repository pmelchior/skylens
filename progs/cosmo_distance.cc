#include <skylens/Cosmology.h>
#include <iostream>
#include <tclap/CmdLine.h>

int main(int argc, char* argv[]) {

  TCLAP::CmdLine cmd("Compute cosmological distances", ' ', "0.2");
  TCLAP::UnlabeledValueArg<double> z("z","Emission redshift",true,1,"double", cmd);
  TCLAP::ValueArg<double> o("o","o","Observation redshift",false,0,"double", cmd);
  std::vector<TCLAP::Arg*> v;
  TCLAP::SwitchArg a("a","angular","Angular diameter distance",false);
  TCLAP::SwitchArg l("l","luminosity","Luminosity distance",false);
  TCLAP::SwitchArg c("c","comoving","Comoving distance",false);
  TCLAP::SwitchArg p("p","proper","Proper distance",false);
  v.push_back(&a);
  v.push_back(&l);
  v.push_back(&c);
  v.push_back(&p);
  cmd.xorAdd(v);
  cmd.parse(argc,argv);

  skylens::Cosmology cosmo;
  // D in units [c/H0] = [cm] -> [Mpc/h]
  double c_H0 = cosmo.getc()/cosmo.getH0()*cosmo.h100/cosmo.getMpc();
  if (a.isSet())
    std::cout << cosmo.Dang(z.getValue(),o.getValue())*c_H0 << " Mpc/h" << std::endl;
  if (l.isSet())
    std::cout << cosmo.Dlum(z.getValue(),o.getValue())*c_H0 << " Mpc/h" << std::endl;
  if (c.isSet())
    std::cout << cosmo.Dcom(z.getValue(),o.getValue())*c_H0 << " Mpc/h" << std::endl;
  if (p.isSet())
    std::cout << cosmo.Dprop(z.getValue(),o.getValue())*c_H0 << " Mpc/h" << std::endl;

}
