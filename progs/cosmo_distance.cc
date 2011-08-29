#include <astro/cosmology.h>
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

  astro::cosmology cosmo;
  const astro::constants& consts = cosmo.getConstants();
  // D in units [c/H0] = [cm] -> [Mpc/h]
  double h = 0.7;
  double c_H0 = consts.get_lightspeed()/consts.get_Hubble()*h/consts.get_Megaparsec();
  if (a.isSet())
    std::cout << cosmo.angularDist(o.getValue(),z.getValue())*c_H0 << " Mpc/h" << std::endl;
  if (l.isSet())
    std::cout << cosmo.luminosityDist(o.getValue(),z.getValue())*c_H0 << " Mpc/h" << std::endl;
  if (c.isSet())
    std::cout << cosmo.comovingDist(o.getValue(),z.getValue())*c_H0 << " Mpc/h" << std::endl;
  if (p.isSet())
    std::cout << cosmo.properDist(o.getValue(),z.getValue())*c_H0 << " Mpc/h" << std::endl;

}
