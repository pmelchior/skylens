#include <iostream>
#include <Telescope.h>

using namespace skylens;

void printInfo(const Telescope& t) {
  std::cout << "Telescope name\t" << t.name << std::endl;
  std::cout << "Telescope diameter\t" << t.diameter << "\t[meter]" << std::endl;
  std::cout << "Telescope pixel scale\t" << t.pixsize << "\t[arcsec/pixel]" << std::endl;
  std::cout << "Telescope FOV\t" << t.fov_x << "x" << t.fov_y << "\t[arcsec^2]" << std::endl;
  std::cout << "Telescope RON\t" << t.ron << "\t[???]" << std::endl;
  std::cout << "Telescope gain\t" << t.gain << "\t[???]" << std::endl;
  std::cout << "Filter band\t" << t.band << std::endl;
  std::cout << "Filter eff. wavelength\t" << t.total.lambdaEff() << "\t[Angstrom]" << std::endl;
}


int main() {
  // simply test telescope constructor
  SUBARU s("B");
  printInfo(s);
  return 0;
}
