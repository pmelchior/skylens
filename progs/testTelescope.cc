#include <iostream>
#include <Telescope.h>

using namespace skylens;

int main() {
  // simply test telescope constructor
  double d;
  double pix;
  double rn;
  std::string nick;
  SUBARU s("B");
  d=s.getDiameter();
  pix=s.getPixelScale();
  rn=s.getReadOutNoise();
  nick=s.getName();
  std::cout << "Telescope NickName=" << nick << "\n";
  std::cout << "Telescope diameter=" << d << "\n";
  std::cout << "Telescope pixel scale=" << pix << "\n";
  std::cout << "Telescope RON=" << rn << "\n";
  return 0;
}
