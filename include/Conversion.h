#ifndef SKYLENS_CONVERION_H
#define SKYLENS_CONVERION_H

#include <sed.h>
#include "Telescope.h"

namespace skylens {
  class Conversion {
  public:
    static double mag2flux(double mag);
    static double flux2mag(double flux);
    static double flux2photons(double flux, double time, const Telescope& tel);
    static double photons2flux(double photons, double time, const Telescope& tel);
    static double emission2photons(const sed& emission, double time, const Telescope& tel);
    static double photons2ADU(double photons, double gain);
    static double ADU2photons(double ADU, double gain);
    static double ADU2flux(double ADU, double zeropoint);
    static double flux2ADU(double flux, double zeropoint);
    static double zeroPoint(double time, const Telescope& tel);
  };
} // end namespace
#endif
