#ifndef SKYLENS_CONVERION_H
#define SKYLENS_CONVERION_H

#include "Filter.h"
#include "SED.h"
#include "Telescope.h"

namespace skylens {
  /// Flux conversion class.
  /// The formulae uastro::sed are taken from Grazian et al. (2005).
  class Conversion {
  public:
    /// Compute flux \f$f\ [\SI{}{ergs/(s~cm^2~Hz)}]\f$
    /// from AB magnitudes: \f$f=10^{-0.4(mag + 48.6)}\f$.
    static double mag2flux(double mag);
    /// Invert mag2flux().
    static double flux2mag(double flux);
    /// Compute number of photons from a source radiating with \p flux \f$f\f$
    /// when integrated for \p time \f$t\ [\SI{}{s}]\f$ with a telescope \p tel:
    /// \f[n=\frac{\pi D^2 t f}{4h}\int{\frac{T(\lambda)}{\lambda}d\lambda},\f]
    /// where \f$T(\lambda)\f$ and \f$D\ [\SI{}{cm}]\f$ describes the total 
    /// transmission and diameter of the \p tel 
    /// and \f$h = \SI{6.6260}{ergs~s}\f$ is Planck's constant.
    static double flux2photons(double flux, double time, const Telescope& tel, const Filter& transmission);
    /// Invert flux2photons.
    static double photons2flux(double photons, double time, const Telescope& tel, const Filter& transmission);
    /// Compute number of photons from a source characterized by \p emission \f$\sigma\f$:
    /// \f[n=\frac{\pi D^2 t}{4h}\int{\frac{T(\lambda)\sigma(\lambda)}{\lambda}d\lambda}\f]
    /// See flux2photons() for details.
    static double emission2photons(const SED& emission, double time, const Telescope& tel, const Filter& transmission);
    /// Convert \p photons \f$n\f$ to ADUs via \p gain \f$g\f$:
    /// \f$\text{ADU} = \frac{n}{g}\f$.
    static double photons2ADU(double photons, double gain);
    /// Invert photons2ADU().
    static double ADU2photons(double ADU, double gain);
    /// Compute flux \f$f\ [\SI{}{ergs/(s~cm^2~Hz)}]\f$ from \p ADU and
    /// the magnitude \p zeropoint \f$m_0\f$ of the telescope:
    /// \f[f = \text{ADU}\cdot 10^{-0.4(m_0 + 48.6)}\f]
    static double ADU2flux(double ADU, double zeropoint);
    /// Invert ADU2flux().
    static double flux2ADU(double flux, double zeropoint);
    /// Compute zeropoint \f$m_0\f$ of the telescope \p tel for a 
    /// exposure of \p time seconds:
    /// \f[2.5\log \Bigl(\frac{4.3035 (D/100)^2 t}{g} \int{\frac{T(\lambda)}{\lambda}d\lambda} \Bigr) + 25\f]
    static double zeroPoint(const Telescope& tel, const Filter& transmission, double time);
  };
} // end namespace
#endif
