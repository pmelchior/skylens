#ifndef SKYLENS_OBSERVATION_H
#define SKYLENS_OBSERVATION_H

#include "Layer.h"
#include "Telescope.h"
#include <sed.h>
#include <shapelens/frame/Image.h>
#include <gsl/gsl_rng.h>

namespace skylens {
  /// Observation class.
  /// This class allows to raytrace through the existing LayerStack and
  /// thus to perform the actual observation.\n\n
  /// The constructors also create the SkyFluxLayer and set
  /// the noise characterisitic according to eq. (31) in Meneghetti 
  /// et al. (2008):\n
  /// \p ron is the term \f$n_{exp}\bigl(RON/GAIN/bigr)\f$,
  /// \p flat_field the term \f$(f+ a^2/n_{exp}^2)\f$  
  class Observation {
  public:
    /// Default constructor.
    Observation(const Telescope& tel, double exptime);
    /// Draw the observation onto \p im.
    /// Shoots rays from the center of each pixel of \p im through
    /// all layers of the LayerStack.\n\n
    /// If \p auto_adjust is set to \p true, \p im is resized as to 
    /// describe the entire FoV as seen be the Telescope \p tel; otherwise,
    /// the pixel positions are taken from the Grid of \p im, such that one 
    /// can specify the particular field of the observation. This requires
    /// a proper setup of WCS of \p im.
    void makeImage(shapelens::Image<double>& im, bool auto_adjust = true) const;
    /// Get total transmission \f$T(\lambda)\f$, including atmospheric
    /// extinction.
    /// According to // according to Grazian et al. (2004), eq. 3
    const filter& getTotalTransmittance() const;

    void createSkyFluxLayer(const sed& sky);
    void createSkyFluxLayer(double sky_mag);
    void computeTransmittance(const filter& atmosphere, double airmass);
    void computeTransmittance(double absorption, double airmass);
    void setNoise(int nexp=1);
  private:
    const Telescope& tel;
    bool hasNoise;
    filter total_air;
    double time, ron, flat_field;
    void addNoise(const gsl_rng* r, double& signal) const;
    
  };
} // end namespace

#endif
