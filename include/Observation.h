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
    /// Constructor.
    /// \p sky contains a spectrum of the sky from which the sky brightness
    /// in the observational band is determined, \p time is the integration
    /// time, and \p n_exposures is the number of identical exposures of the
    /// same field, which affects the noise variance.
    Observation(const Telescope& t, const sed& sky, double time, int n_exposures = 1);
    /// Constructor.
    /// \p sky_mag gives magintude of the sky (assuming flat spectrum) 
    /// in the observational band, \p time is the integration
    /// time, and \p n_exposures is the number of identical exposures of the
    /// same field, which affects the noise variance.
    Observation(const Telescope& t, double sky_mag, double time, int n_exposures = 1);
    /// Destructor.
    ~Observation();
    /// Draw the observation onto \p im.
    /// Shoots rays from the center of each pixel of \p im trhough
    /// all layers of the LayerStack.\n\n
    /// If \p auto_adjust is set to \p true, \p im is resized as to 
    /// describe the entire FoV as seen be the Telescope \p tel; otherwise,
    /// the pixel positions are taken from the Grid of \p im, such that one 
    /// can specify the particular field of the observation.
    void makeImage(shapelens::Image<double>& im, bool auto_adjust = true);
    /// The subpixel sampling.
    /// \p SUBPIX is the number of samplings per pixel in each direction:\n
    /// * <tt>SUBPIX=1</tt>: centered sampling (default)
    /// * <tt>SUBPIX=2</tt>: centered sampling within 2x2 subpixels
    /// * ...
    static unsigned int SUBPIX;
  private:
    const Telescope& tel;
    gsl_rng * r;
    double time, ron, flat_field;
    int nexp;
    void setNoise();
    void addNoise(double& signal);
  };
} // end namespace

#endif
