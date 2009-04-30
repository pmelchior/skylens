#ifndef SKYLENS_OBSERVATION_H
#define SKYLENS_OBSERVATION_H

#include <Layer.h>
#include <Telescope.h>
#include <filter.h>
#include <frame/Image.h>

namespace skylens {
  /// Observation class.
  /// This class allows to raytrace through the existing LayerStack and
  /// thus to perform the actual observation.
  class Observation {
  public:
    /// Constructor.
    /// \p sky contains a spectrum of the sky from which the sky brightness
    /// in the observational band is determined, \p time is the integration
    /// time, and \p n_exposures is the number of identical exposures of the
    /// same field, which affects the noise variance.
    Observation(const Telescope& t, const filter& sky, double time, int n_exposures = 1);
    /// Draw the observation onto \p im.
    /// Shoots rays from the center of each pixel of \p im trhough
    /// all layers of the LayerStack.\n\n
    /// If \p auto_adjust is set to \p true, \p im is resized as to 
    /// describe the entire FoV as seen be the Telescope \p tel; otherwise,
    /// the pixel positions are taken from the Grid of \p im, such that one 
    /// can specify the particular field of the observation.
    void makeImage(shapelens::Image<double>& im, bool auto_adjust = true);
  private:
    const Telescope& tel;
    const filter& sky;
    double time;
    int nexp;
    double computeSkyFlux();
  };
} // end namespace

#endif
