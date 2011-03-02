#ifndef SKYLENS_TELESCOPE_H
#define SKYLENS_TELESCOPE_H

#include "PSF.h"
#include <libastro/filter.h>
#include <string>

/// Namespace for SkyLens++
namespace skylens {

  /// Class for telescopes.
  class Telescope {
  public:
    /// Constructor.
    Telescope();
    /// Argumented constructor.
    /// \p configfile denotes the name of the telescope's config file
    /// and \p band denotes the filter FITS file \f$F(\lambda)\f$
    /// to be used during observation.\n\n
    Telescope(std::string configfile, std::string bandfile);
    /// total filter shape, ignoring airmass extinction.
    /// \f[T(\lambda) = C(\lambda)\,F(\lambda)\,M(\lambda)\,O(\lambda)\f]
    astro::filter total;
    /// mirror diameter
    double diameter;
    /// flat-field accuracy.
    double flat_acc;
    /// pixel size
    double pixsize;
    /// Field of view in x direction
    double fov_x;
    /// Field of view in y direction
    double fov_y;
    /// gain
    double gain;
    /// read-out noise
    double ron;
    /// nickname
    std::string name;
    /// name of filter band
    std::string bandname;
    /// PSF shape.
    //PSF psf;
  private:
    /// Read configuration file \p telescope.conf from \p path.
    void readConfig(std::string path, std::string configfile);
  };
} // end namespace
#endif
