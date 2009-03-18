#ifndef PSF_H
#define PSF_H

#include <string>
#include <frame/Object.h>

namespace skylens {
  class PSF {
  public:
    /// Default constructor.
    PSF();
    /// Constructor from a file.
    /// \p filename can be either a FITS file or a SIF file.
    PSF(std::string filename);
    /// Get pixelized shape of PSF at given position.
    const shapelens::Object& getShape() const;
  private:
    shapelens::Object psf;
  };
}

#endif
