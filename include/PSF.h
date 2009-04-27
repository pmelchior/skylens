#ifndef PSF_H
#define PSF_H

#include <string>
#include <frame/Object.h>
#include <shapelets/ShapeletObject.h>

namespace skylens {
  class PSF {
  public:
    /// Default constructor.
    PSF();
    /// Constructor from a SIF file.
    PSF(std::string siffile);
    /// Get shapelet model of PSF
    const shapelens::ShapeletObject& getShape() const;
  private:
    shapelens::ShapeletObject psf;
  };
}

#endif
