#ifndef SKYLENS_PSF_H
#define SKYLENS_PSF_H

#include <string>
#include <shapelens/frame/Object.h>
#include <shapelens/shapelets/ShapeletObject.h>

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
