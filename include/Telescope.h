#ifndef TELESCOPE_H
#define TELESCOPE_H

#include <PSF.h>
#include <filter.h>
#include <string>

/// Namespace for SkyLens++
namespace skylens {

  /// Abstract base class for all telescopes.
  class Telescope {
  public:
    /// Get mirror diameter.
    double getDiameter() const;
    /// Get pixel scale of CCD.
    double getPixelScale() const;
    /// Get gain of CCD.
    double getGain() const;
    /// Get read-out noise of CCD.
    double getReadOutNoise() const;
    /// Get `nickname` of the telescope
    std::string getName() const;
    /// Get the name of the filter band.
    std::string getFilterName() const;
    /// Get filter shape.
    const filter& getFilter() const;
    /// Get PSF shape;
    const PSF& getPSF() const;

  protected:
    /// total filter shape.
    filter total;
    /// mirror diameter
    double d;
    /// flat-field accuracy.
    double flat_acc;
    /// pixel size
    double px;
    /// gain
    double gain;
    /// read-out noise
    double ron;
    /// nickname
    std::string name;
    /// name of filter band
    std::string band;
    /// PSF shape.
    PSF psf;
    /// Read configuration file \p telescope.conf from \p path.
    void readConfig(std::string path);
  };


  // derived class from here on...

  class LBT : public Telescope {
  public:
    LBT(std::string band, std::string instrument);
  };

  /// HST telescope.
  class HST : public Telescope {
  public:
    /// Constructor for given \p band and \p instrument.
    /// Available instruments: 
    /// - ACS_WFC
    /// - NICMOS_NIC1
    /// - NICMOS_NIC2
    /// - NICMOS_NIC3
    /// - WFC3_IR
    /// - WFC3_UVIS
    HST(std::string band, std::string instrument);
  };
  class VLT : public Telescope {
  public:
    VLT(std::string band, std::string instrument);
  };
  class Euclid : public Telescope {
  public:
    Euclid(std::string band);
  };
  class WFI : public Telescope {
  public:
    WFI(std::string band);
  };
  /// SUBARU telescope.
  class SUBARU : public Telescope {
  public:
    /// Constructor for given \p band.
    SUBARU(std::string band);
  };
  class CFHT : public Telescope {
  public:
    CFHT(std::string band);
  };
}

#endif
