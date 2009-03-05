#ifndef TELESCOPE_H
#define TELESCOPE_H

#include <filter.h>
#include <string>

/// Namespace for SkyLens++
namespace skylens {

  /// Abstract base class for all telescopes.
  class Telescope {
  public:
    // no constructor, as we want to use this as abstract base class
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
    std::string filter_name;
    /// Read configuration file \p telescope.conf from \p path.
    void readConfig(std::string path);
  };


  // derived class from here on...

  class LBT : public Telescope {
  public:
    LBT(std::string band, std::string instrument);
  };
  class HST : public Telescope {
  public:
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
