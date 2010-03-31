#ifndef SKYLENS_SOURCECATALOG_H
#define SKYLENS_SOURCECATALOG_H

#include <map>
#include <set>
#include <string>
#include "../include/Telescope.h"
#include "../include/Helpers.h"
#include <shapelens/frame/Point.h>
#include <libastro/filter.h>
#include <libastro/sed.h>


namespace skylens {
  struct GalaxyInfo {
    std::map<std::string, std::pair<double,double> > mags;
    double redshift;
    std::string sed;
    double radius;
    double ellipticity;
    double n_sersic;
    unsigned int model_type;
    unsigned long object_id;
    double mag;
    /// map: band name -> ADU of source per second
    std::map<std::string, double> adus;
    shapelens::Point<double> centroid;
    double redshift_layer;
    double rotation;
  };
  
  /// Source catalogs for use in SkyLens.
  class SourceCatalog : public std::map<unsigned long, GalaxyInfo> {
  public:
    /// Default constructor.
    SourceCatalog() {}
    /// Constructor.
    /// \p configfile must obey the standards of a SourceCatalog configuration file.
    SourceCatalog(std::string configfile);
    /// Constructor from save catalog file.
    /// \p configfile must obey the standards of a SourceCatalog configuration file,
    /// \p catalogfile must by written obeying the format used by save().
    SourceCatalog(std::string configfile, std::string catalogfile);
    /// Adjust galaxy numbers to account for FoV change from reference catalog
    /// to \p tel.
    void adjustNumber(const Telescope& tel);
    /// Distribute sources randomly in FoV and assign it to closest
    /// GalaxyLayer.
    void distribute(const Telescope& tel);
    /// Choose the bands from reference catalog which have at least \p fraction
    /// overlap with \p tel.total .
    void selectOverlapBands(const Telescope& tel, double fraction=0.1);
    /// Compute ADU per second for each source in each band indentified by
    /// selectOverlapBands().
    void computeADUinBands(const Telescope& tel);
    /// Create GalaxyLayer instances of all sources, based on the redshift list given
    /// in configuration file at construction time.
    void createGalaxyLayers() const;
    /// Save catalog in specific format to \p filename.
    void save(std::string filename) const;

    /// Container for filter curves and database details for each band
    class Band {
    public:
      /// Name of filter.
      std::string name;
      /// Filter curve.
      filter curve;
      /// Overlap with observation filter.
      double overlap;
      /// map: model_type -> DB details
      std::map<char, std::string> dbdetails;
      bool operator<(const Band& b) const;
    };
    /// Details of reference catalog.
    class ImagingReference {
    public:
      /// Pixel size.
      double pixsize;
      /// Field-of-view.
      double fov;
      /// set of Band, sorted by central wavelength.
      std::set<Band> bands;
      /// map: source SED number -> sed
      std::map<std::string, sed> seds;
    };
  private:
    ImagingReference imref;
    std::string tablename, query;
    std::set<double> redshifts;
    double getRedshiftNearestLayer(double z);
    void computeADU(GalaxyInfo& info, const Telescope& tel);
    void readConfig(std::string configfile);
  };
} // end namespace

#endif
