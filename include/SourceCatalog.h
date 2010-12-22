#ifndef SKYLENS_SOURCECATALOG_H
#define SKYLENS_SOURCECATALOG_H

#include <list>
#include <map>
#include <set>
#include <string>
#include "../include/Telescope.h"
#include "../include/Helpers.h"
#include <shapelens/frame/Point.h>
#include <shapelens/modelfit/SourceModel.h>
#include <shapelens/utils/Property.h>
#include <shapelens/utils/SQLiteDB.h>
#include <libastro/filter.h>
#include <libastro/sed.h>


namespace skylens {
  /// Container for galactic information.
  struct GalaxyInfo {
    /// Object id in reference catalog
    unsigned long object_id;
    /// map: band name -> (magnitude/error)
    std::map<std::string, std::pair<double,double> > mags;
    /// redshift
    double redshift;
    /// SED name/identifier
    std::string sed;
    /// Model type
    unsigned int model_type;
    /// SED normalization.
    double sed_norm;
    /// Magnitude in simulation
    double mag;
    /// map: band name -> simulated ADU per second
    std::map<std::string, double> adus;
    /// Position in simulation
    shapelens::Point<double> centroid;
    /// Nearest GalaxyLayer in simulation
    double redshift_layer;
    /// Rotation w.r.t. to observed model (negative value indicates reflection)
    double rotation;
  };
  
  /// Source catalogs for use in SkyLens.
  /// This class loads information about source galaxies from a reference 
  /// catalog in the DB. It can modify these informations (numbers in FoV,
  /// flux in ADU, positions) to accord the simulation.\n\n
  /// The catalog entries can be written to an ASCII file, or
  /// reloaded from such a file.\n\n
  /// The class furthermore has/needs access to the DB tables which contain
  /// the model information. With this, it can create the GalaxyLayer instances.
  class SourceCatalog : public std::list<GalaxyInfo> {
  public:
    /// Default constructor.
    SourceCatalog();
    /// Constructor.
    /// \p configfile must obey the standards of a SourceCatalog configuration file.
    SourceCatalog(std::string configfile);
    /// Constructor from saved catalog data.
    /// \p i denotes a running number to discriminate different
    /// catalogs in one database.
    SourceCatalog(shapelens::SQLiteDB& db, int i=0);
    /// Adjust galaxy numbers to account for FoV change from reference catalog
    /// to \p fov.
    void adjustNumber(const shapelens::Point<double>& fov);
    /// Get replication ratio.
    /// Defined as \f$r\equiv N_{sim}/N_{ref}\f$, the ratio between the number
    /// of simulated galaxies and the number of galaxies in the reference
    /// catalog. \f$r>1\f$ indicates source replication.
    double getReplicationRatio() const;
    /// Distribute sources randomly in FoV and assign it to closest
    /// GalaxyLayer.
    /// If \p keepPosition is true, only it leaves the centroid unchanged.
    void distribute(const shapelens::Point<double>& fov, bool keepPosition = false);
    /// Compute ADU per second for each source in each band indentified by
    /// selectOverlapBands().
    void computeADUinBands(const Telescope& tel, const filter& transmittance);
    /// Create GalaxyLayer instances of all sources, based on the redshift list given
    /// in configuration file at construction time.
    /// \p exptime is the exposure time in seconds.
    void createGalaxyLayers(double exptime);
    /// Save catalog to SQLite database.
    /// \p i denotes a running number to discriminate different
    /// catalogs in one database.
    void save(shapelens::SQLiteDB& db, int i=0) const;

    /// Container for filter curves and database details for each band
    class Band {
    public:
      /// Name of filter.
      std::string name;
      /// Filter curve.
      filter curve;
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
    ImagingReference imref;
    shapelens::Property config;
  private:
    std::string tablename, query, where;
    std::map<double, shapelens::SourceModelList> layers;
    double getRedshiftNearestLayer(double z);
    void parseConfig(std::string configfile = "");
    void setRotationMatrix(NumMatrix<double>& O, double rotation) const;
    std::set<char> need_model;
    std::map<char, std::map<std::string, std::string> > model_band_tables;
    double replication_ratio;
  };
} // end namespace

#endif
