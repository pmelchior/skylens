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
#include <libastro/filter.h>
#include <libastro/sed.h>


namespace skylens {
  /// Container for galactic information.
  struct GalaxyInfo {
    /// map: band name -> (magnitude/error)
    std::map<std::string, std::pair<double,double> > mags;
    /// redshift
    double redshift;
    /// SED name/identifier
    std::string sed;
    /// Effective radius
    double radius;
    /// Ellipticity
    double ellipticity;
    /// Sersic index
    double n_sersic;
    /// Model type
    unsigned int model_type;
    /// Object id in reference catalog
    unsigned long object_id;
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
  
  typedef std::map<unsigned long, boost::shared_ptr<shapelens::ShapeletObject> > ShapeletObjectCat;

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
    /// Constructor from save catalog file.
    /// \p configfile must obey the standards of a SourceCatalog configuration file,
    /// \p catalogfile must by written obeying the format used by save().
    SourceCatalog(std::string configfile, std::string catalogfile);
    /// Adjust galaxy numbers to account for FoV change from reference catalog
    /// to \p tel.
    void adjustNumber(const Telescope& tel);
    /// Get replication ratio.
    /// Defined as \f$r\equiv N_{sim}/N_{ref}\f$, the ratio between the number
    /// of simulated galaxies and the number of galaxies in the reference
    /// catalog. \f$r>1\f$ indicates source replication.
    double getReplicationRatio() const;
    /// Distribute sources randomly in FoV and assign it to closest
    /// GalaxyLayer.
    void distribute(const Telescope& tel);
    /// Choose the bands from reference catalog which have at least \p fraction
    /// overlap with \p transmittance of Observation.
    void selectOverlapBands(const filter& transmittance, double fraction=0.1);
    /// Compute ADU per second for each source in each band indentified by
    /// selectOverlapBands().
    void computeADUinBands(const Telescope& tel, const filter& transmittance);
    /// Create GalaxyLayer instances of all sources, based on the redshift list given
    /// in configuration file at construction time.
    /// \p exptime is the exposure time in seconds.
    void createGalaxyLayers(double exptime);
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
    std::string tablename, query, where;
    std::map<double, shapelens::SourceModelList> layers;
    double getRedshiftNearestLayer(double z);
    void readConfig(std::string configfile);
    void setRotationMatrix(NumMatrix<double>& O, double rotation) const;
    std::map<std::string, ShapeletObjectCat> getShapeletModels();
    double replication_ratio;
  };
} // end namespace

#endif
