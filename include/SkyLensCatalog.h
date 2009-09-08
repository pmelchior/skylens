#ifndef SKYLENS_SKYLENSCATALOG_H
#define SKYLENS_SKYLENSCATALOG_H

#include <map>
#include <string>
#include <shapelens/frame/Shapes.h>
#include <shapelens/frame/Catalog.h>


namespace skylens {
  struct GalaxyInfo {
    std::map<std::string, double> mags;
    double redshift;
    double sed;
    double radius;
    double ellipticity;
    double n_sersic;
    unsigned int model_type;
    unsigned long object_id;
    double flux;
    shapelens::Point<double> centroid;
    shapelens::Rectangle<double> bb;
  };
  
  /// Galaxy catalogs for use in SkyLens.
  class SkyLensCatalog : public std::map<unsigned long, GalaxyInfo> {
  public:
    /// Constructor
    /// \p where is the content of a <tt> SQL WHERE</tt> statement
    /// and each entry thereof is joined by a logical \p AND.
    SkyLensCatalog(const std::map<std::string, std::string>& where);
    /// Adjust galaxy numbers to account for FoV change between
    /// \p FoV_ref and \p FoV_out.
    void adjustGalaxyNumber(double FoV_ref, double FoV_out);
    const std::string dbname;
    const std::string tablename;
    const std::map<std::string, std::string>& where;
    /// Save catalog in SExtractor-like ASCII format.
    shapelens::Catalog getCatalog();
  };
} // end namespace

#endif
