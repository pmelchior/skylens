#ifndef SKYLENS_LAYER_H
#define SKYLENS_LAYER_H

#include <string>
#include <map>
#include <shapelens/Image.h>
#include <complex>
#include <list>
#include "PSF.h"
#include "LensingInformation.h"
#include "Singleton.h"
#include "SourceModel.h"
#include "Cosmology.h"
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>

namespace skylens {

  /// The global cosmological model.
  /// As long as it is not changed, it's a vanilla LCMD model.
  typedef skylens::Singleton<Cosmology> SingleCosmology;

  /// Abstract base class for all Layer types
  class Layer {
  public:
    /// Virtual destructor.
    virtual ~Layer(); 
    /// Get flux at position \p P from the Layer.
    /// If \p z is set, only the source layer at the specified redshift 
    /// will contribute flux, while transformation layers act normally.
    virtual double getFlux(const shapelens::Point<double>& P, double* z=NULL) const = 0;
    /// Get type of the Layer.
    /// The first character contains the layer type, the second the subtype
    /// - \p T: TransformationLayer
    ///   - \p TL: LensingLayer
    ///   - \p TS: ShearLayer
    ///   - \p TF: FlexionLayer
    ///   - \p TD: DitherLayer
    ///   - \p TM: MaskLayer
    ///   - \p TC: ConvolutionLayer
    ///   - \p T0: NullLayer
    /// - \p S: SourceLayer
    ///   - \p SG: GalaxyLayer
    ///   - \p SC: ClusterMemberLayer
    ///   - \p S*: StarLayer
    ///   - \p SS: SkyFluxLayer
    virtual std::string getType() const = 0;
    /// Get redshift of this Layer.
    double getRedshift() const;
    /// The transparency of this Layer.
    /// If <tt>transparent == true</tt>, this Layer behaves as if
    /// if is not present.
    bool transparent;
  protected:
    /// Redshift.
    double z;
  };

  /// Stack of all Layers (ordered by redshift) in simulation.
  class LayerStack : public std::multimap<double,Layer*> {
  public:
    /// Constructor.
    LayerStack();
    /// Destructor.
    /// Deletes all Object with pointers in LayerStack.
    ~LayerStack();
    /// Write to any \p std::ostream
    friend std::ostream& operator<<(std::ostream& os, const LayerStack& ls);
  };
  /// Type for ensuring a single LayerStack in any simulation.
  typedef skylens::Singleton< LayerStack > SingleLayerStack;

  /// LensingLayer class.
  class LensingLayer : public Layer {
  public:
    /// Constructor from a deflection angle FITS file.
    /// The units of angles are radians, compute for the lens at redshift
    /// \p ZLENS and sources at redshift \p ZSOURCE; \p SIDEL is the 
    /// horizontal size of the lens plane in Mpc/h.\n
    /// \b CAUTION: These parameters need to be set in the FITS header.
    LensingLayer(double z, std::string deflection_file, const shapelens::Point<double>* center = NULL);
    /// Get flux at position \p P from the Layer.
    /// If \p z is set, only the source layer at the specified redshift 
    /// will contribute flux, while transformation layers act normally.
    virtual double getFlux(const shapelens::Point<double>& P, double* z=NULL) const;
    /// Get type of the Layer.
    /// Returns \p TL
    virtual std::string getType() const;
    /// Get the position of the center of the lens.
    shapelens::Point<double> getCenter() const;
    /// Get set of critical points for a source layer at redshift \p zs.
    /// If \p det_sign is non-zero, method returns tangential (> 0) or radial 
    /// (< 0) critical points only. 
    std::vector<shapelens::Point<double> > findCriticalPoints(double zs, int det_sign=0) const;
    /// Get shear at lens-plane position \f$\theta\f$ for source redshift 
    /// \f$z_s\f$.
    std::complex<double> getShear(const shapelens::Point<double>& theta, double zs, bool reduced) const;
    /// Get convergence at lens-plane position \f$\theta\f$ for source 
    /// redshift \f$z_s\f$.
    double getConvergence(const shapelens::Point<double>& theta, double zs) const;
    /// Set second derivatives of potential \f$\phi\f$ at \f$\theta\f$
    /// for source redshift \f$z_s\f$.
    void set_Dphi(const Point<double>& theta, double zs, double& phixx, double& phixy, double& phiyx, double& phiyy) const;
    /// Get source-plane position \f$\beta\f$ for given image position 
    /// \f$\theta\f$ and source redshift \f$z_s\f$.
    Point<double> getBeta(const Point<double>& theta, double zs) const;
    /// Find (multiple) image positions for source-plane position \f$\beta\f$
    /// and redshift \f$z_s\f$.
    std::vector<Point<double> > findImages(const Point<double>& beta, double zs, const Rectangle<double>& area) const;

  private:
    shapelens::Image<std::complex<float> > a;
    float scale0, theta0;
    LayerStack& ls;
    LayerStack::iterator me;
    LensingInformation& li;
    Cosmology& cosmo;
    void setDistances(const LayerStack::iterator& iter) const;
    std::list<Rectangle<double> > getCellsEnclosing(const Point<double>& beta, double zs, const Rectangle<double>& area, int C) const;
    void finiteDifferences(const Point<int>& P0, double& phixx, double& phixy, double& phiyx, double& phiyy) const;
  };
 
  /// ShearLayer class.
  /// Implements transformation of a constant shear.
  class ShearLayer : public Layer {
  public:
    /// Constructor.
    ShearLayer(double z, std::complex<double> gamma);
    /// Get flux at position \p P from the Layer.
    /// If \p z is set, only the source layer at the specified redshift 
    /// will contribute flux, while transformation layers act normally.
    virtual double getFlux(const shapelens::Point<double>& P, double* z=NULL) const;
    /// Get type of the Layer.
    /// Returns \p TS
    virtual std::string getType() const;
  private:
    std::complex<double> gamma;
    LayerStack& ls;
    LayerStack::iterator me;
  };

  /// NullLayer class.
  /// The Layer is made for ensuring that all Layers are connected,
  /// no matter what configuration the LayerStack has.
  /// It is set to <tt>z=-1000</tt> to always be at the the beginning of
  /// the stack and acts as an empty TransformationLayer.
  class NullLayer : public Layer {
  public:
    /// Constructor.
    /// The NullLayer will be inserted in the LayerStack at redshift 
    /// <tt>z = -1000</tt>.
    NullLayer();
    /// Get flux at position \p P from the Layer.
    /// If \p z is set, only the source layer at the specified redshift 
    /// will contribute flux, while transformation layers act normally.
    virtual double getFlux(const shapelens::Point<double>& P, double* z=NULL) const;
    /// Get type of the Layer.
    /// Returns \p T0
    virtual std::string getType() const;
  private:
    LayerStack& ls;
    LayerStack::iterator me;
  };

  /// DitherLayer class.
  class DitherLayer : public Layer {
  public:
    /// Constructor.
    /// The DitherLayer will be inserted in the LayerStack at redshift 
    /// <tt>z = -3</tt>.
    DitherLayer(double dx, double dy);
    /// Get flux at position \p P from the Layer.
    /// If \p z is set, only the source layer at the specified redshift 
    /// will contribute flux, while transformation layers act normally.
    virtual double getFlux(const shapelens::Point<double>& P, double* z=NULL) const;
    /// Get type of the Layer.
    /// Returns \p TD
    virtual std::string getType() const;
    /// Set new displacement values.
    void setDisplacement(double dx, double dy);
  private:
    double dx, dy;
    LayerStack& ls;
    LayerStack::iterator me;
  };

  /// MaskLayer class.
  /// MaskLayer is used to mask out regions of the final image by defining
  /// a set of simple polygons (polygons with no overlapping or crossing edges).\n
  /// The coordinates of the mask polygons are given in units of \p arcsec
  /// such that <tt>(0,0)</tt> defines the bottom-left corner and 
  /// <tt>(FOV_X, FOV_Y)</tt> the top-right corner, where \p FOV_X/FOV_Y
  /// are given in the telescope configuration file.
  class MaskLayer : public Layer {
  public:
    /// Constructor.
    /// Polygon coordinates have to be given in units of \p arcsec.
    /// The MaskLayer will be inserted in the LayerStack at redshift 
    /// <tt>z = -2</tt>.
    MaskLayer(const std::list<shapelens::Polygon<double> >& masks);
    /// Constructor from a mask file.
    MaskLayer(std::string maskfile);
    /// Get flux at position \p P from the Layer.
    /// If \p z is set, only the source layer at the specified redshift 
    /// will contribute flux, while transformation layers act normally.
    virtual double getFlux(const shapelens::Point<double>& P, double* z=NULL) const;
    /// Get type of the Layer.
    /// Returns \p TM
    virtual std::string getType() const;
    /// Clear set of masked areas.
    void clearMasks();
  private:
    std::list<shapelens::Polygon<double> > masks;
    LayerStack& ls;
    LayerStack::iterator me;
  };

  /// ConvolutionLayer class.
  /// The ConvolutionLayer contains an (super-)image of the unconvolved
  /// Layers at higher redhift and convolves it with the PSF.\n
  /// Sampling with getFlux() uses pixel interpolation.
  class ConvolutionLayer : public Layer {
  public:
    /// Constructor.
    /// \p FoV is given in \p arcsec, \p pixsize in <tt>arsec/pixel</tt>.
    /// The ConvolutionLayer will be inserted in the LayerStack at redshift 
    /// <tt>z = 0</tt>.
    ConvolutionLayer(double FOV, double pixsize, const PSF& psf);
    /// Get flux at position \p P from the Layer.
    /// If \p z is set, only the source layer at the specified redshift 
    /// will contribute flux, while transformation layers act normally.
    virtual double getFlux(const shapelens::Point<double>& P, double* z=NULL) const;
    /// Get type of the Layer.
    /// Returns \p TC.
    virtual std::string getType() const;
    /// Clear the superimage.
    void clear();
    shapelens::Image<double> im;
  private:
    int L, PAD;
    double pixsize;
    const PSF& psf;
    LayerStack::iterator me;
    LayerStack& ls;
    void convolveImage(shapelens::Image<double>& im, shapelens::Object& kernel);
  };

  
  typedef boost::geometry::model::point<float, 2, boost::geometry::cs::cartesian> BPoint;
  typedef boost::geometry::model::box<BPoint> BBox;
  typedef std::pair<BBox, size_t> BBoxIndex;
  typedef boost::geometry::index::rtree<BBoxIndex, boost::geometry::index::rstar<16, 4> > RTree;

  /// GalaxyLayer class.
  class GalaxyLayer : public Layer {
  public:
    /// Constructor.
    GalaxyLayer(double z, const SourceModelList& galaxies);
    /// Get flux at position \p P from the Layer.
    /// If \p z is set, only the source layer at the specified redshift 
    /// will contribute flux, while transformation layers act normally.
    virtual double getFlux(const shapelens::Point<double>& P, double* z=NULL) const;
    /// Get type of the Layer.
    /// Returns \p SG
    virtual std::string getType() const;
  private:
    LayerStack& ls;
    SourceModelList galaxies;
    RTree rtree;
  };

  /// ClusterMemberLayer class.
  class ClusterMemberLayer : public Layer {
  public:
    /// Constructor.
    // FIXME: what arguments for constructor???
    ClusterMemberLayer(double z);
    /// Get flux at position \p P from the Layer.
    /// If \p z is set, only the source layer at the specified redshift 
    /// will contribute flux, while transformation layers act normally.
    virtual double getFlux(const shapelens::Point<double>& P, double* z=NULL) const;
    /// Get type of the Layer.
    /// Returns \p SC
    virtual std::string getType() const;
  private:
    LayerStack& ls;
  };

  /// StarLayer class.
  class StarLayer : public Layer {
  public:
    /// Constructor.
    /// The StarLayer will be inserted in the LayerStack at redshift 
    /// <tt>z = 0</tt>.
    StarLayer(const SourceModelList& stars);
    /// Get flux at position \p P from the Layer.
    /// If \p z is set, only the source layer at the specified redshift 
    /// will contribute flux, while transformation layers act normally.
    virtual double getFlux(const shapelens::Point<double>& P, double* z=NULL) const;
    /// Get type of the Layer.
    /// Returns \p S*
    virtual std::string getType() const;
  private:
    LayerStack& ls;
    const SourceModelList& stars;
    RTree rtree;
  };

  /// SkyFluxLayer class.
  /// This layer creates a constant sheat of light to account for
  /// sky background brightness.
  class SkyFluxLayer : public Layer {
  public:
    /// Constructor.
    /// The SkyFluxLayer will be inserted in the LayerStack at redshift 
    /// <tt>z = -1</tt>.
    SkyFluxLayer(double flux);
    /// Get flux at position \p P from the Layer.
    /// If \p z is set, only the source layer at the specified redshift 
    /// will contribute flux, while transformation layers act normally.
    virtual double getFlux(const shapelens::Point<double>& P, double* z=NULL) const;
    /// Get type of the Layer.
    /// Returns \p SS
    virtual std::string getType() const;
    /// Set flux.
    void setFlux(double flux);
  private:
    LayerStack& ls;
    double flux;
  };
  

} // end namespace

#endif
