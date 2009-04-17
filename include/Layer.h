#ifndef LAYER_H
#define LAYER_H

#include <string>
#include <map>
#include <Singleton.h>
#include <frame/Image.h>
#include <PSF.h>
#include <complex>
#include <RTree.h>
#include <list>

namespace skylens {

  /// Abstract base class for all Layer types
  class Layer {
  public:
    /// Virtual destructor.
    virtual ~Layer();
    /// Get flux at position <tt>(x,y)</tt> from this Layer.
    virtual double getFlux(double x, double y) const = 0;
    /// Get type of the Layer.
    /// The first character contains the layer type, the second the subtype
    /// - \p T: TransformationLayer
    ///   - \p TL: LensingLayer
    ///   - \p TS: ShearLayer
    ///   - \p TD: DitherLayer
    ///   - \p TM: MaskLayer
    /// - \p S: SourceLayer
    ///   - \p SG: GalaxyLayer
    ///   - \p SC: ClusterMemberLayer
    ///   - \p S*: StarLayer
    ///   - \p S=: ConstFluxLayer
    virtual std::string getType() const = 0;
    /// Get redshift of this Layer.
    double getRedshift() const;
  protected:
    /// Redshift.
    double z;
  };

  /// Stack of all Layers (ordered by redshift) in simulation.
  typedef std::map<double,Layer*> LayerStack;
  /// Type for ensuring a single LayerStack in any simulation.
  typedef skydb::Singleton< LayerStack > SingleLayerStack;

  /// LensingLayer class.
  class LensingLayer : public Layer {
  public:
    /// Constructor.
    LensingLayer(double z, std::string deflection_file);
    /// Get flux at position <tt>(x,y)</tt> from this Layer.
    virtual double getFlux(double x, double y) const;
    /// Get type of the Layer.
    /// Returns \p TL
    virtual std::string getType() const;
  private:
    shapelens::Image<complex<double> > a;
    LayerStack& ls;
    LayerStack::iterator me;
  };
 
  /// ShearLayer class.
  /// Implements transformation of a constant shear.
  class ShearLayer : public Layer {
  public:
    /// Constructor.
    ShearLayer(double z, complex<double> gamma);
    /// Get flux at position <tt>(x,y)</tt> from this Layer.
    virtual double getFlux(double x, double y) const;
    /// Get type of the Layer.
    /// Returns \p TS
    virtual std::string getType() const;
  private:
    complex<double> gamma;
    LayerStack& ls;
    LayerStack::iterator me;
  };

  /// DitherLayer class.
  class DitherLayer : public Layer {
  public:
    /// Constructor.
    DitherLayer(double z, double dx, double dy);
    /// Get flux at position <tt>(x,y)</tt> from this Layer.
    virtual double getFlux(double x, double y) const;
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
  /// The coordinates of the mask polygons are given in units of the total
  /// FoV, i.e. <tt>(0,0)</tt> defines the bottom-left corner, <tt>(1,1)</tt>
  /// the position outside the top-right corner.
  class MaskLayer : public Layer {
  public:
    /// Constructor.
    /// \p FoV is given in \p arcsec and polygon coordinates in units of \p FoV.
    MaskLayer(double z, double FoV, const std::list<shapelens::Polygon<double> >& masks);
    /// Constructor from a mask file.
    MaskLayer(double z, double FoV, std::string maskfile);
    /// Get flux at position <tt>(x,y)</tt> from this Layer.
    virtual double getFlux(double x, double y) const;
    /// Get type of the Layer.
    /// Returns \p TM
    virtual std::string getType() const;
    /// Clear set of masked areas.
    void clearMasks();
    /// Add masked area.
    void addMask(const SpatialIndex::IShape& shape);
  private:
    std::list<shapelens::Polygon<double> > masks;
    LayerStack& ls;
    LayerStack::iterator me;
  };

  /// GalaxyLayer class.
  class GalaxyLayer : public Layer {
  public:
    /// Constructor.
    // FIXME: what arguments for constructor???
    GalaxyLayer(double z);
    /// Get flux at position <tt>(x,y)</tt> from this Layer.
    virtual double getFlux(double x, double y) const;
    /// Get type of the Layer.
    /// Returns \p SG
    virtual std::string getType() const;
  private:
    LayerStack& ls;
  };

  /// ClusterMemberLayer class.
  class ClusterMemberLayer : public Layer {
  public:
    /// Constructor.
    // FIXME: what arguments for constructor???
    ClusterMemberLayer(double z);
    /// Get flux at position <tt>(x,y)</tt> from this Layer.
    virtual double getFlux(double x, double y) const;
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
    /// FIXME: what arguments for constructor???
    StarLayer(double z, const PSF& psf);
    /// Get flux at position <tt>(x,y)</tt> from this Layer.
    virtual double getFlux(double x, double y) const;
    /// Get type of the Layer.
    /// Returns \p S*
    virtual std::string getType() const;
  private:
    LayerStack& ls;
    const PSF& psf;
  };

  /// ConstFluxLayer class.
  /// This layer is mainly thought for testing purpose, as it creates a
  /// constant sheat of light.
  class ConstFluxLayer : public Layer {
  public:
    /// Constructor.
    ConstFluxLayer(double z, double flux);
    /// Get flux at position <tt>(x,y)</tt> from this Layer.
    virtual double getFlux(double x, double y) const;
    /// Get type of the Layer.
    /// Returns \p S=
    virtual std::string getType() const;
  private:
    LayerStack& ls;
    double flux;
  };

  

} // end namespace

#endif
