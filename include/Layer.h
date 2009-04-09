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
  class MaskLayer : public Layer {
  public:
    /// Constructor.
    MaskLayer(double z, const std::list<SpatialIndex::IShape>& ls);
    /// Get flux at position <tt>(x,y)</tt> from this Layer.
    virtual double getFlux(double x, double y) const;
    /// Get type of the Layer.
    /// Returns \p TM
    virtual std::string getType() const;
  private:
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

} // end namespace

#endif
