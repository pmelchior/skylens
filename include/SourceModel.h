#ifndef SKYLENS_SOURCEMODEL_H
#define SKYLENS_SOURCEMODEL_H

#include <vector>
#include <boost/shared_ptr.hpp>
#include <shapelens/Object.h>
#include <shapelens/Catalog.h>
#include <shapelens/Shapes.h>

namespace skylens {
  using namespace shapelens;
  /// Base class for idealized source models.
  /// A galaxy model is a idealized representation of a two-dimensional shape.
  /// It's main advantage: It can be sampled at any resolution.
  class SourceModel {
  public:
    /// Destructor.
    virtual ~SourceModel();
    /// Sample model at \p P.
    virtual double getValue(const Point<double>& P) const = 0;
    /// Get rectangluar support area of the model.
    const Rectangle<double>& getSupport() const;
    /// Whether model contains this Point.
    virtual bool contains(const Point<double>& P) const; 
    /// Get centroid of model.
    const Point<double>& getCentroid() const;
    /// Get total flux of model.
    virtual double getFlux() const = 0;
    /// Get type of model.
    virtual char getModelType() const = 0;
    /// Get reference ID of model.
    unsigned long getID() const;
  protected:
    /// Support area.
    Rectangle<double> support;
    /// Centroid.
    Point<double> centroid;
    /// Coordinate transformation for all calls to getValue().
    boost::shared_ptr<CoordinateTransformation> ct;
    /// Reference id.
    unsigned long id;
    /// Compute rectangular SourceModel::support for elliptical sources.
    /// Considers position of SourceModel::centroid.
    void setEllipticalSupport(double radius, const std::complex<double>& eps);
  };

  /// Populated an object by sampling the SourceModel.
  /// The coordinates of \obj are taken from its Grid, hence for a perfectly
  /// centered image of size \p N, the grid should have this form:
  /// \code
  /// obj.grid.setSize(-N/2,-N/2,N,N);
  /// \endcode
  /// The model can be oversampled by the factor \p S, such that each pixel
  /// is split into \f$S\times S\f$ subpixels, and the model is evaluated at the
  /// center of each subpixel.
  void setObject(const SourceModel& model, Object& obj, int S = 1);
 
  /// Collection of SourceModel entities.
  class SourceModelList : public std::vector<boost::shared_ptr<SourceModel> > {
  public:
    /// Create catalog from all entries of SourceModelList.
    Catalog getCatalog() const;
  };

  /// Sersic model class.
  /// The model has the form
  /// \f[I_S\bigl((x,y)\bigl) = exp\Bigl\lbrace -b_n\Bigl[\Bigl(\frac{r}{R_e}\Bigr)^{1/n} -1\Bigr]\Bigr\rbrace\f]
  /// with \f$r=\sqrt{x^2 + y^2}\f$ and \f$b_n\f$ defined by the relation between the complete and the incompete Gamma function,
  /// \f[\Gamma(2n)=2\gamma(2n,b_n)\ \Rightarrow\ b_n\approx 2 n - 1/3 + 0.009876/n.\f]
  /// The ensure vanishing flux at large radii, the profile is truncated at 
  /// \f$5R_e\f$ and the appropriate value at that position is subtracted from 
  /// \f$I_S\f$.\n\n
  /// Details can be seen in Graham & Driver, 2005, PASA, 22, 118-127 and
  /// Prugniel & Simien, 1997, A&A, 321, 111.
  class SersicModel : public SourceModel {
  public:
    /// Constructor with Sersic index \p n, effective radius \p Re, \p flux, and
    /// intrinsic ellipticity \p eps.
    /// If <tt>truncation!=0</tt>, the profile is truncated at is this number 
    /// of radii.
    SersicModel(double n, double Re, double flux, std::complex<double> eps, double truncation = 0, const CoordinateTransformation* CT = NULL, unsigned long id=0);
    /// Sample model at \p P.
    virtual double getValue(const Point<double>& P) const;
    /// Get total flux of model.
    virtual double getFlux() const;
    /// Get type of model.
    virtual char getModelType() const;
    /// Whether model contains this Point.
    virtual bool contains(const Point<double>& P) const; 
  private:
    double n, Re, b,limit,flux,flux_limit,shear_norm,flux_scale;
    std::complex<double> eps;
  };

  /// Moffat model class.
  /// The model has the form
  ///\f[I_M\bigl((x,y)\bigl) = \bigl(1+\alpha r^2\bigr)^{-\beta}\ \text{with}\ r=\sqrt{x^2 + y^2}\ \text{and}\ \alpha = \frac{2^{1/\beta}-1}{(FWHM/2)^2}.\f]
  /// The ensure vanishing flux at large radii, the profile can be truncated at
  /// some \p FWHM and the appropriate value at that position is subtracted from 
  /// \f$I_M\f$.
  class MoffatModel : public SourceModel {
  public:
    /// Constructor with Moffat index \p beta, width \p FWHM, \p flux, and
    /// intrinsic ellipticity \p eps.
    /// If <tt>truncation!=0</tt>, the profile is truncated at is this number of \p FWHM.
    MoffatModel(double beta, double FWHM, double flux, std::complex<double> eps, double truncation = 0, const CoordinateTransformation* CT = NULL, unsigned long id=0);
    /// Sample model at \p P.
    virtual double getValue(const Point<double>& P) const;
    /// Get total flux of model.
    virtual double getFlux() const;
    /// Get type of model.
    virtual char getModelType() const;
    /// Whether model contains this Point.
    virtual bool contains(const Point<double>& P) const; 
  private:
    double beta, alpha, limit, flux_limit, flux,shear_norm,flux_scale;
    std::complex<double> eps;
  };

  /// Pseudo-Airy model class.
  /// The model has the form
  /// \f[I_M\bigl((x,y)\bigl) = \begin{cases}\frac{\sin^2(r/r_d)}{(r/r_d)^2} & \text{if}\ \ r < r_d\\\frac{\sin^2(r/r_d)}{(r/r_d)^3} & \text{if}\ \ r > r_d\end{cases} \text{with}\ r=\sqrt{x^2 + y^2}\ \text{and}\ r_d= 0.5 * FWHM / 1.203\f]
  class AiryModel : public SourceModel {
  public:
    /// Constructor.
    /// If <tt>truncation!=0</tt>, the profile is truncated at is this number of \p FWHM.
    AiryModel(double FWHM, double flux, std::complex<double> eps, double truncation = 0, const CoordinateTransformation* CT = NULL, unsigned long id=0);
    /// Sample model at \p P.
    virtual double getValue(const Point<double>& P) const;
    /// Get total flux of model.
    virtual double getFlux() const;
    /// Get type of model.
    virtual char getModelType() const;
    /// Whether model contains this Point.
    virtual bool contains(const Point<double>& P) const; 
  private:
    double r_d, limit, flux_limit, flux, shear_norm, flux_scale;
    std::complex<double> eps;
  };


  /// Model from interpolated pixel data.
  /// The class provides several interpolation types for the Image given 
  /// at construction time.
  class InterpolatedModel : public SourceModel {
  public:
    /// Constructor.
    /// \p order defines order of interpolation.
    /// - <tt>1</tt>: bi-linear
    /// - <tt>n > 1</tt>: polynomial
    /// - <tt>-3</tt>: bi-cubic
    ///
    /// For more details, see Interpolation.
    InterpolatedModel(const boost::shared_ptr<Object>& obj, double flux, const CoordinateTransformation* CT = NULL, int order = 1, unsigned long id=0);
    /// Sample model at \p P.
    virtual double getValue(const Point<double>& P) const;
    /// Get total flux of model.
    virtual double getFlux() const;
    /// Get type of model.
    virtual char getModelType() const;
  private:
    boost::shared_ptr<Object> obj;
    int order;
    double flux,flux_scale;
    Point<double> reference;
  };


} // end namespace
#endif
