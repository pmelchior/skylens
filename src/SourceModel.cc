#include "../include/SourceModel.h"
#include "../include/Interpolation.h"
#include <shapelens/Moments.h>
#include <shapelens/MathHelper.h>

namespace skylens {

  SourceModel::~SourceModel() {
  }

  const Rectangle<double>& SourceModel::getSupport() const {
    return support;
  }

  const Point<double>& SourceModel::getCentroid() const {
    return centroid;
  }

  unsigned long SourceModel::getID() const {
    return id;
  }

  void SourceModel::setEllipticalSupport(double radius, const complex<double>& eps) {
    // compute orientation angle
    double theta;
    if (real(eps)!=0)
      theta = 0.5*atan(imag(eps)/real(eps));
    else
      theta = 0;
    // theta can not distinguish between vertical and horizontal orientation:
    // theta = 0 in both cases
    // since we want to have angle to x-Axis, map theta onto unbroken range of 180 deg
    if (real(eps) < 0)
      theta += M_PI_2;
    // compute size of semi-major and semi-minor axis (axis-parallel system)
    double a = (1 + abs(eps))*radius;
    double b = (1 - abs(eps))*radius;
    // compute curve parameter t which maximizes x or y
    double tx = atan(-b*tan(theta)/a);
    double ty = atan(b/(tan(theta)*a));
    // insert in parametric equation for rotated ellipse
    double max_x = fabs(a*cos(tx)*cos(theta) - b*sin(tx)*sin(theta));
    double max_y = fabs(a*cos(ty)*sin(theta) + b*sin(ty)*cos(theta));

    // lower-left
    support.ll(0) = -max_x;
    support.ll(1) = -max_y;
    support.tr(0) = max_x;
    support.tr(1) = max_y;
  }

  bool SourceModel::contains(const Point<double>& P) const {
    return support.contains(P);
  }

  void setObject(const SourceModel& model, Object& obj, int S) {
    double offset = 1./S; // subpixel offsets
    for (unsigned int i=0; i < obj.size(); i++) {
      obj(i) = 0;
      Point<double> P = obj.grid(i), P_;
      for (int n1 = 0; n1 < S; n1++) {
	P_(0) = P(0) + (0.5+n1)*offset - 0.5;
	for (int n2 = 0; n2 < S; n2++) {
	  P_(1) = P(1) + (0.5+n2)*offset - 0.5;
	  obj(i) += model.getValue(P_);
	}
      }
      obj(i) /= S*S;
    }
  }

  // ##### SourceModelList ##### //
  Catalog SourceModelList::getCatalog() const {
    Catalog cat;
    CatObject co;
    co.FLAGS = 0;
    for (unsigned long i=0; i < SourceModelList::size(); i++) {
      const SourceModel& sm = *SourceModelList::operator[](i);
      const Rectangle<double>& support = sm.getSupport();
      const Point<double>& centroid = sm.getCentroid();
      co.XMIN = (int) floor(support.ll(0));
      co.YMIN = (int) floor(support.ll(1));
      co.XMAX = (int) ceil(support.tr(0));
      co.YMAX = (int) ceil(support.tr(1));
      co.XCENTROID = centroid(0);
      co.YCENTROID = centroid(1);
      co.OPT["MODEL_TYPE"] = int(sm.getModelType());
      co.OPT["MODEL_ID"] = int(sm.getID());
      cat[i] = co;
    }
    return cat;
  }

  double fasterPow(double x, double y) {
    return exp(y*log(x));
  }

  // ##### Sersic Model ##### //
  SersicModel::SersicModel(double n, double Re, double flux_eff, complex<double> eps, double truncation, const CoordinateTransformation* ct_, unsigned long id) : 
    n(n), Re(Re), eps(eps) {
    if (truncation > 0)
      limit = truncation*Re;
    else
      limit = 1000*Re; // much larger than he image
    shear_norm = 1 - pow2(abs(eps));
    // compute WC of centroid (which is 0/0 in image coords)
    SourceModel::centroid(0) = 0;
    SourceModel::centroid(1) = 0;
    // compute support size from Re and eps
    SourceModel::setEllipticalSupport(limit,eps);
    // set the WCS from CT
    if (ct_!=NULL) {
      ct = ct_->clone();
      ct->transform(SourceModel::centroid);
      SourceModel::support.apply(*ct);
    }
    SourceModel::id = id;

    // Prugniel & Simien (1997), eq. A3a
    // limit for b > 0: n > 0.1282, reinforcing this limit
    if (n < 0.1282) {
      std::cerr << "SersicModel run with n_s = " << n << ", adjusted to n_s = 0.1282" << std::endl;
      n = 0.1282;
    }
    b = 2*n - 1./3 + 0.009876/n;
    double RRe1n = pow(limit/Re,1./n);
    // flux at limit
    flux_limit = exp(-b*(RRe1n -1));
    // compute total flux of model (considering the truncation at limit)
    flux = (Re*Re*2*M_PI*n*exp(b)/pow(b,2*n) * (gsl_sf_gamma(2*n) - gsl_sf_gamma_inc(2*n,b*RRe1n)));
    // subtract level at limit such that the profile vanishes there
    flux -= M_PI*limit*limit*flux_limit;
    // correct for shearing
    flux *= shear_norm;
    // compute rescaling factor for flux
    flux_scale = flux_eff/flux;
  }

  double SersicModel::getValue(const Point<double>& P) const {
    // get image coords from WC
    Point<double> P_ = P;
    if (ct.use_count() != 0)
      ct->inverse_transform(P_);

    // additionally apply shear transformation for an elliptical profile
    double x_ = (1-real(eps))*P_(0) - imag(eps)*P_(1);
    double y_ = -imag(eps)*P_(0) + (1+real(eps))*P_(1);
    
    double radius = sqrt(x_*x_ + y_*y_)/shear_norm; // shear changes size
    if (radius < limit)
      return flux_scale*(exp(-b*(fasterPow(radius/Re,1./n) -1)) - flux_limit);
    else
      return 0;
  }

  double SersicModel::getFlux() const {
    return flux_scale*flux;
  }

  char SersicModel::getModelType() const {
    return 0;
  }

  bool SersicModel::contains(const Point<double>& P) const {
    //std::cout << "contains " << P << "?\t";
    // get image coords from WC
    Point<double> P_ = P;
    if (ct.use_count() != 0)
      ct->inverse_transform(P_);

    // additionally apply shear transformation for an elliptical profile
    double x_ = (1-real(eps))*P_(0) - imag(eps)*P_(1);
    double y_ = -imag(eps)*P_(0) + (1+real(eps))*P_(1);
    
    double radius = sqrt(x_*x_ + y_*y_)/shear_norm; // shear changes size
    if (radius < limit) {
      //std::cout << "YES" << std::endl;
      return true;
    }
    else {
      //std::cout << "NO" << std::endl;
      return false;
    }
  }

  // ##### Moffat Model ##### //
  MoffatModel::MoffatModel(double beta, double FWHM, double flux_eff, complex<double> eps, double truncation, const CoordinateTransformation* ct_, unsigned long id) :
    beta(beta), eps(eps) {

    alpha = (pow(2.,1./beta)-1)/(FWHM*FWHM/4);
    if (truncation > 0) 
      limit = truncation*FWHM;
    else
      limit = 1000*FWHM; // much larger than any image
    shear_norm = 1 - abs(eps)*abs(eps);

    // compute WC of centroid (which is 0/0 in image coords)
    SourceModel::centroid(0) = 0;
    SourceModel::centroid(1) = 0;
    // compute support size from Re and eps
    SourceModel::setEllipticalSupport(limit,eps);
    // set the WCS from CT
    if (ct_!=NULL) {
      ct = ct_->clone();
      ct->transform(SourceModel::centroid);
      SourceModel::support.apply(*ct);
    }
    SourceModel::id = id;

    // flux at limit
    flux_limit = pow(1+alpha*limit*limit,-beta);
    // compute total flux of model (considering the truncation at limit)
    flux = 2*M_PI*(-1 + pow(1+alpha*limit*limit,1-beta))/(2*alpha - 2*alpha*beta);
    // subtract level at R such that the profile vanishes there
    flux -= M_PI*limit*limit*flux_limit;
    // correct for shearing
    flux *= shear_norm;
    // compute rescaling factor for flux
    flux_scale = flux_eff/flux;
  }

  double MoffatModel::getValue(const Point<double>& P) const {
    // get image coords from WC
    Point<double> P_ = P;
    if (ct.use_count() != 0)
      ct->inverse_transform(P_);

    // additionally apply shear transformation for an elliptical profile
    double x_ = (1-real(eps))*P_(0) - imag(eps)*P_(1);
    double y_ = -imag(eps)*P_(0) + (1+real(eps))*P_(1);
  
    double radius = sqrt(x_*x_ + y_*y_)/shear_norm;
    if (radius < limit)
      return flux_scale*(pow(1+alpha*radius*radius,-beta) - flux_limit);
    else
      return 0;
  }

  double MoffatModel::getFlux() const {
    return flux_scale*flux;
  }

  char MoffatModel::getModelType() const {
    return 3;
  }

  bool MoffatModel::contains(const Point<double>& P) const {
    // get image coords from WC
    Point<double> P_ = P;
    if (ct.use_count() != 0)
      ct->inverse_transform(P_);

    // additionally apply shear transformation for an elliptical profile
    double x_ = (1-real(eps))*P_(0) - imag(eps)*P_(1);
    double y_ = -imag(eps)*P_(0) + (1+real(eps))*P_(1);
    
    double radius = sqrt(x_*x_ + y_*y_)/shear_norm; // shear changes size
    if (radius < limit)
      return true;
    else
      return false;
  }

  // ##### Psuedo-Airy Model ##### //
  AiryModel::AiryModel(double FWHM, double flux_eff, complex<double> eps, double truncation, const CoordinateTransformation* ct_, unsigned long id) :
  eps(eps) {
    
    if (truncation > 0) 
      limit = truncation*FWHM;
    else
      limit = 1000*FWHM; // much larger than any image
    shear_norm = 1 - abs(eps)*abs(eps);
    
    // compute WC of centroid (which is 0/0 in image coords)
    SourceModel::centroid(0) = 0;
    SourceModel::centroid(1) = 0;
    // compute support size from Re and eps
    SourceModel::setEllipticalSupport(limit,eps);
    // set the WCS from CT
    if (ct_!=NULL) {
      ct = ct_->clone();
      ct->transform(SourceModel::centroid);
      SourceModel::support.apply(*ct);
    }
    SourceModel::id = id;
    // scale radius, see GREAT10 documentation
    r_d = 0.5 * FWHM / 1.203;
    // FIXME: flux and flux_limit not correctly set
    flux_limit = 0.;
    flux = 2 * M_PI;
    flux *= shear_norm;
    flux_scale = flux_eff/flux;
  }
  
  double AiryModel::getValue(const Point<double>& P) const {
    // get image coords from WC
    Point<double> P_ = P;
    if (ct.use_count() != 0)
      ct->inverse_transform(P_);
    
    // additionally apply shear transformation for an elliptical profile
    double x_ = (1-real(eps))*P_(0) - imag(eps)*P_(1);
    double y_ = -imag(eps)*P_(0) + (1+real(eps))*P_(1);
    
    double radius = sqrt(x_*x_ + y_*y_)/shear_norm/r_d;
    double sr = sin(radius);
    double temp = sr/radius;
    if (radius < 1e-3) // singularity at zero
      temp = 1;

    if (radius < 1.)
      return flux_scale*temp*temp;
    else
      return flux_scale*temp*temp/radius;
  }
  
  double AiryModel::getFlux() const {
    return 2 * M_PI;
  }
  
  char AiryModel::getModelType() const {
    return 4;
  }

  bool AiryModel::contains(const Point<double>& P) const {
    // get image coords from WC
    Point<double> P_ = P;
    if (ct.use_count() != 0)
      ct->inverse_transform(P_);

    // additionally apply shear transformation for an elliptical profile
    double x_ = (1-real(eps))*P_(0) - imag(eps)*P_(1);
    double y_ = -imag(eps)*P_(0) + (1+real(eps))*P_(1);
    
    double radius = sqrt(x_*x_ + y_*y_)/shear_norm; // shear changes size
    if (radius < limit)
      return true;
    else
      return false;
  }


  // ##### Interpolated Model ##### //
  InterpolatedModel::InterpolatedModel(const boost::shared_ptr<Object>& obj_, double flux, const CoordinateTransformation* ct_, int order, unsigned long id) : 
    obj(obj_), order(order),flux(flux) {

    // compute WC of centroid (which is 0/0 in image coords)
    SourceModel::centroid(0) = 0;
    SourceModel::centroid(1) = 0;
    SourceModel::support = obj->grid.getSupport().getBoundingBox();
    // account of centroid offset of obj
    SourceModel::support -= obj->centroid;
    if (ct_ != NULL) {
      ct = ct_->clone();
      ct->transform(SourceModel::centroid);
      SourceModel::support.apply(*ct);
    }
    SourceModel::id = id;
  
    // compute rescaling factor for flux
    Moments objflux(*obj,FlatWeightFunction(), 0); 
    flux_scale = flux/objflux(0);
  }

  double InterpolatedModel::getValue(const Point<double>& P) const {
    // no check here if P is in support, because interpolation returns
    // zero anyway in this case...
    // get image coords from WC
    Point<double> P_ = P;
    if (ct.use_count() != 0)
      ct->inverse_transform(P_);

    // account of centroid offset of obj
    P_ += obj->centroid;
    switch (order) {
    case 1: // simple bi-linear interpolation
      return flux_scale*obj->interpolate(P_);
    case -3: // bi-cubic interpolation
      return flux_scale*Interpolation::bicubic(*obj,P_);
    default: // nth-order polynomial interpolation
      return flux_scale*Interpolation::polynomial(*obj,P_,order);
    }
  }

  double InterpolatedModel::getFlux() const {
    return flux;
  }

  char InterpolatedModel::getModelType() const {
    return 2;
  }

} // end namespace
