#include "../include/Layer.h"
#include "../include/RNG.h"
#include <shapelens/FITS.h>
#include <shapelens/MathHelper.h>

#ifdef HAS_OpenMP
#include <omp.h>
#endif

using namespace skylens;
using std::complex;

LensingLayer::LensingLayer(double z_, std::string angle_file, const shapelens::Point<double>* center) :
  // automatically creates a single instance of LayerStack
  ls(SingleLayerStack::getInstance()),
  li(SingleLensingInformation::getInstance()),
  cosmo(SingleCosmology::getInstance())
{
  Layer::z = z_;
  Layer::transparent = false;
  me = ls.insert(std::pair<double,Layer*>(z,this));

#pragma omp critical
  {
  // check if this is the first lensing layer
  li.z_first_lens = std::min(li.z_first_lens, z);
  // compute c/H0: units of D
  li.c_H0 = cosmo.getc()/cosmo.getH0()*cosmo.h100/cosmo.getMpc();  
  }

  // open file with two real-valued images of first and second component
  fitsfile* fptr = shapelens::FITS::openFile(angle_file);
  // read in lens parameters
  double sidelength, z_lens, z_source, h, omega, lambda;
  shapelens::FITS::readKeyword(fptr,"SIDEL",sidelength);
  shapelens::FITS::readKeyword(fptr,"ZLENS",z_lens);
  shapelens::FITS::readKeyword(fptr,"ZSOURCE",z_source);
  shapelens::FITS::readKeyword(fptr,"OMEGA",omega);
  shapelens::FITS::readKeyword(fptr,"LAMBDA",lambda);
  shapelens::FITS::readKeyword(fptr,"H",h);
  // compute rescaling factor of deflection angle
  // in the cosmology from the FITS file
  Cosmology cosmo_l(omega,lambda,h);
  double Dl, Dls, Ds, c_H0;
  // D in units [c/H0] = [cm] -> [Mpc/h]
  c_H0 = cosmo_l.getc()/cosmo_l.getH0()*h/cosmo_l.getMpc();
  Dl = cosmo_l.Dang(z_lens)*c_H0;
  Ds = cosmo_l.Dang(z_source)*c_H0;
  Dls = cosmo_l.Dang(z_source,z_lens)*c_H0;
  // take out lensing efficiency factor (which is reinserted in getFlux())
  // and convert to arcsec
  scale0 = Ds/Dls * 180/M_PI * 3600;
  
  // Massimo's convention for deflection angles contains 
  // rescaling factor of sidelength/Dl
  // Do we need to apply it: Yes if keyword LENSRESC is set
  bool lensresc;
  try {
    shapelens::FITS::readKeyword(fptr,"LENSRESC",lensresc);
    if (lensresc)
      scale0 *= sidelength/Dl;
  } catch (std::exception) {}
  
  // read in complex deflection angle field
  Image<float> a1;
  shapelens::FITS::readImage(fptr, a1);
  a.resize(a1.size());
  a.grid = a1.grid;
  a.realPart() = a1;
  shapelens::FITS::moveToExtension(fptr, 2);
  shapelens::FITS::readImage(fptr, a1);
  a.imagPart() = a1;

  // compute angular rescaling factor: 
  // arcsec -> pixel position in angle map
  // if the lens needs to be moved, we must set it here!
  theta0 = (sidelength/Dl)*(180/M_PI)*3600/a.grid.getSize(0);
  shapelens::ScalarTransformation S(theta0);

  shapelens::Point<double> center_image(-0.5*a.grid.getSize(0),-0.5*a.grid.getSize(1));
  shapelens::ShiftTransformation Z(center_image);
  if (center != NULL) {
    shapelens::ShiftTransformation ZF(*center);
    S *= ZF;
    Z *= S;
    a.grid.setWCS(Z);
  } else {
    Z *= S;
    a.grid.setWCS(Z);
  }
  shapelens::FITS::closeFile(fptr);
}

// iter points to current lensing layer
void LensingLayer::setDistances(const LayerStack::iterator& iter_) const {
  LayerStack::iterator iter = iter_;
  iter++;
  std::string type;
  li.Dls.insert(std::pair<double, std::map<double,double> >(z, std::map<double,double>()));
  for (iter; iter != ls.end(); iter++) {
    type = iter->second->getType();
    if (type[0] == 'S') {
      li.Ds[iter->first] = cosmo.Dang(iter->first)*li.c_H0;
      li.Dls[iter_->first][iter->first] = cosmo.Dang(iter->first, iter_->first)*li.c_H0;
    }
    else if (type == "TL") {
      // recursively call setDistances
      setDistances(iter);
    }
  }
}

// sum all fluxes until the next transformation layer is found
double LensingLayer::getFlux(const shapelens::Point<double>& P, double* z_) const {
  double flux = 0;
  if (z_ != NULL)
    if (*z_ < z)
      return flux;

  LayerStack::iterator iter = me;
  iter++;
  std::string type;

  // first lensing layer must selectively switch on source layers
  if (li.z_first_lens == z) {
    if (z_ != NULL)
      throw std::invalid_argument("LensingLayer: getFlux() source redshift set before first layer!");

    // first time only: 
    // find source layers and compute their distances to redshift 0 and z_l
#pragma omp critical
    if (li.Ds.size() == 0) {
      setDistances(me);
    }

    // raytrace thru lensed source layers: one at a time
    for (std::map<double, double>::iterator source_iter = li.Ds.begin(); source_iter != li.Ds.end(); source_iter++) {
      // active source layer 
      double zs = source_iter->first;

      // apply lens equation: beta = theta - Dls/Ds*alpha(theta)
      complex<float> p(P(0),P(1));
      if (!transparent)
	p -= a.interpolate(P) * scale0* float(li.Dls[z][zs] / source_iter->second);
      for (iter; iter != ls.end(); iter++) { 
	type = iter->second->getType();
	flux += iter->second->getFlux(shapelens::Point<double>(real(p),imag(p)), &zs);
	// since source_iter is the last and only shining
	// source layer, we can stop here
	if (type[0] == 'T' || iter->first == zs)
	  break;
      }
    }
  }
  
  // any farther lensing layer just computes the deflection angle for the
  // specified source redshify
  else {
    // apply lens equation:
    // beta = theta - Dls/Ds*alpha(theta)
    complex<float> p(P(0),P(1));
    if (!transparent)
      p -= a.interpolate(P) * scale0 * float(li.Dls[z][*z_] / li.Ds[*z_]);

    for (iter; iter != ls.end(); iter++) {
      std::string type = iter->second->getType();
      flux += iter->second->getFlux(shapelens::Point<double>(real(p),imag(p)), z_);
      // since given source is the last and only shining
      // source layer, we can stop here
      if (type[0] == 'T' || iter->first == *z_)
	break;
    }
  }
 
  return flux;
}

std::string LensingLayer::getType() const {
  return "TL";
}

shapelens::Point<double> LensingLayer::getCenter() const {
  shapelens::Point<int> pc(a.grid.getSize(0)/2,a.grid.getSize(1)/2);
  return a.grid(a.grid.getPixel(pc));
}


inline double linInt(const std::vector<double> phi, double tx, double ty) {
  return phi[0]*(1-tx)*(1-ty) + phi[1]*tx*(1-ty) + phi[2]*(1-tx)*ty + phi[3]*tx*ty;
}

std::vector<shapelens::Point<double> > LensingLayer::findCriticalPoints(double zs, int det_sign) const {
  double D_ls = cosmo.Dang(zs,z);
  double D_s = cosmo.Dang(zs);
  std::vector<shapelens::Point<double> > cpoints;
  Point<int> P;
  double phixx, phiyy, phixy, phiyx, kappa, gamma1, gamma2, gamma, lambda_t, lambda_r;
  double acc = 1e-1;
  int C = 10;

  // compute finite differences -> kappa, gamma -> eigenvalues of A
  // if below rather soft threshold of acc:
  // span a CxC subpixel grid and search for lambda < acc/C**2
  for (long i = 3; i < a.grid.getSize(0) - 3; i++) {
    for (long j = 3; j < a.grid.getSize(1) - 3; j++) {
      P(0) = i;
      P(1) = j;
      finiteDifferences(P, phixx, phixy, phiyx, phiyy);
      kappa = 0.5*(phixx + phiyy) * scale0 * D_ls / D_s;
      gamma1 = 0.5*(phixx - phiyy) * scale0 * D_ls / D_s;
      gamma2 = phixy * scale0 * D_ls / D_s;
      gamma = sqrt(gamma1*gamma1+gamma2*gamma2);
      // double jacdet = (1 - kappa) * (1 - kappa) - gamma*gamma;
      lambda_t = 1 - kappa - gamma; // tangential eigenvalue
      lambda_r = 1 - kappa + gamma; // radial eigenvalue
      if (fabs(lambda_t) < acc || fabs(lambda_r) < acc) {
	
	std::vector<double> phixx_(4), phixy_(4), phiyy_(4);
	phixx_[0] = phixx;
	phixy_[0] = phixy;
	phiyy_[0] = phiyy;
	P(0) += 1;
	finiteDifferences(P, phixx_[1], phixy_[1], phixy_[1], phiyy_[1]);
	P(0) -= 1;
	P(1) += 1;
	finiteDifferences(P, phixx_[2], phixy_[2], phixy_[2], phiyy_[2]);
	P(0) += 1;
	finiteDifferences(P, phixx_[3], phixy_[3], phixy_[3], phiyy_[3]);

	Point<double> critical;
	for (double tx=1./(2*C); tx < 1; tx += 1./C) { 
	  for (double ty=1./(2*C); ty < 1; ty += 1./C) { 

	    phixx = linInt(phixx_, tx, ty);
	    phixy = linInt(phixy_, tx, ty);
	    phiyy = linInt(phiyy_, tx, ty);

	    kappa = 0.5*(phixx + phiyy) * scale0 * D_ls / D_s;
	    gamma1 = 0.5*(phixx - phiyy) * scale0 * D_ls / D_s;
	    gamma2 = phixy * scale0 * D_ls / D_s;
	    gamma = sqrt(gamma1*gamma1+gamma2*gamma2);
	    lambda_t = 1 - kappa - gamma;
	    lambda_r = 1 - kappa + gamma;
	    if ((lambda_t*lambda_r < acc*acc/pow4(C) && det_sign == 0) || (fabs(lambda_t) < acc/C*C && det_sign > 0) || (fabs(lambda_r) < acc/C*C && det_sign < 0)) {
	      critical(0) = i + tx;
	      critical(1) = j + ty;
	      a.grid.getWCS().transform(critical);
	      cpoints.push_back(critical);
	    }
	  }
	}
      }
    }
  }
  return cpoints;
}

void LensingLayer::finiteDifferences(const Point<int>& P0, double& phixx, double& phixy, double& phiyx, double& phiyy) const {
  // centered finite difference for first derivatives of 6th order
  int i = P0(0), j = P0(1);
  phixx = (real(a(i+3,j)) - real(a(i-3,j))  - 9*(real(a(i+2,j)) - real(a(i-2,j))) + 45*(real(a(i+1,j)) - real(a(i-1,j))))/(60*theta0);
  phiyx = (imag(a(i+3,j)) - imag(a(i-3,j))  - 9*(imag(a(i+2,j)) - imag(a(i-2,j))) + 45*(imag(a(i+1,j)) - imag(a(i-1,j))))/(60*theta0);
  phiyy = (imag(a(i,j+3)) - imag(a(i,j-3))  - 9*(imag(a(i,j+2)) - imag(a(i,j-2))) + 45*(imag(a(i,j+1)) - imag(a(i,j-1))))/(60*theta0);

  // for efficiency, undo this if you want to test accuracy of interpolation
  phixy = phiyx;
}

void LensingLayer::set_Dphi(const Point<double>& theta, double zs, double& phixx, double& phixy, double& phiyx, double& phiyy) const {
  double D_ls = cosmo.Dang(zs,z);
  double D_s = cosmo.Dang(zs);

  Point<int> P0 = a.grid.getCoords(theta); // lower-left point in grid of a
  Point<double> P = theta;
  a.grid.getWCS().inverse_transform(P);    // theta in pixel coordinates of a
  double tx = (P(0)-P0(0));
  double ty = (P(1)-P0(1));

  std::vector<double> phixx_(4), phixy_(4), phiyy_(4);
  finiteDifferences(P0, phixx_[0], phixy_[0], phixy_[0], phiyy_[0]);
  P0(0) += 1;
  finiteDifferences(P0, phixx_[1], phixy_[1], phixy_[1], phiyy_[1]);
  P0(0) -= 1;
  P0(1) += 1;
  finiteDifferences(P0, phixx_[2], phixy_[2], phixy_[2], phiyy_[2]);
  P0(0) += 1;
  finiteDifferences(P0, phixx_[3], phixy_[3], phixy_[3], phiyy_[3]);

  phixx = linInt(phixx_, tx, ty) * scale0 * D_ls / D_s;
  phiyy = linInt(phiyy_, tx, ty) * scale0 * D_ls / D_s;
  phixy = phiyx = linInt(phixy_, tx, ty) * scale0 * D_ls / D_s;
}


std::complex<double> LensingLayer::getShear(const Point<double>& theta, double zs, bool reduced) const {
  double phixx, phiyy, phixy, phiyx;
  set_Dphi(theta, zs, phixx, phixy, phiyx, phiyy);
  complex<double> gamma(0.5*(phixx - phiyy), phixy);
  if (reduced) {
    double kappa = 0.5*(phixx + phiyy);
    gamma /= 1-kappa;
  }
  return gamma;
}

double LensingLayer::getConvergence(const Point<double>& theta, double zs) const {
 double phixx, phiyy, phixy, phiyx;
 set_Dphi(theta, zs, phixx, phixy, phiyx, phiyy);
 return 0.5*(phixx + phiyy);
}


Point<double> LensingLayer::getBeta(const Point<double>& theta, double zs) const {
  // speed this up eventually
  double D_ls = cosmo.Dang(zs,z);
  double D_s = cosmo.Dang(zs);
  complex<float> p(theta(0),theta(1));
  p -= a.interpolate(theta) * scale0 * float(D_ls / D_s);
  return Point<double>(real(p), imag(p));
}


// find cells in a image-plane grid, whose source-plane outline encloses
// the desired point beta. To avoid the ambiguity of non-simple polygons
// (to determine whether point is inside), cells are split into triangles.
// However, since the mapping can distort the shape of those triangles,
// it is possible that no match is found, especcially if the point is close
// to the edge. Therefore iterate with slighly shifted search grids 
// if necessary.
// Heavily inspired by from Matthias Bartelmann's libastro
std::list<Rectangle<double> > LensingLayer::getCellsEnclosing(const Point<double>& beta, double zs, const Rectangle<double>& area, int C) const {
  Point<double> theta, theta_, beta_;
  std::list<Rectangle<double> > cells;
  Rectangle<double> cell;
  double cellsize0 = (area.tr(0)-area.ll(0))/C, cellsize1 = (area.tr(1)-area.ll(1))/C;
  for (theta(0) = area.ll(0); theta(0) < area.tr(0); theta(0) += cellsize0) {
    for (theta(1) = area.ll(1); theta(1) < area.tr(1); theta(1) += cellsize1) {
      theta_ = theta;
      beta_ = getBeta(theta_, zs);
      double d11 = beta(0) - beta_(0), d21 = beta(1) - beta_(1);
      theta_(0) += cellsize0;
      beta_ = getBeta(theta_, zs);
      double d12 = beta(0) - beta_(0), d22 = beta(1) - beta_(1);
      theta_(0) -= cellsize0;
      theta_(1) += cellsize1;
      beta_ = getBeta(theta_, zs);
      double d13 = beta(0) - beta_(0), d23 = beta(1) - beta_(1);
      theta_(0) += cellsize0;
      beta_ = getBeta(theta_, zs);
      double d14 = beta(0) - beta_(0), d24 = beta(1) - beta_(1);

      // cross product of two rectangles inside search cell
      double c11=d11*d23-d13*d21;
      double c12=d13*d24-d14*d23;
      double c13=d14*d21-d11*d24;
      
      double c21=d11*d22-d12*d21;
      double c22=d12*d24-d14*d22;
      double c23=d14*d21-d11*d24;

      double p11=c11*c12;
      double p12=c12*c13;
      double p13=c13*c11;

      double p21=c21*c22;
      double p22=c22*c23;
      double p23=c23*c21;

      bool l1=p11>0.0 && p12>0.0 && p13>0.0;
      bool l2=p21>0.0 && p22>0.0 && p23>0.0;

      // if beta is in either triangle spanned by the new search cell: bingo!
      if (l1 || l2) {
	cell.ll = theta;
	cell.tr = theta_;
	cells.push_back(cell);
      }
    }
  }
  return cells;
}

std::vector<Point<double> > LensingLayer::findImages(const Point<double>& beta, double zs, const Rectangle<double>& area) const {
  // set up initial search grid of 10x10 cells
  //Rectangle<double> bbox = a.grid.getSupport().getBoundingBox();
  int C = 10, level = 1;
  std::list<Rectangle<double> > cells = getCellsEnclosing(beta, zs, area, C);
  std::vector<Point<double> > thetas; // multiple solutions possible

  // for random displacements
  RNG rng;// = Singleton<RNG>::getInstance();
  const gsl_rng * r = rng.getRNG();

  if (cells.size()) { // if not: deflection angles too large: not in bbox
                      // thetas empty then...
    // cell size larger than lens map pixels
    while (pow_int(C, level) < a.grid.getSize(0)) { 

      // next level: C x C cells in each matching cell before
      std::list<Rectangle<double> > cells_;
      for (std::list<Rectangle<double> >::const_iterator iter = cells.begin(); iter!= cells.end(); iter++) {
	Rectangle<double> area_ = *iter;
	std::list<Rectangle<double> > cells__;
	do { // if point cannot be located, it's close to boundary of cell:
	     // move the cell around
	  cells__ = getCellsEnclosing(beta, zs, area_, C);
	  Point<double> delta((-0.5 + gsl_rng_uniform (r))*(iter->tr(0) - iter->ll(0)),
			      (-0.5 + gsl_rng_uniform (r))*(iter->tr(1) - iter->ll(1)));
	  area_.ll = iter->ll + delta;
	  area_.tr = iter->tr + delta;
	} while (cells__.size() == 0); 
	cells_.insert(cells_.end(), cells__.begin(), cells__.end());
      }
      level += 1;
      cells = cells_;
    }

    // get centers of found cells
    for (std::list<Rectangle<double> >::const_iterator iter = cells.begin(); iter!= cells.end(); iter++)
      thetas.push_back(Point<double> ((iter->tr(0) + iter->ll(0))/2,
				      (iter->tr(1) + iter->ll(1))/2));
  }
  return thetas;
}
