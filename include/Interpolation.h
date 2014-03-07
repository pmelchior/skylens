#ifndef SKYLENS_INTERPOLATION_H
#define SKYLENS_INTERPOLATION_H

#include "Singleton.h"
#include <TMV.h>
#include <gsl/gsl_sf.h>
#include <stdexcept>
#include <shapelens/ShapeLens.h>
#include <shapelens/Image.h>
#include <shapelens/MathHelper.h>

namespace skylens {
  // helper class: uses a singleton for matrix w in bicubic
  // hence this matrix is only constructed once.
  template <class T>
  class BicubicW : public tmv::Matrix<T> {
  public:
  BicubicW() : tmv::Matrix<T>(16,16) {
      T wt[256] = 
	{ 1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	  0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,
	  -3,0,0,3,0,0,0,0,-2,0,0,-1,0,0,0,0,
	  2,0,0,-2,0,0,0,0,1,0,0,1,0,0,0,0,
	  0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,
	  0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,
	  0,0,0,0,-3,0,0,3,0,0,0,0,-2,0,0,-1,
	  0,0,0,0,2,0,0,-2,0,0,0,0,1,0,0,1,
	  -3,3,0,0,-2,-1,0,0,0,0,0,0,0,0,0,0,
	  0,0,0,0,0,0,0,0,-3,3,0,0,-2,-1,0,0,
	  9,-9,9,-9,6,3,-3,-6,6,-6,-3,3,4,2,1,2,
	  -6,6,-6,6,-4,-2,2,4,-3,3,3,-3,-2,-1,-1,-2,
	  2,-2,0,0,1,1,0,0,0,0,0,0,0,0,0,0,
	  0,0,0,0,0,0,0,0,2,-2,0,0,1,1,0,0,
	  -6,6,-6,6,-3,-3,3,3,-4,4,2,-2,-2,-2,-1,-1,
	  4,-4,4,-4,2,2,-2,-2,2,-2,-2,2,1,1,1,1};
      // copying this is not the fastest thing to do...
      // .. but we do it only once (w is static)
      for (int i = 0; i < 256; i++)
	*(tmv::Matrix<T>::ptr()+i) = *(wt+i);
    }
  };

  // helper class: uses singleton for weights w_j in polynomial
  class PolynomialW : public tmv::Vector<data_t> {
  public:
  PolynomialW() : tmv::Vector<data_t>(2) {
      setOrder(1);
    }
    void setOrder(int n_) {
      if (n_ != n) {
	n = n_;
	tmv::Vector<data_t>::resize(n+1);
	for (int j=0; j<=n; j++)
	  tmv::Vector<data_t>::operator()(j) = pow_int(-1,j)*(gsl_sf_fact(n)/(gsl_sf_fact(j)*gsl_sf_fact(n-j)));
      }
    }
    private:
	int n;
  };

  /// Class for various types of Interpolation.
  /// Implemented are 
  /// - polynomial interpolation (using the barycentric formula
  ///   with complexity \f$O(n)\f$ 
  ///   (see http://web.comlab.ox.ac.uk/people/Nick.Trefethen/barycentric.pdf)
  ///   of tmv::Vector and Image data.
  /// - Bi-cubic interpolation of Image data (see Numerical Recipes, chapter 3).
  ///
  class Interpolation {
  public:
    /// Polynomial interpolation of tmv::Vector data.
    /// Interpolates \p data at \p x using a polynomial of order \p n
    /// with barycentric formula
    /// \f[ p(x)=\frac{\sum_{j=0}^n \frac{w_j}{x-x_j} f_j}{\sum_{j=0}^n \frac{w_j}{x-x_j}}, \f]
    /// where \f$x_j\f$ are equidistantly spaced around \p x. In this case,
    /// \f$w_j=(-1)^j\binom{n}{j}\f$.
    template <class T>
      static T polynomial(const tmv::Vector<T>& data, data_t x, int n) {
      int x_ = (int) floor (x);
      if (x_ >= 0 && x_ < data.size()) {
	if (fabs(x_ - x) < 1e-10) // if x is too close to sampling point
	  return data(x_);        // use data from there instead
	else {
	  data_t k_j, denom = 0;
	  int x_j;
	  T num(0);
	  PolynomialW& w = Singleton<PolynomialW>::getInstance();
	  w.setOrder(n);
	  for (int j = 0; j <= n; j++) {
	    x_j = x_ + j - n/2; // x_j equidistantly sampled around x
	    if (x_j >=0 && x_j < data.size()) {
	      k_j = w(j)/(x-x_j);
	      num += k_j*data(x_j);
	      denom += k_j;
	    }
	  }
	  return num/denom;
	}
      } else
	throw std::invalid_argument("Interpolation: x is not within data");
    }
    /// Polynomial interpolation of tmv::Vector \p data with arbitrary 
    /// sampling points \p xx.
    /// \p xx must be an ordered list of sampling points, and \p data
    /// must be sampled there (and thus have the same ordering).\n
    /// Analogous to equidistant sampling apart from
    /// \f[ w_j = \frac{1}{\Pi_{k\neq j} (x_j-x_k)} \f].
    template <class T>
      static T polynomial(const tmv::Vector<T>& data, const tmv::Vector<data_t>& xx, data_t x, int n) {
      if (xx.size() != data.size())
	throw std::invalid_argument("Interpolation: data and xx have different size");
      if (x >= xx.minElement() && x <= xx.maxElement()) {
	int index = 0;
	data_t diff = x - xx(index);
	for (int i=1; i <= xx.size(); i++) {
	  if (fabs(x-xx(i)) < diff) {
	    diff = fabs(x-xx(i));
	    index = i;
	  }
	}
	if (fabs(xx(index) - x) < 1e-10)
	  return data(index);
	if (index < n/2 || data.size() - index < n/2)
	  throw std::invalid_argument("Interpolation: x is too close to boundary");
	else {
	  data_t denom = 0, k_j = 1; // k_j = w_j/(x-x_j);
	  T num(0);
	  for (int j = 0; j <= n; j++) {
	    for (int k = 0; k <= n; k++)
	      if (k != j)
		k_j*= xx(index + j-n/2)-xx(index + k-n/2);
	    k_j = 1./k_j;
	    k_j /= (x-xx(index + j-n/2));
	    num += k_j*data(index + j-n/2);
	    denom += k_j;
	  }
	  return num/denom;
	}
      } else
	throw std::invalid_argument("Interpolation: x is not within xx");
    }
    /// 2D polynomial interpolation of order \p n of \p im at 
    /// a World coordinate point \p P.
    /// Uses a series of 1D interpolations of order \p n.
    template <class T>
      static T polynomial(const Image<T>& im, const Point<data_t>& P, int n) {
      const Grid& grid = im.grid;
      Point<int> IC = grid.getCoords(P);
      if (grid.getPixel(IC) != -1) { // P inside im
	data_t x,y;
	const CoordinateTransformation& ct = grid.getWCS();
	Point<data_t> P_ = P;
	ct.inverse_transform(P_); // World -> pixel
	x = P_(0);
	y = P_(1);
	tmv::Vector<T> first_run(n+1),row(n+1);
	int x_ = (int) floor(x);
	int y_ = (int) floor(y);
	for (int i=0; i <= n; i++) {
	  int y__ = y_ + i - n/2;
	  for (int j=0; j <= n; j++)
	    row(j) = im.get(Point<int>(x_ + j - n/2, y__));
	  first_run(i) = polynomial(row,x-(x_ - n/2),n);
	}
	return polynomial(first_run,y-(y_ - n/2),n);
      } else
	throw std::invalid_argument("Interpolation: (x/y) is not within image");
    }
    /// Bi-cubic interpolation of \p im at a World coordinate point \p P.
    /// Uses Numerical Recipes \p bcucof and finite differences for gradients
    template <class T>
      static T bicubic(const Image<T>& im, const Point<data_t>& P) {
      const Grid& grid = im.grid;
      Point<int> IC = grid.getCoords(P);
      if (grid.getPixel(IC) != -1) { // P inside im
	data_t x,y;
	const CoordinateTransformation& ct = grid.getWCS();
	Point<data_t> P_ = P;
	ct.inverse_transform(P_); // World -> pixel
	x = P_(0);
	y = P_(1);
	int x_ = (int) floor(x);
	int y_ = (int) floor(y);
	tmv::Matrix<T>& w = Singleton<BicubicW<T> >::getInstance();

	tmv::Vector<T> yy(4), y1(4), y2(4), y12(4);
	for (int i=1; i < 5; i++) {
	  int x__, y__;
	  switch (i) { // numering scheme from Fig. 3.6.1
	  case 1: x__ = x_; y__ = y_; break;
	  case 2: x__ = x_ + 1; y__ = y_; break;
	  case 3: x__ = x_ + 1; y__ = y_ + 1; break;
	  case 4: x__ = x_; y__ = y_ + 1; break;
	  }
	  yy(i-1) = im.get(Point<int>(x__,y__));
	  // finite differences (centered)
	  y1(i-1) = 0.5*(im.get(Point<int>(x__+1,y__)) - 
			 im.get(Point<int>(x__-1,y__)));
	  y2(i-1) = 0.5*(im.get(Point<int>(x__,y__+1)) - 
			 im.get(Point<int>(x__,y__-1)));
	  y12(i-1) = 0.25*(im.get(Point<int>(x__+1,y__+1)) - 
			   im.get(Point<int>(x__+1,y__-1)) - 
			   im.get(Point<int>(x__-1,y__+1)) - 
			   im.get(Point<int>(x__-1,y__-1)));
	}
	tmv::Vector<T> xx(16);
	for (int i=0; i < 4; i++) {
	  xx(i)=yy(i);
	  xx(i+4)=y1(i);
	  xx(i+8)=y2(i);
	  xx(i+12)=y12(i);
	}
	tmv::Vector<T> cl = w*xx;
	int l=0;
	T p=0;
	for (int i=0; i < 4; i++) {
	  for (int j=0; j < 4; j++) {
	    // cl(l) = c(i,j)
	    p += cl(l)*pow_int(x-x_,i)*pow_int(y-y_,j); 
	    l++;
	  }
	}
	return p;
      } else
	return 0;//throw std::invalid_argument("Interpolation: (x/y) is not within image");
    }
  };
} // end namespace
#endif
