#include "../include/Filter.h"
#include "../include/Integrator.h"
#include <shapelens/FITS.h>
#include <shapelens/MathHelper.h>
#include <stdexcept>

namespace skylens {
  using shapelens::FITS;

  Filter::Filter() : prefactor(1), z(0) {}
  Filter::Filter(const std::string& filename, double threshold) :  prefactor(1), z(0) {
    fitsfile* fptr = FITS::openTable(filename);
    long nrows = FITS::getTableRows(fptr);
    int nu_col = FITS::getTableColumnNumber(fptr, "frequency");
    int fn_col = FITS::getTableColumnNumber(fptr, "transmission");
    double nu, fn;
    for (long i=0; i < nrows; i++) {
      FITS::readTableValue(fptr, i, nu_col, nu);
      FITS::readTableValue(fptr, i, fn_col, fn);
      Filter::insert(std::pair<double, double>(nu, fn));
    }
    removeZeros(threshold);
  }

  void Filter::save(const std::string& filename) const {
    // don't forget the write out prefactor*nu
  }
  
  Filter& Filter::operator*= (double c) {
    prefactor *= c;
    return *this;
  }

  Filter& Filter::operator/= (double c) {
    prefactor /= c;
    return *this;
  }

  Filter& Filter::operator*= (const Filter& f) {
    // non-overlapping filter: null filter
    if (f.getNuMin() > Filter::getNuMax() || f.getNuMax() < Filter::getNuMin()) {
      Filter fnull;
      fnull.insert(std::pair<double, double>(Filter::getNuMin(), 0));
      fnull.insert(std::pair<double, double>(Filter::getNuMax(), 0));
      *this = fnull;
    }
    else {
      // when filters overlap we need to multiply with f.avg() at each nu
      Filter::iterator iter = Filter::begin(), next, last = --Filter::end();
      while (iter != last) {
	next = iter;
	next++;
	iter->second *= f.avg(iter->first/(1+z), next->first/(1+z));
	iter++;
      }
    }
    return *this;
  }

  Filter& Filter::operator+= (const Filter& f) {
    // non-overlapping filter: just merge
    if (f.getNuMin() > Filter::getNuMax() || f.getNuMax() < Filter::getNuMin()) {
      for (Filter::const_iterator iter = f.begin(); iter != f.end(); iter++)
	Filter::insert(std::pair<double,double>(iter->first, f.prefactor*iter->second / prefactor)); // need to normalize the prefactors
    }
    else {
      // when filters overlap we need to add transmissions
      throw std::runtime_error("Filter::operator+= not implemented for overalpping filters!");
    }
    return *this;
  }

  double Filter::getNuMin() const {
    return Filter::begin()->first / (1+z);
  }

  double Filter::getNuMax() const {
    return Filter::rbegin()->first / (1+z);
  }
  
  // linear interpolation between elements bracketing nu
  double Filter::operator() (double nu) const {
    if (nu > getNuMax() || nu < getNuMin())
      return 0;
    Filter::const_iterator i = Filter::upper_bound(nu * (1+z));
    Filter::const_iterator l = i; 
    --l;
    const double delta = (nu - l->first/(1+z))/(i->first/(1+z) - l->first/(1+z));
    return prefactor*(delta*i->second + (1-delta)*l->second);
  }

  double Filter::avg(double numin, double numax) const {
    if (numin > numax) { // wrong order
      double tmp = numax;
      numax = numin;
      numin = tmp;
    }
    if (numin > getNuMax() || numax < getNuMin())
      return 0;
    else {
      double avg = 0, nu;
      // first part: between numin and its next upper sample i
      Filter::const_iterator i = Filter::upper_bound(numin*(1+z));
      Filter::const_iterator l = i; 
      --l;
      double delta = (numin - l->first/(1+z))/(i->first/(1+z) - l->first/(1+z));
      double fn = delta*i->second + (1-delta)*l->second;
      avg += 0.5*(i->first/(1+z) - numin)*(fn + i->second);
      // middle part: between i and lower sample of numax
      l = i;
      i = Filter::upper_bound(numax*(1+z));
      i--;
      while (l!=i) {
	nu = l->first/(1+z);
	fn = l->second;
	l++;
	avg += 0.5*(l->first/(1+z) - nu)*(fn + l->second);
      }
      // last part
      i++;
      delta = (numax - l->first/(1+z))/(i->first/(1+z) - l->first/(1+z));
      fn = delta*i->second + (1-delta)*l->second;
      avg += 0.5*(numax - l->first/(1+z))*(fn + i->second);
      return prefactor * avg / (numax - numin);
    }
  }

  // computes the integral of fn/nu
  // in typical range of nu, 1/nu still well-behaved enough
  // to use trapezoid rule (which would be appropriate for fn alone)
  double Filter::computeNorm() const {
    double norm = 0;
    Filter::const_iterator iter=Filter::begin(), last = --Filter::end();
    double nu, fn, nu2, fn2;
    while (iter != last) {
      nu = iter->first / (1+z);
      fn = iter->second;
      iter++;
      nu2 = iter->first / (1+z);
      fn2 = iter->second;
      norm += 0.5*(nu2 - nu)*(fn/nu + fn2/nu2);
    }
    return norm * prefactor;
  }

  // computes integral of fn*log(3e3/nu)/nu
  double Filter::computeLambdaEff() const {
    double t = 0, norm = 0;
    Filter::const_iterator iter=Filter::begin(), last = --Filter::end();
    double nu, fn, nu2, fn2;
    while (iter != last) {
      nu = iter->first / (1+z);
      fn = iter->second;
      iter++;
      nu2 = iter->first / (1+z);
      fn2 = iter->second;
      norm += 0.5*(nu2 - nu)*(fn/nu + fn2/nu2);
      t += 0.5*(nu2 - nu)*(fn/nu*log(3e3/nu) + fn2/nu2*log(3e3/nu2));
    }
    return exp(t/norm);
  }

  // width integrand only mildly nonlinear: use trapezoid rule
  double Filter::computeWidth() const {
    double width = 0;
    double le = computeLambdaEff();
    Filter::const_iterator iter=Filter::begin(), last = --Filter::end();
    double nu, fn, nu2, fn2;
    while (iter != last) {
      nu = iter->first / (1+z); 
      fn = iter->second;
      iter++;
      nu2 = iter->first / (1+z); 
      fn2 = iter->second;
      width += 0.5*(nu2 - nu)*
	(fn*shapelens::pow2(log(3e3/nu/le))/nu +
	 fn2*shapelens::pow2(log(3e3/nu2/le))/nu2);
    }
    return width;
  }

  // remove consecutive zeros, but retain one on either side to provide
  // bounds for interpolation/integration
  void Filter::removeZeros(double threshold) {
    Filter::iterator iter = Filter::begin(), next, last = --Filter::end();
    while (iter != last) {
      if (iter->second <= threshold) {
	next = iter;
	next++;
	while (next != Filter::end()) {
	  if (next->second <= threshold)
	    next++;
	  else {
	    break;
	  }
	}
	if (iter != Filter::begin())
	  iter++;
	if (next != Filter::end())
	  next--;
	Filter::erase(iter, next++);
	iter = next;
      }
      else
	iter++;
    }
  }

} // end namespace
