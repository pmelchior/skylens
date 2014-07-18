#ifndef SKYLENS_FILTER_H
#define SKYLENS_FILTER_H

#include <map>
#include <string>

namespace skylens {
  
  /// Class for filter curves.
  /// Filters are maps from frequency to transmission efficiency 
  /// used to describe observation bands. 
  /// 
  /// Derived from Matthias Bartelmann's libastro,
  /// substantially modified to work as a std::map by Peter Melchior
  class Filter : public std::map<double, double> {
  public:
    Filter();
    Filter(const std::string& filename, double threshold=1e-3);
    Filter& operator*= (double c);
    Filter& operator/= (double c);
    Filter& operator*= (const Filter& f);
    Filter& operator+= (const Filter& f);
    double operator() (double nu) const;
    void save(const std::string& filename) const;
    void removeZeros(double threshold);
    double getNuMin() const;
    double getNuMax() const;
    double computeNorm() const;
    double avg(double numin, double numax) const;
    double computeLambdaEff() const;
    double computeWidth() const;
  protected:
    double prefactor, z;
  };

} // end namespace

#endif
