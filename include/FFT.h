#ifndef SKYLENS_FFT_H
#define SKYLENS_FFT_H

#include <fftw3.h>
#include <shapelens/Image.h>

using shapelens::data_t;

namespace skylens {

/// Class for storing result of one-dimensional FFT.
/// This class exploits that the Fourier transform of a real vector
/// obeys the Hermiticity condition
/// \f[F(i) = F^*(N-i)\f]
/// for a vector index \f$i\f$ and the length of the real vector \f$N\f$.

class FourierTransform1D : public tmv::Vector<std::complex<data_t> > {
 public:
  /// Constructor.
  FourierTransform1D();
  /// Constructor for real vector of size \f$N\f$.
  FourierTransform1D(unsigned int N);
  /// Access operator.
  /// \b CAUTION: If \f$i>N/2\f$, the result is not conjugated
  /// as required by the Hermiticity condition.
  std::complex<data_t>& operator()(unsigned int i);
  /// const access operator.
  /// \b CAUTION: If \f$i>N/2\f$, the result is not conjugated
  /// as required by the Hermiticity condition.
  const std::complex<data_t>& operator()(unsigned int i) const;
  /// Get wavenumber of index \p i of the transform.
  data_t getWavenumber(int i) const;
  /// Resize transform for real vector of size \f$N\f$.
  void resize(unsigned int N);
  /// Get size of real vector.
  int getRealSize() const;
 private:
  int N;
};

/// Class for storing result of one-dimensional FFT.
/// This class exploits that the Fourier transform of a real matrix
/// obeys the Hermiticity condition
/// \f[F(i,j) = F^*(i,J-j)\f]
/// for vector indices \f$i, j\f$ and the column number  of the real matrix \f$J\f$.
class FourierTransform2D : public tmv::Vector<std::complex<data_t> > {
 public:
  /// Constructor.
  FourierTransform2D();
  /// Constructor for real matrix of size \f$N\times J\f$.
  FourierTransform2D(unsigned int N, unsigned int J);
  /// Access operator.
  /// \b CAUTION: If \f$j>J/2\f$, the result is not conjugated
  /// as required by the Hermiticity condition.
  std::complex<data_t>& operator()(unsigned int i, unsigned int j);
  /// const access operator.
  /// \b CAUTION: If \f$j>J/2\f$, the result is not conjugated
  /// as required by the Hermiticity condition.
  const std::complex<data_t>& operator()(unsigned int i, unsigned int j) const;
  /// Get wavenumber of index \p i of the transform.
  std::complex<data_t> getWavenumber(int i, int j) const;
  /// Get vector index for given matrix indices
  unsigned int getIndex(unsigned int i, unsigned int j) const;
  /// Resize transform for real matrix of size \f$N\times J\f$.
  void resize(unsigned int N, unsigned int J);
  /// Get size of real matrix in \p dimension.
  int getRealSize(bool dimension) const;
  /// Copy operator for base-class.
  FourierTransform2D& operator=(const tmv::Vector<std::complex<data_t> >& v);
 private:
  int N,J;
  data_t wavenumber(int k, bool dimension) const;
};

/// class for one- and two-dimensional Fourier Transforms.
class FFT {
 public:
  /// Transform \f$f(x)\rightarrow F(k)\f$.
  static void transform(const tmv::Vector<data_t>& f, FourierTransform1D& F);
  /// Transform \f$F(k)\rightarrow f(x)\f$.
  static void transform(const  FourierTransform1D& F, tmv::Vector<data_t>& f);
  /// Transform two-dimensional \f$f(\vec{x})\rightarrow F(\vec{k})\f$.
  static void transform(const tmv::Matrix<data_t>& f,  FourierTransform2D& F);
  /// Transform two-dimensional \f$F(\vec{k})\rightarrow f(\vec{x})\f$.
  static void transform(const  FourierTransform2D& F, tmv::Matrix<data_t>& f);
  /// Transform Image \f$f(\vec{x})\rightarrow F(\vec{k})\f$.
  static void transform(const Image<data_t>& f,  FourierTransform2D& F);
  /// Transform Image \f$F(\vec{k})\rightarrow f(\vec{x})\f$.
  static void transform(const  FourierTransform2D& F, Image<data_t>& f);
  /// Convolve \p data with kernel
  /// It is assumed that both \p data and \p kernel have the same sizes.
  static void convolve(Image<data_t>& data, const Image<data_t>& kernel);

  static void convolve(Image<data_t>& data, FourierTransform2D& data_transformed, const Image<data_t>& kernel);
  static void conv_multiply(const FourierTransform2D& f1, const FourierTransform2D& f2, FourierTransform2D& target);
};
} // end namespace
#endif // FFT_H
