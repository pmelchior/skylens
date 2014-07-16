#ifndef SKYLENS_INTEGRATOR_H
#define SKYLENS_INTEGRATOR_H

#include <gsl/gsl_integration.h>

namespace skylens {
  
  /// Class to integrate member functions of arbitrary classes using GSL 
  /// integration routines. 
  /// Definite, semi-indefinite and indefinite integrations are supported.
  ///
  /// Originally written by Emanuel Ziegler, 
  /// slightly extended and adapted by Matthias Bartlemann.
  template <class T, double (T::*function) (double)>
    class Integrator {
    protected:
      double epsabs, epsrel, error;
      size_t limit;
      int errorcode;
      gsl_function integrand;
      gsl_integration_workspace * work;
      static double wrapper (double x, void *object);
    public:
      
      /// The constructor is initialised with the object itself, absolute and
      /// relative error bounds, and the refinement limit for the integration
      /// interval.
      Integrator (T &object, double abserror = 0.0, double relerror = 1e-4,
		  size_t refinementlimit = 256);
      /// Default destructor.
      ~Integrator ();
      /// Given two boundaries, the definite integral between these boundaries
      /// is carried out.
      double operator () (double min, double max);
      /// With a single boundary, the semi-indefinite integral is carried out
      /// between this boundary and infinity.
      double operator () (double min);
      /// Without boundaries, the indefinite integral is carried out, i.e. the
      /// integral between minus and plus infinity.
      double operator () ();
      /// Returns the integration error.
      inline double getError () const;
      /// Returns the GSL integration error code.
      inline int getErrorCode () const;
  };

  template <class T, double (T::*function) (double)>
  Integrator<T, function>::Integrator
    (T &object, double abserror, double relerror, size_t refinementlimit):
    epsabs (abserror), epsrel (relerror), limit (refinementlimit)
  {
    integrand.function = &wrapper;
    integrand.params = &object;
    work = gsl_integration_workspace_alloc(refinementlimit);
    error = 0.0;
    errorcode = 0;
  }

  template <class T, double (T::*function) (double)>
  Integrator<T,function>::~Integrator ()
  { gsl_integration_workspace_free (work); }

  template <class T, double (T::*function) (double)>
  double Integrator<T,function>::operator () (double min, double max)
  {
    double result;
    errorcode = gsl_integration_qag
      (&integrand, min, max, epsabs, epsrel, limit,
       GSL_INTEG_GAUSS15, work, &result, &error);
    return result;
  }

  template <class T, double (T::*function) (double)>
  double Integrator<T,function>::operator () (double min)
  {
    double result;
    errorcode = gsl_integration_qagiu
      (&integrand, min, epsabs, epsrel, limit, work, &result, &error);
    return result;
  }

  template <class T, double (T::*function) (double)>
  double Integrator<T,function>::operator () ()
  {
    double result;
    errorcode = gsl_integration_qagi
      (&integrand, epsabs, epsrel, limit, work, &result, &error);
    return result;
  }

  template <class T, double (T::*function) (double)>
  double Integrator<T,function>::getError () const
  { return error; }

  template <class T, double (T::*function) (double)>
  int Integrator<T,function>::getErrorCode () const
  { return errorcode; }

  template <class T, double (T::*function) (double)>
  double Integrator<T,function>::wrapper (double x, void *object)
  { return (reinterpret_cast<T *> (object)->*function) (x); };

}

#endif
