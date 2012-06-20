#ifndef DUNE_FEM_REDUCEDBASISSPACE_SINESPACE_HH
#define DUNE_FEM_REDUCEDBASISSPACE_SINESPACE_HH

#include "../reducedbasisspace.hh"
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/operator/lagrangeinterpolation.hh>

namespace Dune {
template< class FunctionSpaceImp >
class ExactFunction
  : public Fem::Function< FunctionSpaceImp, ExactFunction< FunctionSpaceImp > >
{
private:
  typedef ExactFunction< FunctionSpaceImp >      ThisType;
  typedef Function< FunctionSpaceImp, ThisType > BaseType;

public:
  typedef FunctionSpaceImp FunctionSpaceType;

public:
  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType  RangeType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType  RangeFieldType;

  static const unsigned int dimDomain = FunctionSpaceType::dimDomain;
  static const unsigned int dimRange = FunctionSpaceType::dimRange;

public:
  inline void evaluate(const DomainType& x, RangeType& y) const {
    y = 1;
    for (unsigned int i = 0; i < dimDomain; ++i)
    {
      const DomainFieldType& xi = x[i];
      y *= xi - xi * xi;
    }
  } // evaluate

  inline void evaluate(const DomainType& x, const RangeFieldType t, RangeType& y) const {
    evaluate(x, y);
  }
};

template< class FunctionSpaceImp >
class SineBaseFunction
  : public Fem::Function< FunctionSpaceImp, SineBaseFunction< FunctionSpaceImp > >
{
private:
  typedef SineBaseFunction< FunctionSpaceImp >   ThisType;
  typedef Function< FunctionSpaceImp, ThisType > BaseType;

public:
  typedef FunctionSpaceImp FunctionSpaceType;

public:
  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType  RangeType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType  RangeFieldType;

  static const unsigned int dimDomain = FunctionSpaceType::dimDomain;
  static const unsigned int dimRange = FunctionSpaceType::dimRange;

  typedef FieldVector< int, dimDomain > CoefficientType;

public:
  explicit SineBaseFunction(const CoefficientType coefficient)
    : coefficient_(coefficient) {
    base_func_counter_ += 1;
    base_func_id_ = base_func_counter_;
    // std :: cout << "base_func_id_ = " << base_func_id_ << std :: endl;
  }

  inline void evaluate(const DomainType& x, RangeType& y) const {
    y = 1;

    for (unsigned int i = 0; i < dimDomain; ++i)
    {
      y *= sqrt(2) * sin(M_PI * coefficient_[i] * x[i]);
    }

    // the fifth basis function - counting starts with zero
    if ( (coefficient_[0] == 2) && (coefficient_[1] == 2) )
    { y = 1; }
  } // evaluate

  inline void evaluate(const DomainType& x, const RangeFieldType t, RangeType& y) const {
    evaluate(x, y);
  }

protected:
  const CoefficientType coefficient_;
  int base_func_id_;

  static int base_func_counter_;
};

template< class FunctionSpaceImp >
int SineBaseFunction< FunctionSpaceImp >::base_func_counter_ = -1;

// number of basis functions = maxCoefficient^worlddim
template< class FEMSpaceImp, unsigned int maxCoefficient >
class SineReducedBasisSpace
  : public Dune::Multiscale::ReducedBasisSpace< AdaptiveDiscreteFunction< FEMSpaceImp > >
{
private:
  typedef SineReducedBasisSpace< FEMSpaceImp,
                                 maxCoefficient > ThisType;
  typedef Dune::Multiscale
    ::ReducedBasisSpace
  < AdaptiveDiscreteFunction
    < FEMSpaceImp > > BaseType;

public:
  typedef FEMSpaceImp FEMSpaceType;

  typedef AdaptiveDiscreteFunction< FEMSpaceType > FEMFunctionType;

  typedef SineBaseFunction< typename FEMSpaceType
                              ::FunctionSpaceType > ContinuousBaseFunctionType;

private:
  typedef typename ContinuousBaseFunctionType::CoefficientType CoefficientType;

public:
  inline explicit SineReducedBasisSpace(FEMSpaceType& FEMSpace)
    : BaseType(FEMSpace) {
    // CoefficientType coefficient( -maxCoefficient );
    CoefficientType coefficient(1.0);

    for (unsigned int i = 0; i < maxCoefficient; ++i)
    {
      for (unsigned int j = 0; j < maxCoefficient; ++j)
      {
        coefficient[0] = i + 1;
        coefficient[1] = j + 1;
        addBaseFunction(coefficient);
        // std :: cout << "coefficient = " << coefficient << std :: endl;
      }
    }
  }

private:
  static inline int abs(const CoefficientType& coefficient) {
    int value = 0;

    for (unsigned int i = 0; i < CoefficientType::dimension; ++i)
      value += (coefficient[i] < 0 ? -coefficient[i] : coefficient[i]);
    return value;
  } // abs

  inline void addBaseFunction(const CoefficientType& coefficient) {
    FEMFunctionType discreteBaseFunction("base function", baseFunctionSpace_);
    ContinuousBaseFunctionType continuousBaseFunction(coefficient);

    LagrangeInterpolation< FEMFunctionType >
      ::interpolateFunction(continuousBaseFunction, discreteBaseFunction);
    BaseType::addBaseFunction(discreteBaseFunction);
  } // addBaseFunction

public:
  #if 1
  inline void getBaseFunction(unsigned int i, FEMFunctionType& base_func) const {
    BaseType::baseFunction(i, base_func);
  }

  #endif // if 1

  using BaseType::baseFunctionSpace_;
};
}

#endif // ifndef DUNE_FEM_REDUCEDBASISSPACE_SINESPACE_HH
