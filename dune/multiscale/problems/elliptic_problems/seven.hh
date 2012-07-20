#ifndef DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_SEVEN
#define DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_SEVEN

#include <dune/fem/function/common/function.hh>
#include <dune/multiscale/problems/constants.hh>
#include <dune/multiscale/problems/base.hh>

// ! USE THIS ONE FOR MSFEM TESTS! PURELY LINEAR ELLIPTIC!

// !############################## Elliptic Problem 6 ###################################


// eps = 0.001 => H Ref = 16
// eps = 0.002 => H Ref = 14
// eps = 0.004 => H Ref = 12
// eps = 0.008 => H Ref = 10
// eps = 0.016 => H Ref = 8
// eps = 0.032 => H Ref = 6
// eps = 0.064 => H Ref = 4
// eps = 0.128 => H Ref = 2
// eps = 0.256 => H Ref = 0

// NOTE that (delta/epsilon_est) needs to be a positive integer!

// is an exact solution available?
// #define EXACTSOLUTION_AVAILABLE

// Note that in the following, 'Imp' abbreviates 'Implementation'
namespace Problem {
namespace Seven {
// description see below 0.05
CONSTANTSFUNCTION(0.01, 0.01, 0.01)

// model problem information
struct ModelProblemData
  : public IModelProblemData
{
  ModelProblemData(const std::string filename = "no_name")
    : IModelProblemData(constants(), filename) {
  }

  inline int get_Number_of_Model_Problem() const {
    return 7;
  }

  //! \copydoc IModelProblemData::getMacroGridFile()
  inline std::string getMacroGridFile() const {
    return("../dune/multiscale/grids/macro_grids/elliptic/cube_three.dgf");
  }

  //! \copydoc IModelProblemData::getRefinementLevelReferenceProblem()
  inline int getRefinementLevelReferenceProblem() const {
    return 10;
  }
};

// !FirstSource defines the right hand side (RHS) of the governing problem (i.e. it defines 'f').
// The value of the right hand side (i.e. the value of 'f') at 'x' is accessed by the method 'evaluate'. That means 'y
// := f(x)' and 'y' is returned. It is only important that 'RHSFunction' knows the function space ('FuncSpace') that it
// is part from. (f \in FunctionSpace)

template< class FunctionSpaceImp >
class FirstSource
  : public Dune::Fem::Function< FunctionSpaceImp, FirstSource< FunctionSpaceImp > >
{
public:
  typedef FunctionSpaceImp FunctionSpaceType;

private:
  typedef FirstSource< FunctionSpaceType >                   ThisType;
  typedef Dune::Fem::Function< FunctionSpaceType, ThisType > BaseType;

public:
  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType  RangeType;

  static const int dimDomain = DomainType::dimension;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType  RangeFieldType;

  typedef DomainFieldType TimeType;

public:
  inline void evaluate(const DomainType& /*x*/,
                       RangeType& y) const {
    y = 1.0;
  }

  inline void evaluate(const DomainType& x,
                       const TimeType& /*time*/,
                       RangeType& y) const {
    evaluate(x, y);
  }
};

  /** \brief default class for the second source term G.
   * Realization: set G(x) = 0: **/
  NULLFUNCTION(SecondSource)

// the (non-linear) diffusion operator A^{\epsilon}(x,\xi)
// A^{\epsilon} : R^d -> R^d

template< class FunctionSpaceImp >
class Diffusion
  : public Dune::Fem::Function< FunctionSpaceImp, Diffusion< FunctionSpaceImp > >
{
public:
  typedef FunctionSpaceImp FunctionSpaceType;

private:
  typedef Diffusion< FunctionSpaceType >                     ThisType;
  typedef Dune::Fem::Function< FunctionSpaceType, ThisType > BaseType;

public:
  typedef typename FunctionSpaceType::DomainType        DomainType;
  typedef typename FunctionSpaceType::RangeType         RangeType;
  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType  RangeFieldType;

  typedef DomainFieldType TimeType;

public:
  // in the linear setting, use the structure
  // A^{\epsilon}_i(x,\xi) = A^{\epsilon}_{i1}(x) \xi_1 + A^{\epsilon}_{i2}(x) \xi_2
  // the usage of an evaluate method with "evaluate ( i, j, x, y, z)" should be avoided
  // use "evaluate ( i, x, y, z)" instead and return RangeType-vector.

  // the following method generates arbitrary numbers, with a log-normal distribution
  // the expected value (the value with the highest probability) is:
  // E = exp( m + (s²/2))
  inline double rand_log_normal(float& m, float& s) const {
    // float m = 0.0;
    // float s = 0.1;

    // m is a real number
    // s is a postiv real number. s is a meassure for the variance. The smaller s, the smaller the variance from the
    // expected value

    // we use the box-muller method to generate the number:

    float x1, x2, w;

    float random_number_1, random_number_2, random_number_3;

    do {
      do {
        random_number_1 = std::rand();
        random_number_2 = std::rand();
      } while ( (random_number_1 == 0.0) || (random_number_2 == 0.0) );

      if (random_number_1 > random_number_2)
      { x1 = ( 2.0 * (random_number_2 / random_number_1) ) - 1.0; } else
      { x1 = ( 2.0 * (random_number_1 / random_number_2) ) - 1.0; }

      random_number_3 = std::rand();

      if (random_number_3 > random_number_2)
      { x2 = ( 2.0 * (random_number_2 / random_number_3) ) - 1.0; } else
      { x2 = ( 2.0 * (random_number_3 / random_number_2) ) - 1.0; }

      w = x1 * x1 + x2 * x2;
    } while (w >= 1.0);

    w = sqrt( ( -2.0 * log(w) ) / w );

    // the log-normal arbitrary number:
    return exp( m + (x1 * w * s) );
  }   // end method

  // instantiate all possible cases of the evaluate-method:

  // (diffusive) flux = A^{\epsilon}( x , gradient_of_a_function )
  void diffusiveFlux(const DomainType& x,
                     const JacobianRangeType& gradient,
                     JacobianRangeType& flux) const {
// in case of a stochastic perturbation:
    #ifdef STOCHASTIC_PERTURBATION

    float m = 0.0;
    float s = VARIANCE;
    // the expected value in case of a log-normal distribution:
    float expected_value = exp( m + (pow(s, 2.0) / 2.0) );
    double arb_num = rand_log_normal(m, s);
    // std :: cout << "arb_num = " << arb_num << std :: endl;

    float perturbation = arb_num - expected_value;

    #endif // ifdef STOCHASTIC_PERTURBATION

    double coefficient_0 = 1.01 + cos( 2.0 * M_PI * (x[0] / constants().epsilon) );    // !x[ 0 ] + x[ 1 ]; //! cos( 2.0 * M_PI *
                                                                           // (x[0] / constants().epsilon) );
    double coefficient_1 = 1.01 + cos( 2.0 * M_PI * (x[0] / constants().epsilon) );    // !x[ 0 ] + x[ 1 ]; //! cos( 2.0 * M_PI *
                                                                           // (x[0] / constants().epsilon) );

    #ifdef STOCHASTIC_PERTURBATION

    coefficient_0 += perturbation;
    coefficient_1 += perturbation;

    if (coefficient_0 < 0.0001)
    { coefficient_0 = 0.0001; }

    if (coefficient_1 < 0.0001)
    { coefficient_1 = 0.0001; }

    #endif // ifdef STOCHASTIC_PERTURBATION

    flux[0][0] = coefficient_0 * gradient[0][0];
    flux[0][1] = coefficient_1 * gradient[0][1];
  } // diffusiveFlux

  // the jacobian matrix (JA^{\epsilon}) of the diffusion operator A^{\epsilon} at the position "\nabla v" in direction
  // "nabla w", i.e.
  // jacobian diffusiv flux = JA^{\epsilon}(\nabla v) nabla w:

  // jacobianDiffusiveFlux = A^{\epsilon}( x , position_gradient ) direction_gradient
  void jacobianDiffusiveFlux(const DomainType& x,
                             const JacobianRangeType& /*position_gradient*/,
                             const JacobianRangeType& direction_gradient,
                             JacobianRangeType& flux) const {
    #ifdef STOCHASTIC_PERTURBATION

    float m = 0.0;
    float s = VARIANCE;
    // the expected value in case of a log-normal distribution:
    float expected_value = exp( m + (pow(s, 2.0) / 2.0) );
    double arb_num = rand_log_normal(m, s);
    // std :: cout << "arb_num = " << arb_num << std :: endl;

    float perturbation = arb_num - expected_value;

    #endif // ifdef STOCHASTIC_PERTURBATION

    double coefficient_0 = 1.01 + cos( 2.0 * M_PI * (x[0] / constants().epsilon) );    // !x[ 0 ] + x[ 1 ]; //! cos( 2.0 * M_PI *
                                                                           // (x[0] / constants().epsilon) );
    double coefficient_1 = 1.01 + cos( 2.0 * M_PI * (x[0] / constants().epsilon) );    // !x[ 0 ] + x[ 1 ]; //!cos( 2.0 * M_PI *
                                                                           // (x[0] / constants().epsilon) );

    #ifdef STOCHASTIC_PERTURBATION

    coefficient_0 += perturbation;
    coefficient_1 += perturbation;

    if (coefficient_0 < 0.0001)
    { coefficient_0 = 0.0001; }

    if (coefficient_1 < 0.0001)
    { coefficient_1 = 0.0001; }

// std :: cout << "coefficient_0 = " << coefficient_0 << std :: endl;
// std :: cout << "coefficient_0 + perturbation = " << coefficient_0 + perturbation << std :: endl;

    #endif // ifdef STOCHASTIC_PERTURBATION

    flux[0][0] = coefficient_0 * direction_gradient[0][0];
    flux[0][1] = coefficient_1 * direction_gradient[0][1];
  } // jacobianDiffusiveFlux

// deprecated
// ! defaults (not to be used):
// z_i = A^{\epsilon}_i(x,vec)
// instantiate all possible cases of the evaluate-method:
  inline void evaluate(const int i,
                       const int j,
                       const DomainType& x,
                       RangeType& z) const {
    if (i == j)
    { z = 1.01 + cos(2.0 * M_PI * x[0]); }       // !x[ 0 ] + x[ 1 ];} //!cos( 2.0 * M_PI * x[0] ); }
    else
    { z = 0.0; }
  }

  inline void evaluate(const int /*i*/,
                       const int /*j*/,
                       const DomainType& /*x*/,
                       const DomainType& /*y*/,
                       RangeType& /*z*/) const {
    DUNE_THROW(Dune::NotImplemented, "Inadmissible call for 'evaluate'");
  } // evaluate

  inline void evaluate(const int /*i*/,
                       const int /*j*/,
                       const DomainType& /*x*/,
                       const TimeType& /*time*/,
                       RangeType& /*z*/) const {
    DUNE_THROW(Dune::NotImplemented, "Inadmissible call for 'evaluate'");
  } // evaluate

  // dummy implementation
  inline void evaluate(const DomainType& /*x*/,
                       RangeType& /*y*/) const {
    DUNE_THROW(Dune::NotImplemented, "Inadmissible call for 'evaluate'");
  } // evaluate
};

template< class FunctionSpaceImp, class FieldMatrixImp >
class HomDiffusion
  : public Dune::Fem::Function< FunctionSpaceImp, HomDiffusion< FunctionSpaceImp, FieldMatrixImp > >
{
public:
  typedef FunctionSpaceImp FunctionSpaceType;
  typedef FieldMatrixImp   FieldMatrixType;

private:
  typedef HomDiffusion< FunctionSpaceType, FieldMatrixType > ThisType;
  typedef Dune::Fem::Function< FunctionSpaceType, ThisType > BaseType;

public:
  typedef typename FunctionSpaceType::DomainType        DomainType;
  typedef typename FunctionSpaceType::RangeType         RangeType;
  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType  RangeFieldType;

  typedef DomainFieldType TimeType;

public:
  FieldMatrixType* A_hom_;

public:
  inline explicit HomDiffusion(FieldMatrixType& A_hom)
    : A_hom_(&A_hom)
  {}

  // in the linear setting, use the structure
  // A^{\epsilon}_i(x,\xi) = A^{\epsilon}_{i1}(x) \xi_1 + A^{\epsilon}_{i2}(x) \xi_2
  // the usage of an evaluate method with "evaluate ( i, j, x, y, z)" should be avoided
  // use "evaluate ( i, x, y, z)" instead and return RangeType-vector.

  // instantiate all possible cases of the evaluate-method:

  // (diffusive) flux = A^{\epsilon}( x , gradient_of_a_function )
  void diffusiveFlux(const DomainType& /*x*/,
                     const JacobianRangeType& gradient,
                     JacobianRangeType& flux) const {
    flux[0][0] = (*A_hom_)[0][0] * gradient[0][0] + (*A_hom_)[0][1] * gradient[0][1];
    flux[0][1] = (*A_hom_)[1][0] * gradient[0][0] + (*A_hom_)[1][1] * gradient[0][1];
  }

  /** the jacobian matrix (JA^{\epsilon}) of the diffusion operator A^{\epsilon} at the position "\nabla v" in direction
   * "nabla w", i.e.
   * jacobian diffusiv flux = JA^{\epsilon}(\nabla v) nabla w:
   * jacobianDiffusiveFlux = A^{\epsilon}( x , position_gradient ) direction_gradient 
  **/
  void jacobianDiffusiveFlux(const DomainType&,
                             const JacobianRangeType&,
                             const JacobianRangeType&,
                             JacobianRangeType&) const {
    DUNE_THROW(Dune::NotImplemented,"");
  } // jacobianDiffusiveFlux

  template < class... Args >
  void evaluate( Args... ) const
  {
    DUNE_THROW(Dune::NotImplemented, "Inadmissible call for 'evaluate'");
  }
};

// define the mass term:
template< class FunctionSpaceImp >
class MassTerm
  : public Dune::Fem::Function< FunctionSpaceImp, MassTerm< FunctionSpaceImp > >
{
public:
  typedef FunctionSpaceImp FunctionSpaceType;

private:
  typedef MassTerm< FunctionSpaceType >                      ThisType;
  typedef Dune::Fem::Function< FunctionSpaceType, ThisType > BaseType;

public:
  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType  RangeType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType  RangeFieldType;

  typedef DomainFieldType TimeType;

public:
  inline void evaluate(const DomainType& /*x*/,
                       RangeType& y) const {
    y[0] = 0.00001;
  }

  // dummy implementation
  inline void evaluate(const DomainType& x,
                       const TimeType /*time*/,
                       RangeType& y) const {
    DSC_LOG_ERROR << "WARNING! Wrong call for 'evaluate' method of the MassTerm class (evaluate(x,t,y)). Return 0.0."
              << std::endl;
    return evaluate(x, y);
  }
};

//! a dummy function class for functions, vectors and matrices
NULLFUNCTION(DefaultDummyFunction)
// ! Exact solution (typically it is unknown)
template< class FunctionSpaceImp >
class ExactSolution
  : public Dune::Fem::Function< FunctionSpaceImp, ExactSolution< FunctionSpaceImp > >
{
public:
  typedef FunctionSpaceImp FunctionSpaceType;

private:
  typedef ExactSolution< FunctionSpaceType >                 ThisType;
  typedef Dune::Fem::Function< FunctionSpaceType, ThisType > BaseType;

public:
  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType  RangeType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType  RangeFieldType;

  typedef DomainFieldType TimeType;
  // essentially: 'DomainFieldType' is the type of an entry of a domain-element.
  // But: it is also used if 'u' (the exact solution) has a time-dependency ('u = u(x,t)').
  // This makes sense since the time-dependency is a one-dimensional element of the 'DomainType' and is therefor also an
  // entry of a domain-element.

public:
  // in case 'u' has NO time-dependency use the following method:
  inline void evaluate(const DomainType& /*x*/,
                       RangeType& y) const {
    y = 0.0;
  }

  // in case 'u' HAS a time-dependency use the following method:
  // unfortunately GRAPE requires both cases of the method 'evaluate' to be
  // instantiated
  inline void evaluate(const DomainType& x,
                       const TimeType& /*timedummy*/,
                       RangeType& y) const {
    evaluate(x, y);
  }
};
} // namespace Seven {
}


#endif // ifndef DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_SEVEN
