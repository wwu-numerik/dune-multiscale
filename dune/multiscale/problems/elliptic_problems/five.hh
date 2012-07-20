#ifndef DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_FIVE
#define DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_FIVE

// ! IDENTICAL TO MODEL PROBLEM 2!!!!

#include <dune/fem/function/common/function.hh>
#include <dune/multiscale/problems/constants.hh>
#include <dune/multiscale/problems/base.hh>

// !############################## Elliptic Problem 1 ###################################


// is an exact solution available?
// #define EXACTSOLUTION_AVAILABLE

// Note that in the following, 'Imp' abbreviates 'Implementation'
namespace Problem {
namespace Five {
// description see below
// vorher: 0.13
CONSTANTSFUNCTION(0.05, 0.15, 0.3)
// NOTE that (delta/epsilon_est) needs to be a positive integer!

// model problem information
struct ModelProblemData
  : public IModelProblemData
{
  ModelProblemData(const std::string filename = "no_name")
    : IModelProblemData(constants(), filename) {
  }

  inline int get_Number_of_Model_Problem() const {
    return 2;
  }

  //! \copydoc IModelProblemData::getMacroGridFile()
  inline std::string getMacroGridFile() const {
    return("../dune/multiscale/grids/macro_grids/elliptic/earth.dgf");
  }

  //! \copydoc IModelProblemData::getRefinementLevelReferenceProblem()
  inline int getRefinementLevelReferenceProblem() const {
    return 18;
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
  inline void evaluate(const DomainType& x,
                       RangeType& y) const {
    #ifdef LINEAR_PROBLEM

    y = 1.0;

    #else // ifdef LINEAR_PROBLEM

    if (x[1] >= 0.1)
    { y = 1.0; } else
    { y = 0.1; }

    #endif // ifdef LINEAR_PROBLEM
  } // evaluate

  inline void evaluate(const DomainType& x,
                       const TimeType& /*time*/,
                       RangeType& y) const {
    evaluate(x, y);
  }
};

// ! default class for the second source term G.

// Realization: set G(x) = 0:
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

  // instantiate all possible cases of the evaluate-method:

  // (diffusive) flux = A^{\epsilon}( x , gradient_of_a_function )
  void diffusiveFlux(const DomainType& x,
                     const JacobianRangeType& gradient,
                     JacobianRangeType& flux) const {
    double coefficient = 1.0 + (9.0 / 10.0) * sin(2.0 * M_PI * sqrt( fabs(2.0 * x[0]) ) / constants().epsilon) * sin(
      2.0 * M_PI * pow(1.5 * x[1], 2.0) / constants().epsilon);

    if (x[1] <= 0.3)
    {
      coefficient *= 4.0;
    }

    if ( (x[1] > 0.3) && (x[1] < 0.6) )
    {
      coefficient *= 2.0 * ( ( (-5.0 / 3.0) * x[1] ) + (3.0 / 2.0) );
    }

    if (x[1] >= 0.6)
    {
      coefficient *= 1.0;
    }

    #ifdef LINEAR_PROBLEM

    flux[0][0] = coefficient * gradient[0][0];
    flux[0][1] = coefficient * gradient[0][1];

    #else // ifdef LINEAR_PROBLEM

    flux[0][0] = coefficient * ( gradient[0][0] + ( (1.0 / 3.0) * pow(gradient[0][0], 3.0) ) );
    flux[0][1] = coefficient * ( gradient[0][1] + ( (1.0 / 3.0) * pow(gradient[0][1], 3.0) ) );

    #endif // ifdef LINEAR_PROBLEM
  } // diffusiveFlux

  // the jacobian matrix (JA^{\epsilon}) of the diffusion operator A^{\epsilon} at the position "\nabla v" in direction
  // "nabla w", i.e.
  // jacobian diffusiv flux = JA^{\epsilon}(\nabla v) nabla w:

  // jacobianDiffusiveFlux = A^{\epsilon}( x , position_gradient ) direction_gradient
  void jacobianDiffusiveFlux(const DomainType& x,
                             const JacobianRangeType& position_gradient,
                             const JacobianRangeType& direction_gradient,
                             JacobianRangeType& flux) const {
    double coefficient = 1.0 + (9.0 / 10.0) * sin(2.0 * M_PI * sqrt( fabs(2.0 * x[0]) ) / constants().epsilon) * sin(
      2.0 * M_PI * pow(1.5 * x[1], 2.0) / constants().epsilon);

    if (x[1] <= 0.3)
    {
      coefficient *= 4.0;
    }

    if ( (x[1] > 0.3) && (x[1] < 0.6) )
    {
      coefficient *= 2.0 * ( ( (-5.0 / 3.0) * x[1] ) + (3.0 / 2.0) );
    }

    if (x[1] >= 0.6)
    {
      coefficient *= 1.0;
    }

    #ifdef LINEAR_PROBLEM
    flux[0][0] = coefficient * direction_gradient[0][0];
    flux[0][1] = coefficient * direction_gradient[0][1];
    #else // ifdef LINEAR_PROBLEM
    flux[0][0] = coefficient * direction_gradient[0][0]
                 * ( 1.0 + pow(position_gradient[0][0], 2.0) );
    flux[0][1] = coefficient * direction_gradient[0][1]
                 * ( 1.0 + pow(position_gradient[0][1], 2.0) );
    #endif // ifdef LINEAR_PROBLEM
  } // jacobianDiffusiveFlux

  template < class... Args >
  void evaluate( Args... ) const
  {
    DUNE_THROW(Dune::NotImplemented, "Inadmissible call for 'evaluate'");
  }
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
    DUNE_THROW(Dune::InvalidStateException, "No homogenization available");

    #ifdef LINEAR_PROBLEM
    flux[0][0] = (*A_hom_)[0][0] * gradient[0][0] + (*A_hom_)[0][1] * gradient[0][1];
    flux[0][1] = (*A_hom_)[1][0] * gradient[0][0] + (*A_hom_)[1][1] * gradient[0][1];
    #else // ifdef LINEAR_PROBLEM
    flux[0][0] = (*A_hom_)[0][0] * gradient[0][0] + (*A_hom_)[0][1] * gradient[0][1];
    flux[0][1] = (*A_hom_)[1][0] * gradient[0][0] + (*A_hom_)[1][1] * gradient[0][1];
    //! TODO one of the the above is in the wrong branch
    DUNE_THROW(Dune::NotImplemented,"Nonlinear example not yet implemented.");
    #endif // ifdef LINEAR_PROBLEM
  } // diffusiveFlux

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
    std::cout << "WARNING! Wrong call for 'evaluate' method of the MassTerm class (evaluate(x,t,y)). Return 0.0."
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
  inline void evaluate(const DomainType& x,
                       RangeType& y) const {
    // NOT THE EXACT SOLUTION!!!:
    //!TODO why is there an imp below then?
    y = 0.0;
    DUNE_THROW(Dune::NotImplemented, "Exact solution not available");

    double coefficient = 1.0 + (9.0 / 10.0) * sin(2.0 * M_PI * sqrt( fabs(2.0 * x[0]) ) / constants().epsilon) * sin(
      2.0 * M_PI * pow(1.5 * x[1], 2.0) / constants().epsilon);

    if (x[1] <= 0.3)
    {
      coefficient *= 4.0;
    }

    if ( (x[1] > 0.3) && (x[1] < 0.6) )
    {
      coefficient *= 2.0 * ( ( (-5.0 / 3.0) * x[1] ) + (3.0 / 2.0) );
    }

    if (x[1] >= 0.6)
    {
      coefficient *= 1.0;
    }

    y = coefficient;
  } // evaluate

  // in case 'u' HAS a time-dependency use the following method:
  // unfortunately GRAPE requires both cases of the method 'evaluate' to be
  // instantiated
  inline void evaluate(const DomainType& x,
                       const TimeType& /*timedummy*/,
                       RangeType& y) const {
    evaluate(x, y);
  }
};
} // namespace Five {
}

#endif // ifndef DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_FIVE
