#ifndef DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_TWO
#define DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_TWO

#include <dune/fem/function/common/function.hh>
#include <dune/multiscale/problems/constants.hh>
#include <dune/multiscale/problems/base.hh>

namespace Problem {

/** \addtogroup problem_2 Problem::Two
 * @{ **/
//! ------------ Elliptic Problem 2 -------------------

// linear elliptic model problem - heterogeneous setting
// no exact solution available!

// if the diffusion matrix is symmetric, we can use a CG solver, if not, default to BiCGStab.
#define SYMMETRIC_DIFFUSION_MATRIX

// Note that in the following, 'Imp' abbreviates 'Implementation'

namespace Two {
// default value for epsilon (if not sprecified in the parameter file)
CONSTANTSFUNCTION( 0.05 )

// model problem information
struct ModelProblemData
  : public IModelProblemData
{

  static const bool has_exact_solution = false;

  ModelProblemData()
    : IModelProblemData(constants()){
      assert( constants_.epsilon != 0.0);
      if (constants().get("stochastic_pertubation", false) && !(this->problemAllowsStochastics()) )
         DUNE_THROW(Dune::InvalidStateException, "The problem does not allow stochastic perturbations. Please, switch the key off.");
  }

  //! \copydoc IModelProblemData::getMacroGridFile()
  inline std::string getMacroGridFile() const {
    return("../dune/multiscale/grids/macro_grids/elliptic/earth.dgf");
  }

  // are the coefficients periodic? (e.g. A=A(x/eps))
  // this method is only relevant if you want to use a standard homogenizer
  inline bool problemIsPeriodic() const {
    return false; // = problem is not periodic
  }

  // does the problem allow a stochastic perturbation of the coefficients?
  inline bool problemAllowsStochastics() const {
    return false; // = problem does not allow stochastic perturbations
    // (if you want it, you must add the 'perturb' method provided
    // by 'constants.hh' - see model problems 4 to 7 for examples )
  }

};

//! ----------------- Definition of ' f ' ------------------------
CONSTANTFUNCTION(FirstSource, 1.0)
//! ----------------- End Definition of ' f ' ------------------------


//! ----------------- Definition of ' G ' ------------------------
NULLFUNCTION(SecondSource)
//! ----------------- End Definition of ' G ' ------------------------


//! ----------------- Definition of ' A ' ------------------------
// the linear diffusion operator A^{\epsilon}(x,\xi)=A^{\epsilon}(x) \xi
// A^{\epsilon} : \Omega × R² -> R²
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
  Diffusion(){}
  // in the linear setting, use the structure
  // A^{\epsilon}_i(x,\xi) = A^{\epsilon}_{i1}(x) \xi_1 + A^{\epsilon}_{i2}(x) \xi_2

  // (diffusive) flux = A^{\epsilon}( x , direction )
  // (typically direction is some 'gradient_of_a_function')
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

    flux[0][0] = coefficient * gradient[0][0];
    flux[0][1] = coefficient * gradient[0][1];

  } // diffusiveFlux

  // the jacobian matrix (JA^{\epsilon}) of the diffusion operator A^{\epsilon} at the position "\nabla v" in direction
  // "nabla w", i.e.
  // jacobian diffusiv flux = JA^{\epsilon}(\nabla v) nabla w:

  // jacobianDiffusiveFlux = A^{\epsilon}( x , position_gradient ) direction_gradient
  void jacobianDiffusiveFlux(const DomainType& x,
                             const JacobianRangeType& /*position_gradient*/,
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

    flux[0][0] = coefficient * direction_gradient[0][0];
    flux[0][1] = coefficient * direction_gradient[0][1];

  } // jacobianDiffusiveFlux

  /** \deprecated throws Dune::NotImplemented exception **/
  template < class... Args >
  void evaluate( Args... ) const
  {
    DUNE_THROW(Dune::NotImplemented, "Inadmissible call for 'evaluate'");
  }
};
//! ----------------- End Definition of ' A ' ------------------------


//! ----------------- Definition of ' m ' ----------------------------
CONSTANTFUNCTION(MassTerm,  0.0)
//! ----------------- End Definition of ' m ' ------------------------


//! ----------------- Definition of some dummy -----------------------
NULLFUNCTION(DefaultDummyFunction)
//! ----------------- End Definition of some dummy -------------------


//! ----------------- Definition of ' u ' ----------------------------
//! Exact solution is unknown for this model problem
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

  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType  RangeFieldType;

  typedef DomainFieldType TimeType;
  // essentially: 'DomainFieldType' is the type of an entry of a domain-element.
  // But: it is also used if 'u' (the exact solution) has a time-dependency ('u = u(x,t)').
  // This makes sense since the time-dependency is a one-dimensional element of the 'DomainType' and is therefor also an
  // entry of a domain-element.

public:
  ExactSolution(){}

  // in case 'u' has NO time-dependency use the following method:
  inline void evaluate(const DomainType& /*x*/,
                       RangeType& /*y*/) const {
    DUNE_THROW(Dune::NotImplemented, "Exact solution not available!");
  }

  inline void evaluateJacobian(const DomainType& /*x*/, JacobianRangeType& /*grad_u*/) const {
    DUNE_THROW(Dune::NotImplemented, "Exact solution not available!");
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
//! ----------------- End Definition of ' u ' ------------------------

} //! @} namespace Two {
}
//! @}

#endif // ifndef DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_TWO
