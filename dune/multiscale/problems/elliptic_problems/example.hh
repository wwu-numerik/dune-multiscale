#ifndef DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_EXAMPLE
#define DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_EXAMPLE

#include <dune/fem/function/common/function.hh>
#include <dune/multiscale/problems/constants.hh>
#include <dune/multiscale/problems/base.hh>

/**
 * \addtogroup Problem
 *  @{
 * in general we regard problems of the following type:

 * - div ( A (x, ∇u(x) ) ) + m(x) u(x) = f(x) - div G(x)

 *    with a homogenous Dirichlet boundary condition!

 * Here we have:
 * u(x) = exact solution of the problem
 * A(x,·) = diffusion operator (e.g. a scalar function, a matrix or an operator), e.g. with the structure A(x) = A(x, x/eps)
 * m(x) = a mass term (or reaction term), e.g. with the structure m(x) = m(x, x/eps )
 * f(x) = first source term with the structure f = f(x) (=> no micro-scale dependency)
 * G(x) = second source term with the structure G = G(x) (=> no micro-scale dependency).
 * (Note that 'G' is directly implemented! We do not implement '- div G'!)

 * A(x,·) can be a monotone operator (=> use HMM, the MsFEM is not implemented for nonlinear problems)


 * ! we use the following class names:


 * class ExactSolution -> describes u(x)
 * methods:
 *   evaluate  u( x )        --> evaluate
 *   evaluate  ∇u( x )       --> evaluateJacobian


 * class Diffusion -> describes A(x,·)
 * methods:
 *   evaluate A( x , direction )            --> diffusiveFlux
 *   evaluate DA( x , position ) direction  --> jacobianDiffusiveFlux
 *   (DA is the derivative of A with respect to the direction)


 * class MassTerm -> describes m(x)
 * methods:
 *   evaluate m( x )         --> evaluate
 *

 * class FirstSource -> describes f(x)
 * methods:
 *   evaluate f( x )         --> evaluate


 * class SecondSource -> describes G(x)
 * methods:
 *   evaluate G( x )         --> evaluate
 *
 *

 The following example is implemented:
 * Find u with zero boundary condition and
 !! - div ( ( 1 /(8π²) ) ∇u(x) ) = sin( 2 π x_1 ) · sin( 2 π x_2 )

 * In particular:

 * The first source term f(x) is given by:
 *  f(x_1,x_2) :=  sin( 2 π x_1 ) · sin( 2 π x_2 )
 * (since the second source is zero, f will form the right hand side (RHS) of our discrete problem)

 * The second source term G is constantly zero:
 *  G(x) := 0

 * The mass term m is constantly zero:
 *  m(x) := 0

 * The (linear) diffusion operator A(x,d) is given by
 *  A(x,d)=( 1 /(8π²)) d
 *  ( or   A( (x_1, x_2), (d_1 , d_2) ) = ( 1 /(8π²) ) ( d_1 , d_2 )  )
 * i.e. A(x,∇v(x)) = ( 1 /(8π²) ) ∇v(x)
 * (A describes a simple constant matrix)

 * The exact solution u(x) of this problem is given by
 *  u( x_1, x_2 ) = sin( 2 π x_1 ) · sin( 2 π x_2 )

*/
namespace Problem {
namespace Example {
// default value for epsilon (=0.05)
// (in case epsilon is not specified in the paramter file)
CONSTANTSFUNCTION( 0.05 ) // 0.05 is a dummy! no epsilon in our example

// model problem information
struct ModelProblemData
  : public IModelProblemData
{
  // is there an exact solution available? true/false
  // (if 'true' it must be implemented below in the ExactSolution class)
  static const bool has_exact_solution = true;

  ModelProblemData()
    : IModelProblemData(constants()) {
      if (constants().get("stochastic_pertubation", false) && !(this->problemAllowsStochastics()) )
         DUNE_THROW(Dune::InvalidStateException, "The problem does not allow stochastic perturbations. Please, switch the key off.");
  }

  //! \copydoc IModelProblemData::getMacroGridFile()
  // the computational grid in dgf format
  inline std::string getMacroGridFile() const {
    return("../dune/multiscale/grids/macro_grids/elliptic/cube_one.dgf"); // a standard 2d-cube
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


/**
 FirstSource defines the right hand side (RHS) of the governing problem (i.e. it defines 'f').
 The value of the right hand side (i.e. the value of 'f') at 'x' is accessed by the method 'evaluate'.
 That means 'y := f(x)' and 'y' is returned. It is only important that 'FirstSource' knows the function space ('FunctionSpaceImp') that it
 is part from. (f \in FunctionSpace)
*/
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

  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

  static const int dimDomain = DomainType::dimension;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType  RangeFieldType;

  typedef DomainFieldType TimeType;

public:
  //! evaluate f, i.e. return y=f(x) for a given x
  //! the following method defines 'f':
  inline void evaluate(const DomainType& x,
                       RangeType& y) const {

    // f(x_1,x_2) :=  sin( 2 π x_1 ) · sin( 2 π x_2 )
    y = sin(2.0 * M_PI * x[0]) * sin(2.0 * M_PI * x[1]);

  } // end evaluate

};

/** \brief default class for the second source term G.
 * Realization: set G(x) = 0: **/
NULLFUNCTION(SecondSource)

//! the (non-linear) diffusion operator A(x,\xi)
//! A : \Omega × R² -> R²
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

  //! in the linear setting, we have the structure
  //! A(x,\xi) = ( A_1(x,\xi), A_2(x,\xi) ) with
  //!              A_i(x,\xi) ) = A_{i1}(x) \xi_1 + A_{i2}(x) \xi_2
  //! (diffusive) flux = A ( x , direction )
  //! (typically direction is some 'gradient_of_a_function')
  void diffusiveFlux(const DomainType& /*x*/,
                     const JacobianRangeType& direction,
                     JacobianRangeType& flux) const {

    // coefficient = ( 1 /(8π²) )
    double coefficient = 1.0 * ( 1.0 / (8.0 * M_PI * M_PI) );
    //  for 'direction = ∇v(x)', we set 'flux = A(x,∇v(x)) = ( 1 /(8π²) ) ∇v(x)':
    flux[0][0] = coefficient * direction[0][0];
    flux[0][1] = coefficient * direction[0][1];
  } // diffusiveFlux

  /**
      the jacobian matrix (JA) of the diffusion operator A with respect to the direction,
      we evaluate (JA)(x,.) at the position "\nabla v" in direction "nabla w", i.e.
      jacobian diffusiv flux = JA(x,\nabla v) nabla w:
      jacobianDiffusiveFlux = JA( x , position_gradient ) direction_gradient
    */
  void jacobianDiffusiveFlux(const DomainType& /*x*/,
                             const JacobianRangeType& /*position_gradient*/,
                             const JacobianRangeType& direction_gradient,
                             JacobianRangeType& flux) const {
    // Note: in the linear case, we have
    //     'JA( x , position_gradient ) = A( x )'
    // for all directions.

    // coefficient = ( 1 /(8π²) )
    double coefficient = 1.0 * ( 1.0 / (8.0 * M_PI * M_PI) );
    flux[0][0] = coefficient * direction_gradient[0][0];
    flux[0][1] = coefficient * direction_gradient[0][1];
  } // jacobianDiffusiveFlux

  template < class... Args >
  void evaluate( Args... ) const
  {
    DUNE_THROW(Dune::NotImplemented, "Inadmissible call for 'evaluate'");
  }
};

//! ----------------- Definition of ' m ' ----------------------------
CONSTANTFUNCTION(MassTerm,  0.0)

//! ----------------- Definition of some dummy -----------------------
NULLFUNCTION(DefaultDummyFunction)

//! ----------------- Definition of ' u ' ----------------------------
//! Exact solution (typically it is unknown)
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

  //! essentially: 'DomainFieldType' is the type of an entry of a domain-element.
  //! But: it is also used if 'u' (the exact solution) has a time-dependency ('u = u(x,t)').
  //! This makes sense since the time-dependency is a one-dimensional element of the 'DomainType' and is therefor also an
  //! entry of a domain-element.
  typedef DomainFieldType TimeType;

public:
  ExactSolution(){}

  //! evaluate 'u(x)'
  inline void evaluate(const DomainType& x,
                       RangeType& y) const {
    // u( x_1, x_2 ) = sin( 2 π x_1 ) · sin( 2 π x_2 )
    y = sin(2.0 * M_PI * x[0]) * sin(2.0 * M_PI * x[1]);

  } // evaluate

  //! evaluate '∇u(x)'
  inline void evaluateJacobian(const DomainType& x, JacobianRangeType& grad_u) const {
    grad_u[0][0] = 2.0* M_PI* cos(2.0 * M_PI * x[0]) * sin(2.0 * M_PI * x[1]);
    grad_u[0][1] = 2.0* M_PI* sin(2.0 * M_PI * x[0]) * cos(2.0 * M_PI * x[1]);
  } // evaluateJacobian

  //! in case 'u' has a time-dependency use the following method:
  //! (some classes might require this as a default implementation)
  inline void evaluate(const DomainType& x,
                       const TimeType& /*timedummy*/,
                       RangeType& y) const {
    evaluate(x, y);
  }

};

//! ------ Definition of an empty homogenized diffusion matrix -------
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
  const FieldMatrixType& A_hom_;

public:
  //! fill the matrix with given values
  inline explicit HomDiffusion(const FieldMatrixType& A_hom)
    : A_hom_(A_hom)
  {}

  //! in the linear setting, use the structure
  //! A^0_i(x,\xi) = A^0_{i1}(x) \xi_1 + A^{\epsilon}_{i2}(x) \xi_2
  //! instantiate all possible cases of the evaluate-method:
  //! (diffusive) flux = A^0( x , gradient_of_a_function )
  void diffusiveFlux(const DomainType& /*x*/,
                     const JacobianRangeType& gradient,
                     JacobianRangeType& flux) const {
    if ( constants().get("linear", true) )
    {
      flux[0][0] = A_hom_[0][0] * gradient[0][0] + A_hom_[0][1] * gradient[0][1];
      flux[0][1] = A_hom_[1][0] * gradient[0][0] + A_hom_[1][1] * gradient[0][1];
    } else {
      flux[0][0] = A_hom_[0][0] * gradient[0][0] + A_hom_[0][1] * gradient[0][1];
      flux[0][1] = A_hom_[1][0] * gradient[0][0] + A_hom_[1][1] * gradient[0][1];
      //! TODO one of the the above is in the wrong branch
      DUNE_THROW(Dune::NotImplemented,"Nonlinear example not yet implemented.");
    }
  } // diffusiveFlux

  //! the jacobian matrix (JA^0) of the diffusion operator A^0 at the position "\nabla v" in direction
  //! "nabla w", i.e.
  //! jacobian diffusiv flux = JA^0(\nabla v) nabla w:
  //! jacobianDiffusiveFlux = A^0( x , position_gradient ) direction_gradient
  void jacobianDiffusiveFlux(const DomainType& /*x*/,
                             const JacobianRangeType& /*position_gradient*/,
                             const JacobianRangeType& /*direction_gradient*/,
                             JacobianRangeType& /*flux*/) const {
    if ( constants().get("linear", true) )
    {
      DUNE_THROW(Dune::NotImplemented,"linear example not yet implemented.");
    } else {
      DUNE_THROW(Dune::NotImplemented,"Nonlinear example not yet implemented.");
    }
  } // jacobianDiffusiveFlux

  template < class... Args >
  void evaluate( Args... ) const
  {
    DUNE_THROW(Dune::NotImplemented, "Inadmissible call for 'evaluate'");
  }
};


} // namespace Example {
}
//! @} End of Doxygen Groups

#endif // ifndef DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_EXAMPLE
