// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_EXAMPLE
#define DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_EXAMPLE

#include <dune/fem/function/common/function.hh>
#include <dune/multiscale/common/traits.hh>
#include <dune/multiscale/problems/base.hh>

#include <math.h>
#include <string>

/**
 * \addtogroup Problem
 *  @{
 * in general we regard problems of the following type:

 * - div ( A (x, ∇u(x) ) ) + m(x) u(x) = f(x) - div G(x)

 *    with a homogenous Dirichlet boundary condition!

 * Here we have:
 * u(x) = exact solution of the problem
 * A(x,·) = diffusion operator (e.g. a scalar function, a matrix or an operator), e.g. with the structure A(x) = A(x,
 x/eps)
 * m(x) = a mass term (or reaction term), e.g. with the structure m(x) = m(x, x/eps )
 * f(x) = first source term with the structure f = f(x) (=> no micro-scale dependency)
 * G(x) = second source term with the structure G = G(x) (=> no micro-scale dependency).
 * (Note that 'G' is directly implemented! We do not implement '- div G'!)

 * A(x,·) can be a monotone operator (=> use HMM, the MsFEM is not implemented for nonlinear problems)


 * ! we use the following class names:


 * class ExactSolution -> describes u(x)
 * methods:
 *   evaluate  u( x )        --> evaluate
 *   evaluate  ∇u( x )       --> jacobian


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
namespace Dune {
namespace Multiscale {
namespace Problem {
namespace Example {

//! model problem information
struct ModelProblemData : public Dune::Multiscale::Problem::IModelProblemData {
  //! is there an exact solution available? true/false
  //! (if 'true' it must be implemented below in the ExactSolution class)
  static const bool has_exact_solution = true;

  ModelProblemData();

  //! \copydoc IModelProblemData::getMacroGridFile();
  std::string getMacroGridFile() const;

  //! are the coefficients periodic? (e.g. A=A(x/eps))
  //! this method is only relevant if you want to use a standard homogenizer
  bool problemIsPeriodic() const;

  //! does the problem allow a stochastic perturbation of the coefficients?
  bool problemAllowsStochastics() const;
};

/**
 FirstSource defines the right hand side (RHS) of the governing problem (i.e. it defines 'f').
 The value of the right hand side (i.e. the value of 'f') at 'x' is accessed by the method 'evaluate'.
 That means 'y := f(x)' and 'y' is returned. It is only important that 'FirstSource' knows the function space
 ('Dune::Multiscale::CommonTraits::FunctionSpaceType') that it
 is part from. (f \in FunctionSpace)
*/
class FirstSource : public Dune::Multiscale::CommonTraits::FunctionBaseType {
public:
  //! evaluate f, i.e. return y=f(x) for a given x
  //! the following method defines 'f':
  void evaluate(const DomainType& x, RangeType& y) const; // end evaluate
};

/** \brief default class for the second source term G.
 * Realization: set G(x) = 0: **/
MSNULLFUNCTION(SecondSource)

//! the (non-linear) diffusion operator A(x,\xi)
//! A : \Omega × R² -> R²
class Diffusion : public DiffusionBase {
public:
  //! in the linear setting, we have the structure
  //! A(x,\xi) = ( A_1(x,\xi), A_2(x,\xi) ) with
  //!              A_i(x,\xi) ) = A_{i1}(x) \xi_1 + A_{i2}(x) \xi_2
  //! (diffusive) flux = A ( x , direction )
  //! (typically direction is some 'gradient_of_a_function')
  void diffusiveFlux(const DomainType& /*x*/, const JacobianRangeType& direction, JacobianRangeType& flux) const;

  /**
      the jacobian matrix (JA) of the diffusion operator A with respect to the direction,
      we evaluate (JA)(x,.) at the position "\nabla v" in direction "nabla w", i.e.
      jacobian diffusiv flux = JA(x,\nabla v) nabla w:
      jacobianDiffusiveFlux = JA( x , position_gradient ) direction_gradient
    */
  void jacobianDiffusiveFlux(const DomainType& /*x*/, const JacobianRangeType& /*position_gradient*/,
                             const JacobianRangeType& direction_gradient, JacobianRangeType& flux) const;
};

// dummmy for a lower order term F( x , u(x) , grad u(x) ) in a PDE like
// - div ( A grad u ) + F ( x , u(x) , grad u(x) ) = f
// NOTE: the operator describing the pde must be a monotone operator
//! ------- Definition of the (possibly nonlinear) lower term F ---------
class LowerOrderTerm
    //  : public Dune::Multiscale::CommonTraits::FunctionBaseType
    //
    {

public:
  void evaluate(const DomainType& /*x*/, RangeType& /*y*/) const {}
  void evaluate(const DomainType& /*x*/, const TimeType& /*time*/, RangeType& /*y*/) const {}
  void evaluate(const DomainType& /*x*/, const RangeType& /*position*/, const JacobianRangeType& /*direction_gradient*/,
                RangeType& /*y*/) const {}
  void position_derivative(const DomainType& /*x*/, const RangeType& /*position*/,
                           const JacobianRangeType& /*direction_gradient*/, RangeType& /*y*/) const {}
  void direction_derivative(const DomainType& /*x*/, const RangeType& /*position*/,
                            const JacobianRangeType& /*direction_gradient*/, JacobianRangeType& /*y*/) const {}
};

//! ----------------- Definition of ' m ' ----------------------------
MSCONSTANTFUNCTION(MassTerm, 0.0)

//! ------------ Definition of homogeneous boundary conditions ----------
MSNULLFUNCTION(DirichletBoundaryCondition)
MSNULLFUNCTION(NeumannBoundaryCondition)

//! ----------------- Definition of some dummy -----------------------
MSNULLFUNCTION(DefaultDummyFunction)

//! ----------------- Definition of ' u ' ----------------------------
//! Exact solution (typically it is unknown)
class ExactSolution : public Dune::Multiscale::CommonTraits::FunctionBaseType {
public:
public:
  //! essentially: 'DomainFieldType' is the type of an entry of a domain-element.
  //! But: it is also used if 'u' (the exact solution) has a time-dependency ('u = u(x,t)').
  //! This makes sense since the time-dependency is a one-dimensional element of the 'typename
  //FunctionSpaceType::DomainType' and is therefor also
  // an
  //! entry of a domain-element.

public:
  //! evaluate 'u(x)'
  void evaluate(const DomainType& x, RangeType& y) const; // evaluate

  //! evaluate '∇u(x)'
  void jacobian(const DomainType& x, JacobianRangeType& grad_u) const {
    grad_u[0][0] = 2.0 * M_PI * cos(2.0 * M_PI * x[0]) * sin(2.0 * M_PI * x[1]);
    grad_u[0][1] = 2.0 * M_PI * sin(2.0 * M_PI * x[0]) * cos(2.0 * M_PI * x[1]);
  } // jacobian

  //! in case 'u' has a time-dependency use the following method:
  //! (some classes might require this as a default implementation)
  void evaluate(const DomainType& x, const TimeType& /*timedummy*/, RangeType& y) const;
};

// seems completely unused
////! ------ Definition of an empty homogenized diffusion matrix -------
// template< class Dune::Multiscale::CommonTraits::FunctionSpaceType, class FieldMatrixImp >
// class HomDiffusion
//  : public Dune::Multiscale::CommonTraits::FunctionBaseType HomDiffusion<
// Dune::Multiscale::CommonTraits::FunctionSpaceType, FieldMatrixImp > >
//{
// public:
//
//  typedef FieldMatrixImp   FieldMatrixType;

// private:
//  typedef HomDiffusion< FunctionSpaceType, FieldMatrixType > ThisType;
//  typedef Dune::Fem::Function< FunctionSpaceType, ThisType > BaseType;

// public:
//  typedef CommonTraits::FunctionSpaceType::typename FunctionSpaceType::DomainType        typename
// FunctionSpaceType::DomainType;
//  typedef CommonTraits::FunctionSpaceType::typename FunctionSpaceType::RangeType         typename
// FunctionSpaceType::RangeType;
//

//
//  typedef CommonTraits::FunctionSpaceType::RangeFieldType  RangeFieldType;

//

// public:
//  const FieldMatrixType& A_hom_;

// public:
//  //! fill the matrix with given values
//  inline explicit HomDiffusion(const FieldMatrixType& A_hom)
//    : A_hom_(A_hom)
//  {}

//  //! in the linear setting, use the structure
//  //! A^0_i(x,\xi) = A^0_{i1}(x) \xi_1 + A^{\epsilon}_{i2}(x) \xi_2
//  //! instantiate all possible cases of the evaluate-method:
//  //! (diffusive) flux = A^0( x , gradient_of_a_function )
//  void diffusiveFlux(const DomainType& /*x*/,
//                     const JacobianRangeType& gradient,
//                     JacobianRangeType& flux) const {
//    if ( constants().get("linear", true) )
//    {
//      flux[0][0] = A_hom_[0][0] * gradient[0][0] + A_hom_[0][1] * gradient[0][1];
//      flux[0][1] = A_hom_[1][0] * gradient[0][0] + A_hom_[1][1] * gradient[0][1];
//    } else {
//      flux[0][0] = A_hom_[0][0] * gradient[0][0] + A_hom_[0][1] * gradient[0][1];
//      flux[0][1] = A_hom_[1][0] * gradient[0][0] + A_hom_[1][1] * gradient[0][1];
//      //! TODO one of the the above is in the wrong branch
//      DUNE_THROW(Dune::NotImplemented,"Nonlinear example not yet implemented.");
//    }
//  } // diffusiveFlux

//  //! the jacobian matrix (JA^0) of the diffusion operator A^0 at the position "\nabla v" in direction
//  //! "nabla w", i.e.
//  //! jacobian diffusiv flux = JA^0(\nabla v) nabla w:
//  //! jacobianDiffusiveFlux = A^0( x , position_gradient ) direction_gradient
//  void jacobianDiffusiveFlux(const DomainType& /*x*/,
//                             const JacobianRangeType& /*position_gradient*/,
//                             const JacobianRangeType& /*direction_gradient*/,
//                             JacobianRangeType& /*flux*/) const {
//    if ( constants().get("linear", true) )
//    {
//      DUNE_THROW(Dune::NotImplemented,"linear example not yet implemented.");
//    } else {
//      DUNE_THROW(Dune::NotImplemented,"Nonlinear example not yet implemented.");
//    }
//  } // jacobianDiffusiveFlux
//};

} // namespace Example {
}
} // namespace Multiscale {
} // namespace Dune {
//! @} End of Doxygen Groups

#endif // ifndef DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_EXAMPLE
