// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_EIGHT
#define DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_EIGHT

#include <dune/fem/function/common/function.hh>
#include <dune/multiscale/problems/constants.hh>
#include <dune/multiscale/problems/base.hh>
#include <dune/multiscale/common/traits.hh>

namespace Dune {
namespace Multiscale {
namespace Problem {
/** \addtogroup problem_8 Problem::Eight
 * @{ **/
//! ------------ Elliptic Problem 8 -------------------

// a nonlinear model problem - periodic setting

// For details, see 'example.hh'


namespace Eight {
// description see below

// NOTE that (delta/epsilon_est) needs to be a positive integer!

//! model problem information
struct ModelProblemData
  : public IModelProblemData
{
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


//! ----------------- Definition of ' f ' ------------------------


class FirstSource
  : public Dune::Fem::Function< Dune::Multiscale::CommonTraits::FunctionSpaceType,
                                FirstSource >
{
private:
  typedef Dune::Multiscale::CommonTraits::FunctionSpaceType FunctionSpaceType;

public:
  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType  RangeType;

  static const int dimDomain = DomainType::dimension;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType  RangeFieldType;

  typedef DomainFieldType TimeType;

public:
  FirstSource(){}

  void evaluate(const DomainType& x,
                       RangeType& y) const; // evaluate

  void evaluate(const DomainType& x,
                       const TimeType& /*time*/,
                       RangeType& y) const;
};

/** \brief default class for the second source term G.
 * Realization: set G(x) = 0: **/
MSNULLFUNCTION(SecondSource)


//! the (non-linear) diffusion operator A^{\epsilon}(x,\xi)
//! A^{\epsilon} : \Omega × R² -> R²
//!
class Diffusion
  : public Dune::Fem::Function< Dune::Multiscale::CommonTraits::FunctionSpaceType, Diffusion >
{
public:
  typedef Dune::Multiscale::CommonTraits::FunctionSpaceType FunctionSpaceType;

public:
  typedef typename FunctionSpaceType::DomainType        DomainType;
  typedef typename FunctionSpaceType::RangeType         RangeType;
  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType  RangeFieldType;

  typedef DomainFieldType TimeType;

private:

  //! to simplify evaluate
  double additivePart(const DomainType& x, const int i, const int j) const; // additivePart

  // instantiate all possible cases of the evaluate-method:

public:
  Diffusion(){}

  // (diffusive) flux = A^{\epsilon}( x , gradient_of_a_function )
  void diffusiveFlux(const DomainType& x,
                     const JacobianRangeType& gradient,
                     JacobianRangeType& flux) const; // diffusiveFlux

  // the jacobian matrix (JA^{\epsilon}) of the diffusion operator A^{\epsilon} at the position "\nabla v" in direction
  // "nabla w", i.e.
  // jacobian diffusiv flux = JA^{\epsilon}(\nabla v) nabla w:

  // jacobianDiffusiveFlux = A^{\epsilon}( x , position_gradient ) direction_gradient
  void jacobianDiffusiveFlux(const DomainType& x,
                             const JacobianRangeType& position_gradient,
                             const JacobianRangeType& direction_gradient,
                             JacobianRangeType& flux) const; // jacobianDiffusiveFlux
};

// dummmy for a lower order term F( x , u(x) , grad u(x) ) in a PDE like
// - div ( A grad u ) + F ( x , u(x) , grad u(x) ) = f
// NOTE: the operator describing the pde must be a monotone operator
//! ------- Definition of the (possibly nonlinear) lower term F ---------
class LowerOrderTerm
//  : public Dune::Fem::Function< Dune::Multiscale::CommonTraits::FunctionSpaceType,
//                                LowerOrderTerm >
{

public:

  LowerOrderTerm( /*double scaling_factor = 1.0*/ ){} // : scaling_factor_( scaling_factor ) {}

  template< class DomainType , class RangeType >
  void evaluate(const DomainType& x, RangeType& y) const
  {}

  template< class DomainType , class TimeType, class RangeType >
  void evaluate(const DomainType& x, const TimeType& /*time*/, RangeType& y) const
  {}
  
  template< class DomainType , class RangeType, class JacobianRangeType >
  void evaluate(const DomainType& x, const RangeType& position, const JacobianRangeType& direction_gradient, RangeType& y) const
  {}

  template< class DomainType , class RangeType, class JacobianRangeType >
  void position_derivative(const DomainType& x, const RangeType& position, const JacobianRangeType& direction_gradient, RangeType& y) const
  {}

  template< class DomainType , class RangeType, class JacobianRangeType >
  void direction_derivative(const DomainType& x, const RangeType& position, const JacobianRangeType& direction_gradient, JacobianRangeType& y) const
  {}
  
};

//! ----------------- Definition of ' m ' ----------------------------
MSCONSTANTFUNCTION(MassTerm,  0.0)

//! ----------------- Definition of some dummy -----------------------
MSNULLFUNCTION(DefaultDummyFunction)

//! ----------------- Definition of ' u ' ----------------------------
//! Exact solution (typically it is unknown)
class ExactSolution
  : public Dune::Fem::Function< Dune::Multiscale::CommonTraits::FunctionSpaceType, ExactSolution >
{
public:
  typedef Dune::Multiscale::CommonTraits::FunctionSpaceType FunctionSpaceType;

public:
  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType  RangeType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType  RangeFieldType;

  typedef DomainFieldType TimeType;





public:
  ExactSolution(){}

  // in case 'u' has NO time-dependency use the following method:
  void evaluate(const DomainType& x,
                       RangeType& y) const;

  // in case 'u' HAS a time-dependency use the following method:
  // unfortunately GRAPE requires both cases of the method 'evaluate' to be
  // instantiated
  void evaluate(const DomainType& x,
                       const TimeType& /*timedummy*/,
                       RangeType& y) const;

  void jacobian(const DomainType& , typename FunctionSpaceType::JacobianRangeType& ) const;
};

} //! @} namespace Eight {
}
} //namespace Multiscale {
} //namespace Dune {

#endif // ifndef DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_EIGHT
