// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_FIVE
#define DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_FIVE

#include <dune/fem/function/common/function.hh>
#include <dune/multiscale/problems/constants.hh>
#include <dune/multiscale/problems/base.hh>
#include <dune/multiscale/common/traits.hh>

namespace Dune {
namespace Multiscale {
namespace Problem {
/** \addtogroup problem_5 Problem::Five
 * @{ **/
//! ------------ Elliptic Problem 5 -------------------

// nonlinear elliptic model problem - periodic setting
// no exact solution available!

//! For more further details about the implementation see '../base.hh'
//! For details on the classes, see 'example.hh'

namespace Five {

//! model problem information
struct ModelProblemData
  : public IModelProblemData
{
  static const bool has_exact_solution = false;

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
  : public Dune::Multiscale::CommonTraits::FunctionBaseType

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

  void evaluate(const DomainType& x,
                       RangeType& y) const; // evaluate

  void evaluate(const DomainType& x,
                       const TimeType& /*time*/,
                       RangeType& y) const;
};

//! ----------------- Definition of ' G ' ----------------------------
MSNULLFUNCTION(SecondSource)

//! ----------------- Definition of ' A ' ------------------------
//! the (non-linear) diffusion operator A^{\epsilon}(x,\xi)
//! A^{\epsilon} : \Omega × R² -> R²
class Diffusion: public DiffusionBase
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

public:

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
                             JacobianRangeType& flux) const;
};


// dummmy for a lower order term F( x , u(x) , grad u(x) ) in a PDE like
// - div ( A grad u ) + F ( x , u(x) , grad u(x) ) = f
// NOTE: the operator describing the pde must be a monotone operator
//! ------- Definition of the (possibly nonlinear) lower term F ---------
class LowerOrderTerm : public ZeroLowerOrder {};

//! ----------------- Definition of ' m ' ----------------------------
MSCONSTANTFUNCTION(MassTerm,  0.0)


//! ------------ Definition of homogeneous boundary conditions ----------
MSNULLFUNCTION(DirichletBoundaryCondition)
MSNULLFUNCTION(NeumannBoundaryCondition)


//! ----------------- Definition of some dummy -----------------------
MSNULLFUNCTION(DefaultDummyFunction)



//! ----------------- Definition of ' u ' ----------------------------
//! Exact solution is unknown for this model problem
class ExactSolution
  : public Dune::Multiscale::CommonTraits::FunctionBaseType
{
public:
  typedef Dune::Multiscale::CommonTraits::FunctionSpaceType FunctionSpaceType;

public:
  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType  RangeType;

  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType  RangeFieldType;

  typedef DomainFieldType TimeType;





public:
  ExactSolution(){}

  // in case 'u' has NO time-dependency use the following method:
  void evaluate(const DomainType& /*x*/,
                       RangeType& /*y*/) const;

  void jacobian(const DomainType& /*x*/, JacobianRangeType& /*grad_u*/) const;

  // in case 'u' HAS a time-dependency use the following method:
  // unfortunately GRAPE requires both cases of the method 'evaluate' to be
  // instantiated
  void evaluate(const DomainType& x,
                       const TimeType& /*timedummy*/,
                       RangeType& y) const;
};

// set zero dirichlet and neumann-values by default
class DirichletData : public ZeroDirichletData {};
class NeumannData : public ZeroNeumannData {};

} //! @} namespace Five {
}
} //namespace Multiscale {
} //namespace Dune {

#endif // ifndef DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_FIVE
