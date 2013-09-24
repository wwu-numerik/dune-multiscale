// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_FOURTEEN
#define DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_FOURTEEN

#include <dune/fem/function/common/function.hh>
#include <dune/multiscale/problems/constants.hh>
#include <dune/multiscale/problems/base.hh>

namespace Dune {
namespace Multiscale {
namespace Problem {
/** \addtogroup problem_14 Problem::Fourteen
 * @{ **/
//! ------------ Elliptic Problem 14 -------------------

// linear elliptic model problem

//! For more further details about the implementation see '../base.hh'
//! For details on the classes, see 'example.hh'

namespace Fourteen {

//! model problem information
struct ModelProblemData : public IModelProblemData {
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

class FirstSource : public Dune::Multiscale::CommonTraits::FunctionBaseType {
private:
  typedef Dune::Multiscale::CommonTraits::FunctionSpaceType FunctionSpaceType;

public:
  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType RangeType;

  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

  static const int dimDomain = DomainType::dimension;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

  typedef DomainFieldType TimeType;

public:
  FirstSource() {}

  //! evaluate f, i.e. return y=f(x) for a given x
  //! the following method defines 'f':
  void evaluate(const DomainType& x, RangeType& y) const;

  void evaluate(const DomainType& x, const TimeType& /*time*/, RangeType& y) const;

  virtual RangeType evaluate(const DomainType& x) const {
    return Dune::Multiscale::CommonTraits::FunctionBaseType::evaluate(x);
  }
};

/** \brief default class for the second source term G.
 * Realization: set G(x) = 0: **/
MSNULLFUNCTION(SecondSource)

//! the linear diffusion operator A^{\epsilon}(x,\xi)=A^{\epsilon}(x) \xi
//! A^{\epsilon} : \Omega × R² -> R²
class Diffusion : public DiffusionBase {
public:
  typedef Dune::Multiscale::CommonTraits::FunctionSpaceType FunctionSpaceType;

public:
  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType RangeType;
  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

  typedef DomainFieldType TimeType;

public:
  Diffusion();
  //! in the linear setting, use the structure
  //! A^{\epsilon}_i(x,\xi) = A^{\epsilon}_{i1}(x) \xi_1 + A^{\epsilon}_{i2}(x) \xi_2
  //! (diffusive) flux = A^{\epsilon}( x , direction )
  //! (typically direction is some 'gradient_of_a_function')
  void diffusiveFlux(const DomainType& x, const JacobianRangeType& direction, JacobianRangeType& flux) const;

  //! the jacobian matrix (JA^{\epsilon}) of the diffusion operator A^{\epsilon} at the position "\nabla v" in direction
  //! "nabla w", i.e.
  //! jacobian diffusiv flux = JA^{\epsilon}(\nabla v) nabla w:
  //! jacobianDiffusiveFlux = A^{\epsilon}( x , position_gradient ) direction_gradient
  void jacobianDiffusiveFlux(const DomainType& x, const JacobianRangeType& /*position_gradient*/,
                             const JacobianRangeType& direction_gradient, JacobianRangeType& flux) const;
};

//! ----------------- Definition of ' m ' ----------------------------
MSCONSTANTFUNCTION(MassTerm, 0.0)

//! ----------------- Definition of some dummy -----------------------
MSNULLFUNCTION(DefaultDummyFunction)

// dummmy for a lower order term F( x , u(x) , grad u(x) ) in a PDE like
// - div ( A grad u ) + F ( x , u(x) , grad u(x) ) = f
// NOTE: the operator describing the pde must be a monotone operator
//! ------- Definition of the (possibly nonlinear) lower term F ---------
class LowerOrderTerm : public ZeroLowerOrder {};

//! ----------------- Definition of Dirichlet Boundary Condition ------------------------

class DirichletBoundaryCondition : public Dune::Multiscale::CommonTraits::FunctionBaseType {
private:
  typedef Dune::Multiscale::CommonTraits::FunctionSpaceType FunctionSpaceType;

public:
  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType RangeType;

  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

  static const int dimDomain = DomainType::dimension;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

  typedef DomainFieldType TimeType;

public:
  DirichletBoundaryCondition() {}

  void evaluate(const DomainType& x, RangeType& y) const;

  void evaluate(const DomainType& x, const TimeType& /*time*/, RangeType& y) const;

  virtual RangeType evaluate(const DomainType& x) const {
    return Dune::Multiscale::CommonTraits::FunctionBaseType::evaluate(x);
  }
};

//! ----------------- Definition of Neumann Boundary Condition ------------------------

class NeumannBoundaryCondition : public Dune::Multiscale::CommonTraits::FunctionBaseType {
private:
  typedef Dune::Multiscale::CommonTraits::FunctionSpaceType FunctionSpaceType;

public:
  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType RangeType;

  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

  static const int dimDomain = DomainType::dimension;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

  typedef DomainFieldType TimeType;

public:
  NeumannBoundaryCondition() {}

  void evaluate(const DomainType& x, RangeType& y) const;

  void evaluate(const DomainType& x, const TimeType& /*time*/, RangeType& y) const;

  virtual RangeType evaluate(const DomainType& x) const {
    return Dune::Multiscale::CommonTraits::FunctionBaseType::evaluate(x);
  }
};

// Exact solution is unknown:
//! ----------------- Definition of ' u ' ----------------------------
//! Exact solution is unknown for this model problem
class ExactSolution : public Dune::Multiscale::CommonTraits::FunctionBaseType {
public:
  typedef Dune::Multiscale::CommonTraits::FunctionSpaceType FunctionSpaceType;

public:
  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType RangeType;

  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

  typedef DomainFieldType TimeType;

public:
  ExactSolution() {}

  // in case 'u' has NO time-dependency use the following method:
  void evaluate(const DomainType& /*x*/, RangeType& /*y*/) const;

  void jacobian(const DomainType& /*x*/, JacobianRangeType& /*grad_u*/) const;

  // in case 'u' HAS a time-dependency use the following method:
  // unfortunately GRAPE requires both cases of the method 'evaluate' to be
  // instantiated
  void evaluate(const DomainType& x, const TimeType& /*timedummy*/, RangeType& y) const;

  virtual RangeType evaluate(const DomainType& x) const {
    return Dune::Multiscale::CommonTraits::FunctionBaseType::evaluate(x);
  }
};

// set zero dirichlet and neumann-values by default
class DirichletData : public ZeroDirichletData {};
class NeumannData : public ZeroNeumannData {};

} //! @} namespace Fourteen {
}
} // namespace Multiscale {
} // namespace Dune {

#endif // ifndef DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_FOURTEEN
