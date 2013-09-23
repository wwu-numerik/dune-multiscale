// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_SPE10
#define DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_SPE10

#include <dune/fem/function/common/function.hh>
#include <dune/multiscale/problems/constants.hh>
#include <dune/multiscale/problems/base.hh>

namespace Dune {
namespace Multiscale {
namespace Problem {
/** \addtogroup problem_spe10 Problem::SPE10
 * @{ **/
//! ------------ SPE10 Problem -------------------

// linear elliptic model problem - non-periodic setting

//! For more further details about the implementation see '../base.hh'
//! For details on the classes, see 'example.hh'

namespace SPE10 {

//! model problem information
struct ModelProblemData
  : public IModelProblemData
{
  virtual bool hasExactSolution() const {
    return false;
  }

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

  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

  static const int dimDomain = DomainType::dimension;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType  RangeFieldType;

  typedef DomainFieldType TimeType;

public:
  FirstSource();

  //! evaluate f, i.e. return y=f(x) for a given x
  //! the following method defines 'f':
  void evaluate(const DomainType& x, RangeType& y) const;

  void evaluate(const DomainType& x, const TimeType& /*time*/, RangeType& y) const;

  virtual RangeType evaluate(const DomainType& x) const { return Dune::Multiscale::CommonTraits::FunctionBaseType::evaluate(x); }
};


/** \brief default class for the second source term G.
 * Realization: set G(x) = 0: **/
MSNULLFUNCTION(SecondSource)



//! the linear diffusion operator A^{\epsilon}(x,\xi)=A^{\epsilon}(x) \xi
//! A^{\epsilon} : \Omega × R² -> R²
class Diffusion : public DiffusionBase
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
  Diffusion();

  ~Diffusion();
  //! in the linear setting, use the structure
  //! A^{\epsilon}_i(x,\xi) = A^{\epsilon}_{i1}(x) \xi_1 + A^{\epsilon}_{i2}(x) \xi_2
  //! (diffusive) flux = A^{\epsilon}( x , direction )
  //! (typically direction is some 'gradient_of_a_function')
  void diffusiveFlux(const DomainType& x,
                     const JacobianRangeType& direction,
                     JacobianRangeType& flux) const;

  //! the jacobian matrix (JA^{\epsilon}) of the diffusion operator A^{\epsilon} at the position "\nabla v" in direction
  //! "nabla w", i.e.
  //! jacobian diffusiv flux = JA^{\epsilon}(\nabla v) nabla w:
  //! jacobianDiffusiveFlux = A^{\epsilon}( x , position_gradient ) direction_gradient
  void jacobianDiffusiveFlux(const DomainType& x,
                             const JacobianRangeType& /*position_gradient*/,
                             const JacobianRangeType& direction_gradient,
                             JacobianRangeType& flux) const;

private:
  void readPermeability();


  std::vector<double> deltas_;
  double* permeability_;
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

// we have no exact solution
MSNULLFUNCTION(ExactSolution)

// set zero dirichlet and neumann-values by default
class DirichletData : public ZeroDirichletData {};
class NeumannData : public ZeroNeumannData {};

} //! @} namespace SPE10 {
}
} //namespace Multiscale {
} //namespace Dune {


#endif // ifndef DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_NINE
