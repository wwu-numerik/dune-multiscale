// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_ELEVEN
#define DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_ELEVEN

#include <dune/fem/function/common/function.hh>
#include <dune/multiscale/problems/constants.hh>
#include <dune/multiscale/problems/base.hh>

namespace Dune {
namespace Multiscale {
namespace Problem {
/** \addtogroup problem_11 Problem::Eleven
 * @{ **/
//! ------------ Elliptic Problem 11 -------------------

// linear elliptic model problem - periodic setting

//! For more further details about the implementation see '../base.hh'
//! For details on the classes, see 'example.hh'

// if the diffusion matrix is symmetric, we can use a CG solver, if not, default to BiCGStab.
#define SYMMETRIC_DIFFUSION_MATRIX


namespace Eleven {

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
  : public Dune::Fem::Function< Dune::Multiscale::CommonTraits::FunctionSpaceType,
                                FirstSource >
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
  FirstSource( ) {}

  //! evaluate f, i.e. return y=f(x) for a given x
  //! the following method defines 'f':
  void evaluate(const DomainType& x, RangeType& y) const;

  void evaluate(const DomainType& x, const TimeType& /*time*/, RangeType& y) const;
};


/** \brief default class for the second source term G.
 * Realization: set G(x) = 0: **/
MSNULLFUNCTION(SecondSource)



//! the linear diffusion operator A^{\epsilon}(x,\xi)=A^{\epsilon}(x) \xi
//! A^{\epsilon} : \Omega × R² -> R²
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

public:
  Diffusion();
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
};

//! ----------------- Definition of ' m ' ----------------------------
MSCONSTANTFUNCTION(MassTerm,  0.0)

//! ----------------- Definition of some dummy -----------------------
MSNULLFUNCTION(DefaultDummyFunction)



// F( x , u(x) , grad u(x) ) 
// NOTE: the operator describing the pde must be a montone operator
//! ------- Definition of the (possibly nonlinear) lower term F ---------
class LowerOrderTerm
  : public Dune::Fem::Function< Dune::Multiscale::CommonTraits::FunctionSpaceType,
                                LowerOrderTerm >
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

  LowerOrderTerm( /*double scaling_factor = 1.0*/ ){} // : scaling_factor_( scaling_factor ) {}

  // dummy
  void evaluate(const DomainType& x, RangeType& y) const;

  // dummy
  void evaluate(const DomainType& x, const TimeType& /*time*/, RangeType& y) const;
  
  // evaluate y = F(x, u(x), \grad u(x))
  // 'position = u(x)', 'direction_gradient = \grad u(x)'
  void evaluate(const DomainType& x, const RangeType& position, const JacobianRangeType& direction_gradient, RangeType& y) const;

  // evaluate position derivative y = d_1 F (x, u(x), \grad u(x))  (derivative with respect to the second componenent 'u(x)')
  // 'position = u(x)', 'direction_gradient = \grad u(x)'
  void position_derivative(const DomainType& x, const RangeType& position, const JacobianRangeType& direction_gradient, RangeType& y) const;

  // evaluate position derivative y = d_2 F (x, u(x), \grad u(x))  (derivative with respect to the third componenent 'grad u(x)')
  // 'position = u(x)', 'direction_gradient = \grad u(x)'
  void direction_derivative(const DomainType& x, const RangeType& position, const JacobianRangeType& direction_gradient, JacobianRangeType& y) const;
  
  //double scaling_factor_;
};


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

  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType  RangeFieldType;

  typedef DomainFieldType TimeType;

public:
  ExactSolution();

  // evaluate 'u(x)'
  void evaluate(const DomainType& x, RangeType& y) const;

  // evaluate 'grad u(x)'
  void jacobian(const DomainType& x, typename FunctionSpaceType::JacobianRangeType& grad_u) const;

  // in case 'u' HAS a time-dependency use the following method:
  // unfortunately GRAPE requires both cases of the method 'evaluate' to be
  // instantiated
  void evaluate(const DomainType& x, const TimeType&, RangeType& y) const;
};

} //! @} namespace Eleven {
}
} //namespace Multiscale {
} //namespace Dune {


#endif // ifndef DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_ELEVEN
