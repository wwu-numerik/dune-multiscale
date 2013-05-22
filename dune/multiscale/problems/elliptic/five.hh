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

// default value for epsilon (if not sprecified in the parameter file)
CONSTANTSFUNCTION( 0.05 )

// model problem information
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

  static const int dimDomain = DomainType::dimension;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType  RangeFieldType;

  typedef DomainFieldType TimeType;

public:

  inline void evaluate(const DomainType& x,
                       RangeType& y) const {

    // circle of radius 0.2 around the reentrant corner at (0.5,0.5)
    double distance = sqrt( pow(x[0] - 0.5, 2.0) + pow(x[1] - 0.5, 2.0) );

    if (distance < 0.2)
    { y = 1.0; } else
    { y = 0.1; }
  } // evaluate

  inline void evaluate(const DomainType& x,
                       const TimeType& /*time*/,
                       RangeType& y) const {
    evaluate(x, y);
  }
};

//! ----------------- Definition of ' G ' ----------------------------
MSNULLFUNCTION(SecondSource)

//! ----------------- Definition of ' A ' ------------------------
//! the (non-linear) diffusion operator A^{\epsilon}(x,\xi)
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

  // (diffusive) flux = A^{\epsilon}( x , gradient_of_a_function )
  void diffusiveFlux(const DomainType& x,
                     const JacobianRangeType& gradient,
                     JacobianRangeType& flux) const {

    // coeff.first = ( 0.1 + ( 1.0 * pow(cos( 2.0 * M_PI * (x[0] / epsilon) ), 2.0) ) ) + stochastic perturbation
    // coeff.second = ( 0.1 + 1e-3 + ( 0.1 * sin( 2.0 * M_PI * (x[1] / epsilon) ) ) ) + stochastic perturbation
    const auto coeff = constants().coefficients_variant_A(x);

    flux[0][0] = coeff.first * ( gradient[0][0] + ( (1.0 / 3.0) * pow(gradient[0][0], 3.0) ) );
    flux[0][1] = coeff.second * ( gradient[0][1] + ( (1.0 / 3.0) * pow(gradient[0][1], 3.0) ) );
  } // diffusiveFlux

  // the jacobian matrix (JA^{\epsilon}) of the diffusion operator A^{\epsilon} at the position "\nabla v" in direction
  // "nabla w", i.e.
  // jacobian diffusiv flux = JA^{\epsilon}(\nabla v) nabla w:

  // jacobianDiffusiveFlux = A^{\epsilon}( x , position_gradient ) direction_gradient
  void jacobianDiffusiveFlux(const DomainType& x,
                             const JacobianRangeType& position_gradient,
                             const JacobianRangeType& direction_gradient,
                             JacobianRangeType& flux) const {

    // coeff.first = ( 0.1 + ( 1.0 * pow(cos( 2.0 * M_PI * (x[0] / epsilon) ), 2.0) ) ) + stochastic perturbation
    // coeff.second = ( 0.1 + 1e-3 + ( 0.1 * sin( 2.0 * M_PI * (x[1] / epsilon) ) ) ) + stochastic perturbation
    const auto coeff = constants().coefficients_variant_A(x);

    flux[0][0] = coeff.first * direction_gradient[0][0]
                   * ( 1.0 + pow(position_gradient[0][0], 2.0) );
    flux[0][1] = coeff.second * direction_gradient[0][1]
                   * ( 1.0 + pow(position_gradient[0][1], 2.0) );

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
MSCONSTANTFUNCTION(MassTerm,  0.0)
//! ----------------- End Definition of ' m ' ------------------------


//! ----------------- Definition of some dummy -----------------------
MSNULLFUNCTION(DefaultDummyFunction)
//! ----------------- End Definition of some dummy -------------------


//! ----------------- Definition of ' u ' ----------------------------
//! Exact solution is unknown for this model problem
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

} //! @} namespace Five {
}
} //namespace Multiscale {
} //namespace Dune {

#endif // ifndef DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_FIVE
