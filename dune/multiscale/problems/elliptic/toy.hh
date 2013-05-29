// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_TOY
#define DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_TOY

#include <dune/fem/function/common/function.hh>
#include <dune/multiscale/problems/constants.hh>
#include <dune/multiscale/problems/base.hh>

#include <dune/multiscale/tools/misc/linear-lagrange-interpolation.hh>
namespace Dune {
namespace Multiscale {
namespace Problem {
/** \addtogroup problem_0 Problem::Toy
 * @{ **/
// --------- TOY PROBLEM FOR SIMPLE TESTS --------------
//! ------------ Elliptic Problem 0 --------------------
// an epsilon-indepent linear elliptic toy model problem

// we solve for
//  u(x) = x_1 ( 1 - x_1 ) (1 - x_2 ) x_2
// with
//  a_12(x) = a_21(x) = 0 and a_11(x) = a_22(x) = 1 + (x_1)²
// f is defined so that everything fits together

//! For more further details about the implementation see '../base.hh'
//! For details on the classes, see 'example.hh'

// if the diffusion matrix is symmetric, we can use a CG solver, if not, default to BiCGStab.
#define SYMMETRIC_DIFFUSION_MATRIX

// Note that in the following, 'Imp' abbreviates 'Implementation'
namespace Toy {
// default value for epsilon (not required for this toy problem)
CONSTANTSFUNCTION( 1.0 )

// model problem information
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
  typedef typename FunctionSpaceType::DomainType        DomainType;
  typedef typename FunctionSpaceType::RangeType         RangeType;
  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

  static const int dimDomain = DomainType::dimension;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType  RangeFieldType;

  typedef DomainFieldType TimeType;

public:
  FirstSource(){}

  inline void evaluate(const DomainType& x,
                       RangeType& y) const {
    double a_0_x_0 = 1.0 + pow(x[0], 2.0);
    double a_1_x_1 = 1.0 + pow(x[0], 2.0);

    double grad_a_0_x_0 = 2.0 * x[0];
    double grad_a_1_x_1 = 0.0;

    JacobianRangeType grad_u(0.0);

    grad_u[0][0] = (1.0 - x[0]) * (1.0 - x[1]) * x[1] - x[0] * (1.0 - x[1]) * x[1];
    grad_u[0][1] = x[0] * (1.0 - x[0]) * (1.0 - x[1]) - x[0] * (1.0 - x[0]) * x[1];

    const RangeType d_xx_u = (-2.0) * (1.0 - x[1]) * x[1];
    const RangeType d_yy_u = (-2.0) * (1.0 - x[0]) * x[0];

    y = 0.0;
    y -= grad_a_0_x_0 * grad_u[0][0];
    y -= a_0_x_0 * d_xx_u;

    y -= grad_a_1_x_1 * grad_u[0][1];
    y -= a_1_x_1 * d_yy_u;
  } // evaluate

  inline void evaluate(const DomainType& x,
                       const TimeType& /*time*/,
                       RangeType& y) const {
    evaluate(x, y);
  }
};

//! ----------------- End Definition of ' f ' ------------------------


//! ----------------- Definition of ' G ' ------------------------

  /** \brief default class for the second source term G.
   * Realization: set G(x) = 0: **/
  MSNULLFUNCTION(SecondSource)

//! ----------------- End Definition of ' G ' ------------------------



//! ----------------- Definition of ' A ' ------------------------

// the linear diffusion operator A^{\epsilon}(x,\xi)=A^{\epsilon}(x) \xi
// A^{\epsilon} : \Omega × R² -> R²
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
  Diffusion(){}

  // in the linear setting, use the structure
  // A^{\epsilon}_i(x,\xi) = A^{\epsilon}_{i1}(x) \xi_1 + A^{\epsilon}_{i2}(x) \xi_2

  // (diffusive) flux = A^{\epsilon}( x , direction )
  // (typically direction is some 'gradient_of_a_function')
  void diffusiveFlux(const DomainType& x,
                     const JacobianRangeType& gradient,
                     JacobianRangeType& flux) const {
    double a_0 = 1.0 + pow(x[0], 2.0);

    flux[0][0] = a_0 * gradient[0][0];
    flux[0][1] = a_0 * gradient[0][1];
  }

  /** \deprecated throws Dune::NotImplemented exception **/
  template < class... Args >
  void evaluate( Args... ) const
  {
    DUNE_THROW(Dune::NotImplemented, "Inadmissible call for 'evaluate'");
  }

  template < class... Args >
  void jacobianDiffusiveFlux( Args... ) const {
    DUNE_THROW(Dune::NotImplemented, "Dummy body for all-problem compile");
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
// Exact solution (typically it is unknown)
class ExactSolution
  : public Dune::Fem::Function< Dune::Multiscale::CommonTraits::FunctionSpaceType, ExactSolution >
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
  // essentially: 'DomainFieldType' is the type of an entry of a domain-element.
  // But: it is also used if 'u' (the exact solution) has a time-dependency ('u = u(x,t)').
  // This makes sense since the time-dependency is a one-dimensional element of the 'DomainType' and is therefor also an
  // entry of a domain-element.

public:
  ExactSolution(){}

  // in case 'u' has NO time-dependency use the following method:
  inline void evaluate(const DomainType& x,
                       RangeType& y) const {
    y = x[0] * (1.0 - x[0]) * (1.0 - x[1]) * x[1];
  }

  // in case 'u' HAS a time-dependency use the following method:
  // unfortunately GRAPE requires both cases of the method 'evaluate' to be
  // instantiated
  inline void evaluate(const DomainType& x,
                       const TimeType& /*timedummy*/,
                       RangeType& y) const {
    evaluate(x, y);
  }

  inline void evaluateJacobian(const DomainType& x, JacobianRangeType& grad_u) const {
    grad_u[0][0] = (1.0 - x[0]) * (1.0 - x[1]) * x[1] - x[0] * (1.0 - x[1]) * x[1];
    grad_u[0][1] = x[0] * (1.0 - x[0]) * (1.0 - x[1]) - x[0] * (1.0 - x[0]) * x[1];
  }

};
//! ----------------- End Definition of ' u ' ------------------------

} //! @} namespace Toy {
}
} //namespace Multiscale {
} //namespace Dune {

#endif // ifndef DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_TOY
