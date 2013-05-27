// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_TEN
#define DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_TEN

#include <dune/fem/function/common/function.hh>
#include <dune/multiscale/problems/constants.hh>
#include <dune/multiscale/problems/base.hh>

namespace Dune {
namespace Multiscale {
namespace Problem {
/** \addtogroup problem_10 Problem::Ten
 * @{ **/
//! ------------ Elliptic Problem 10 -------------------

// linear elliptic model problem - heterogeneous setting

//! For more further details about the implementation see '../base.hh'
//! For details on the classes, see 'example.hh'

// if the diffusion matrix is symmetric, we can use a CG solver, if not, default to BiCGStab.
#define SYMMETRIC_DIFFUSION_MATRIX

// Note that in the following, 'Imp' abbreviates 'Implementation'
namespace Ten {
// description see below 0.05
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
MSCONSTANTFUNCTION(FirstSource, 1.0)
//! ----------------- End Definition of ' f ' ------------------------


//! ----------------- Definition of ' G ' ------------------------
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
    double diff_coef = 0.0;

    double coefficient
      = ( 1.0
          / (8.0 * M_PI
             * M_PI) ) * ( 1.0 + ( 0.5 * cos( 2.0 * M_PI * (x[0] / constants().epsilon) ) * sin( 2.0 * M_PI * (x[1] / constants().epsilon) ) ) );

    double constant_val = 0.0005;

    double r1 = 0.425;
    double r2 = 0.125;

    // check if x part of the ellipse
    double position = ( pow(x[0] - 0.5, 2.0) / pow(r1, 2.0) ) + ( pow(x[1] - 0.5, 2.0) / pow(r2, 2.0) );

    if (position <= 1.0)
    {
      if (position >= 0.9)
      {
        diff_coef = ( (10.0 * position) - 9.0 ) * coefficient;
        diff_coef += ( 1.0 - ( (10.0 * position) - 9.0 ) ) * constant_val;
      } else {
        diff_coef = constant_val;

        r1 = 0.35;
        r2 = 0.022;

        double new_position = ( pow(x[0] - 0.5, 2.0) / pow(r1, 2.0) ) + ( pow(x[1] - 0.5, 2.0) / pow(r2, 2.0) );

        if (new_position <= 1.0)
          diff_coef *= 100.0;

        r1 = 0.25;
        r2 = 0.022;

        new_position = ( pow(x[0] - 0.5, 2.0) / pow(r1, 2.0) ) + ( pow(x[1] - 0.56, 2.0) / pow(r2, 2.0) );

        if (new_position <= 1.0)
          diff_coef *= 100.0;

        new_position = ( pow(x[0] - 0.5, 2.0) / pow(r1, 2.0) ) + ( pow(x[1] - 0.44, 2.0) / pow(r2, 2.0) );

        if (new_position <= 1.0)
          diff_coef *= 100.0;
      }
    } else {
      diff_coef = coefficient;
    }

    const double stab = 0.0;
    flux[0][0] = diff_coef * gradient[0][0] + stab * gradient[0][1];
    flux[0][1] = diff_coef * gradient[0][1] + stab * gradient[0][0];
  } // diffusiveFlux

  // the jacobian matrix (JA^{\epsilon}) of the diffusion operator A^{\epsilon} at the position "\nabla v" in direction
  // "nabla w", i.e.
  // jacobian diffusiv flux = JA^{\epsilon}(\nabla v) nabla w:

  // jacobianDiffusiveFlux = A^{\epsilon}( x , position_gradient ) direction_gradient
  void jacobianDiffusiveFlux(const DomainType& x,
                             const JacobianRangeType& /*position_gradient*/,
                             const JacobianRangeType& direction_gradient,
                             JacobianRangeType& flux) const {
    double diff_coef = 0.0;

    double coefficient
      = ( 1.0
          / (8.0 * M_PI
             * M_PI) ) * ( 1.0 + ( 0.5 * cos( 2.0 * M_PI * (x[0] / constants().epsilon) ) * sin( 2.0 * M_PI * (x[1] / constants().epsilon) ) ) );

    double constant_val = 0.0005;

    double r1 = 0.425;
    double r2 = 0.125;

    // check if x part of the ellipse
    double position = ( pow(x[0] - 0.5, 2.0) / pow(r1, 2.0) ) + ( pow(x[1] - 0.5, 2.0) / pow(r2, 2.0) );

    if (position <= 1.0)
    {
      if (position >= 0.9)
      {
        diff_coef = ( (10.0 * position) - 9.0 ) * coefficient;
        diff_coef += ( 1.0 - ( (10.0 * position) - 9.0 ) ) * constant_val;
      } else {
        diff_coef = constant_val;

        r1 = 0.35;
        r2 = 0.022;

        double new_position = ( pow(x[0] - 0.5, 2.0) / pow(r1, 2.0) ) + ( pow(x[1] - 0.5, 2.0) / pow(r2, 2.0) );

        if (new_position <= 1.0)
          diff_coef *= 100.0;

        r1 = 0.25;
        r2 = 0.022;

        new_position = ( pow(x[0] - 0.5, 2.0) / pow(r1, 2.0) ) + ( pow(x[1] - 0.56, 2.0) / pow(r2, 2.0) );

        if (new_position <= 1.0)
          diff_coef *= 100.0;

        new_position = ( pow(x[0] - 0.5, 2.0) / pow(r1, 2.0) ) + ( pow(x[1] - 0.44, 2.0) / pow(r2, 2.0) );

        if (new_position <= 1.0)
          diff_coef *= 100.0;
      }
    } else {
      diff_coef = coefficient;
    }
    flux[0][0] = diff_coef * direction_gradient[0][0];
    flux[0][1] = diff_coef * direction_gradient[0][1];
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

} //namespace Ten {
}
} //namespace Multiscale {
} //namespace Dune {

#endif // ifndef DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_TEN