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

// Note that in the following, 'Imp' abbreviates 'Implementation'
namespace Eight {
// description see below

// default value for epsilon
CONSTANTSFUNCTION( 0.001 )
// NOTE that (delta/epsilon_est) needs to be a positive integer!

// model problem information
struct ModelProblemData
  : public IModelProblemData
{
  static const bool has_exact_solution = true;

  ModelProblemData();

  //! \copydoc IModelProblemData::getMacroGridFile();
  inline std::string getMacroGridFile() const;

  //! are the coefficients periodic? (e.g. A=A(x/eps))
  //! this method is only relevant if you want to use a standard homogenizer
  inline bool problemIsPeriodic() const;

  //! does the problem allow a stochastic perturbation of the coefficients?
  inline bool problemAllowsStochastics() const;
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
    y = 0.0;
    y += 2.0 * ( x[0] + x[1] - pow(x[0], 2.0) - pow(x[1], 2.0) );
    y -= 12.0 * pow( (2.0 * x[0]) - 1.0, 2.0 ) * pow( (x[1] * x[1]) - x[1], 3.0 );
    y -= 12.0 * pow( (2.0 * x[1]) - 1.0, 2.0 ) * pow( (x[0] * x[0]) - x[0], 3.0 );
  } // evaluate

  inline void evaluate(const DomainType& x,
                       const TimeType& /*time*/,
                       RangeType& y) const {
    evaluate(x, y);
  }
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
  double additivePart(const DomainType& x, const int i, const int j) const {
    double y = 0.0;

    y -= (x[0] + x[1]) * cos(2.0 * M_PI * x[i] / constants().epsilon) * sin(2.0 * M_PI * x[j] / constants().epsilon);

    double helper1 = 1.0;
    helper1 *= ( 2.0 + sin(2.0 * M_PI * (x[0] + x[1]) / constants().epsilon) );

    double helper2 = 1.0;
    helper2 *= 3.0;
    helper2 *= ( (2.0 * x[i] - 1.0) * (x[j] * x[j] - x[j]) )
               + ( (x[0] + x[1]) * cos(2.0 * M_PI * x[i] / constants().epsilon) * sin(2.0 * M_PI * x[j] / constants().epsilon) );
    helper2 *= (2.0 * x[i] - 1.0) * (x[j] * x[j] - x[j]) * (x[0] + x[1]) * cos(2.0 * M_PI * x[i] / constants().epsilon) * sin(
      2.0 * M_PI * x[j] / constants().epsilon);
    helper2 += pow( (x[0] + x[1]) * cos(2.0 * M_PI * x[i] / constants().epsilon) * sin(2.0 * M_PI * x[j] / constants().epsilon), 3.0 );

    helper1 *= helper2;
    y -= helper1;
    y -= sin(2.0 * M_PI * (x[0] + x[1]) / constants().epsilon) * pow( (2.0 * x[i]) - 1.0, 3.0 ) * pow( (x[j] * x[j]) - x[j], 3.0 );
    return y;
  } // additivePart

  // instantiate all possible cases of the evaluate-method:

public:
  // (diffusive) flux = A^{\epsilon}( x , gradient_of_a_function )
  void diffusiveFlux(const DomainType& x,
                     const JacobianRangeType& gradient,
                     JacobianRangeType& flux) const {
    double coefficient = 2.0 + sin(2.0 * M_PI * (x[0] + x[1]) / constants().epsilon);

    flux[0][0] = gradient[0][0] + ( coefficient * pow(gradient[0][0], 3.0) );
    flux[0][0] -= additivePart(x, 0, 1);
    flux[0][0] *= (-1.0);

    flux[0][1] = gradient[0][1] + ( coefficient * pow(gradient[0][1], 3.0) );
    flux[0][1] -= additivePart(x, 1, 0);
    flux[0][1] *= (-1.0);
  } // diffusiveFlux

  // the jacobian matrix (JA^{\epsilon}) of the diffusion operator A^{\epsilon} at the position "\nabla v" in direction
  // "nabla w", i.e.
  // jacobian diffusiv flux = JA^{\epsilon}(\nabla v) nabla w:

  // jacobianDiffusiveFlux = A^{\epsilon}( x , position_gradient ) direction_gradient
  void jacobianDiffusiveFlux(const DomainType& x,
                             const JacobianRangeType& position_gradient,
                             const JacobianRangeType& direction_gradient,
                             JacobianRangeType& flux) const {
    double coefficient = 2.0 + sin(2.0 * M_PI * (x[0] + x[1]) / constants().epsilon);

    flux[0][0] = direction_gradient[0][0]
                 * ( 1.0 + 3.0 * coefficient * pow(position_gradient[0][0], 2.0) );
    flux[0][1] = direction_gradient[0][1]
                 * ( 1.0 + 3.0 * coefficient * pow(position_gradient[0][1], 2.0) );

    flux[0][0] *= (-1.0);
    flux[0][1] *= (-1.0);
  } // jacobianDiffusiveFlux

  template < class... Args >
  void evaluate( Args... ) const
  {
    DUNE_THROW(Dune::NotImplemented, "Inadmissible call for 'evaluate' method of the Diffusion class! See 'problem_specification.hh' for details.");
  }
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
  // essentially: 'DomainFieldType' is the type of an entry of a domain-element.
  // But: it is also used if 'u' (the exact solution) has a time-dependency ('u = u(x,t)').
  // This makes sense since the time-dependency is a one-dimensional element of the 'DomainType' and is therefor also an
  // entry of a domain-element.

public:

  // in case 'u' has NO time-dependency use the following method:
  inline void evaluate(const DomainType& x,
                       RangeType& y) const {
    y = (-1.0) * ( (x[0] * x[0]) - x[0] ) * ( (x[1] * x[1]) - x[1] );
    y -= constants().epsilon * (x[0] + x[1]) * sin(2.0 * M_PI * x[0] / constants().epsilon) * sin(2.0 * M_PI * x[1] / constants().epsilon);
  }

  // in case 'u' HAS a time-dependency use the following method:
  // unfortunately GRAPE requires both cases of the method 'evaluate' to be
  // instantiated
  inline void evaluate(const DomainType& x,
                       const TimeType& /*timedummy*/,
                       RangeType& y) const {
    evaluate(x, y);
  }

  inline void evaluateJacobian(const DomainType& , typename FunctionSpaceType::JacobianRangeType& ) const {
    DUNE_THROW(Dune::NotImplemented, "Dummy body for all-problem compile");
  }
};

} //! @} namespace Eight {
}
} //namespace Multiscale {
} //namespace Dune {

#endif // ifndef DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_EIGHT
