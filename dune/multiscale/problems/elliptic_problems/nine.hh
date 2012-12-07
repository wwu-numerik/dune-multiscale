#ifndef DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_NINE
#define DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_NINE

#include <dune/fem/function/common/function.hh>
#include <dune/multiscale/problems/constants.hh>
#include <dune/multiscale/problems/base.hh>

//! ------------ Elliptic Problem 9 -------------------

// linear elliptic model problem - periodic setting

//! For more further details about the implementation see '../base.hh'
//! For details on the classes, see 'example.hh'

// if the diffusion matrix is symmetric, we can use a CG solver, if not, default to BiCGStab.
#define SYMMETRIC_DIFFUSION_MATRIX

// Note that in the following, 'Imp' abbreviates 'Implementation'

namespace Problem {
namespace Nine {
// default value for epsilon (if not sprecified in the parameter file)
CONSTANTSFUNCTION( 0.05 )
// model problem information
struct ModelProblemData
  : public IModelProblemData
{
  static const bool has_exact_solution = true;

  ModelProblemData()
    : IModelProblemData(constants()) {
      if (!constants().get("linear", true))
        DUNE_THROW(Dune::InvalidStateException, "problem nine is entirely linear, but problem.linear was false");
      if (constants().get("stochastic_pertubation", false) && !(this->problemAllowsStochastics()) )
         DUNE_THROW(Dune::InvalidStateException, "The problem does not allow stochastic perturbations. Please, switch the key off.");
  }

  //! \copydoc IModelProblemData::getMacroGridFile()
  inline std::string getMacroGridFile() const {
    return("../dune/multiscale/grids/macro_grids/elliptic/cube_three.dgf");
  }

  // are the coefficients periodic? (e.g. A=A(x/eps))
  // this method is only relevant if you want to use a standard homogenizer
  inline bool problemIsPeriodic() const {
    return true; // = problem is periodic
  }
  
  // does the problem allow a stochastic perturbation of the coefficients?
  inline bool problemAllowsStochastics() const {
    return false; // = problem does not allow stochastic perturbations
    // (if you want it, you must add the 'perturb' method provided
    // by 'constants.hh' - see model problems 4 to 7 for examples )
  }
  
};

//! ----------------- Definition of ' f ' ------------------------

template< class FunctionSpaceImp >
// define the (first) source term 'f'
class FirstSource
  : public Dune::Fem::Function< FunctionSpaceImp, FirstSource< FunctionSpaceImp > >
{
public:
  typedef FunctionSpaceImp FunctionSpaceType;

private:
  typedef FirstSource< FunctionSpaceType >                   ThisType;
  typedef Dune::Fem::Function< FunctionSpaceType, ThisType > BaseType;

public:
  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType  RangeType;

  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

  static const int dimDomain = DomainType::dimension;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType  RangeFieldType;

  typedef DomainFieldType TimeType;

public:
  FirstSource(){}

  // evaluate f, i.e. return y=f(x) for a given x
  // the following method defines 'f':
  inline void evaluate(const DomainType& x,
                       RangeType& y) const {

    double coefficient_0 = 2.0 * ( 1.0 / (8.0 * M_PI * M_PI) ) * ( 1.0 / ( 2.0 + cos( 2.0 * M_PI * (x[0] / constants().epsilon) ) ) );
    double coefficient_1 = ( 1.0 / (8.0 * M_PI * M_PI) ) * ( 1.0 + ( 0.5 * cos( 2.0 * M_PI * (x[0] / constants().epsilon) ) ) );

    double d_x0_coefficient_0
      = pow(2.0 + cos( 2.0 * M_PI * (x[0] / constants().epsilon) ), -2.0) * ( 1.0 / (2.0 * M_PI) ) * (1.0 / constants().epsilon) * sin(
      2.0 * M_PI * (x[0] / constants().epsilon) );

    JacobianRangeType grad_u;
    grad_u[0][0] = 2.0* M_PI* cos(2.0 * M_PI * x[0]) * sin(2.0 * M_PI * x[1]);
    grad_u[0][1] = 2.0* M_PI* sin(2.0 * M_PI * x[0]) * cos(2.0 * M_PI * x[1]);

    grad_u[0][0] += (-1.0) * constants().epsilon * M_PI
                    * ( sin(2.0 * M_PI * x[0]) * sin(2.0 * M_PI * x[1]) * sin( 2.0 * M_PI * (x[0] / constants().epsilon) ) );
    grad_u[0][0] += M_PI * ( cos(2.0 * M_PI * x[0]) * sin(2.0 * M_PI * x[1]) * cos( 2.0 * M_PI * (x[0] / constants().epsilon) ) );

    grad_u[0][1] += constants().epsilon * M_PI
                    * ( cos(2.0 * M_PI * x[0]) * cos(2.0 * M_PI * x[1]) * sin( 2.0 * M_PI * (x[0] / constants().epsilon) ) );

    RangeType d_x0_x0_u(0.0);
    d_x0_x0_u -= 4.0 * pow(M_PI, 2.0) * sin(2.0 * M_PI * x[0]) * sin(2.0 * M_PI * x[1]);
    d_x0_x0_u -= 2.0
                 * pow(M_PI,
                       2.0) * ( constants().epsilon + (1.0 / constants().epsilon) ) * cos(2.0 * M_PI * x[0]) * sin(2.0 * M_PI * x[1]) * sin(
      2.0 * M_PI * (x[0] / constants().epsilon) );
    d_x0_x0_u -= 4.0
                 * pow(M_PI, 2.0) * sin(2.0 * M_PI * x[0]) * sin(2.0 * M_PI * x[1]) * cos( 2.0 * M_PI * (x[0] / constants().epsilon) );

    RangeType d_x1_x1_u(0.0);
    d_x1_x1_u -= 4.0 * pow(M_PI, 2.0) * sin(2.0 * M_PI * x[0]) * sin(2.0 * M_PI * x[1]);
    d_x1_x1_u -= 2.0
                 * pow(M_PI,
                       2.0) * constants().epsilon
                 * cos(2.0 * M_PI * x[0]) * sin(2.0 * M_PI * x[1]) * sin( 2.0 * M_PI * (x[0] / constants().epsilon) );

    y = 0.0;
    y -= d_x0_coefficient_0 * grad_u[0][0];
    y -= coefficient_0 * d_x0_x0_u;
    y -= coefficient_1 * d_x1_x1_u;
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
  NULLFUNCTION(SecondSource)

//! ----------------- End Definition of ' G ' ------------------------



//! ----------------- Definition of ' A ' ------------------------

// the linear diffusion operator A^{\epsilon}(x,\xi)=A^{\epsilon}(x) \xi
// A^{\epsilon} : \Omega × R² -> R²

template< class FunctionSpaceImp >
class Diffusion
  : public Dune::Fem::Function< FunctionSpaceImp, Diffusion< FunctionSpaceImp > >
{
public:
  typedef FunctionSpaceImp FunctionSpaceType;

private:
  typedef Diffusion< FunctionSpaceType >                     ThisType;
  typedef Dune::Fem::Function< FunctionSpaceType, ThisType > BaseType;

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
                     const JacobianRangeType& direction,
                     JacobianRangeType& flux) const {
    double coefficient_0 = 2.0 * ( 1.0 / (8.0 * M_PI * M_PI) ) * ( 1.0 / ( 2.0 + cos( 2.0 * M_PI * (x[0] / constants().epsilon) ) ) );
    double coefficient_1 = ( 1.0 / (8.0 * M_PI * M_PI) ) * ( 1.0 + ( 0.5 * cos( 2.0 * M_PI * (x[0] / constants().epsilon) ) ) );

    double stab = 0.0;

    flux[0][0] = coefficient_0 * direction[0][0] + stab * direction[0][1];
    flux[0][1] = coefficient_1 * direction[0][1] + stab * direction[0][0];
  } // diffusiveFlux

  // the jacobian matrix (JA^{\epsilon}) of the diffusion operator A^{\epsilon} at the position "\nabla v" in direction
  // "nabla w", i.e.
  // jacobian diffusiv flux = JA^{\epsilon}(\nabla v) nabla w:

  // jacobianDiffusiveFlux = A^{\epsilon}( x , position_gradient ) direction_gradient
  void jacobianDiffusiveFlux(const DomainType& x,
                             const JacobianRangeType& /*position_gradient*/,
                             const JacobianRangeType& direction_gradient,
                             JacobianRangeType& flux) const {
    double coefficient_0 = 2.0 * ( 1.0 / (8.0 * M_PI * M_PI) ) * ( 1.0 / ( 2.0 + cos( 2.0 * M_PI * (x[0] / constants().epsilon) ) ) );
    double coefficient_1 = ( 1.0 / (8.0 * M_PI * M_PI) ) * ( 1.0 + ( 0.5 * cos( 2.0 * M_PI * (x[0] / constants().epsilon) ) ) );
    flux[0][0] = coefficient_0 * direction_gradient[0][0];
    flux[0][1] = coefficient_1 * direction_gradient[0][1];
  } // jacobianDiffusiveFlux

  template < class... Args >
  void evaluate( Args... ) const
  {
    DUNE_THROW(Dune::NotImplemented, "Inadmissible call for 'evaluate'");
  }
};

//! ----------------- End Definition of ' A ' ------------------------


//! ----------------- Definition of ' m ' ----------------------------
CONSTANTFUNCTION(MassTerm,  0.0)
//! ----------------- End Definition of ' m ' ------------------------


//! ----------------- Definition of some dummy -----------------------
NULLFUNCTION(DefaultDummyFunction)
//! ----------------- End Definition of some dummy -------------------


//! ----------------- Definition of ' u ' ----------------------------
// Exact solution (typically it is unknown)
template< class FunctionSpaceImp >
class ExactSolution
  : public Dune::Fem::Function< FunctionSpaceImp, ExactSolution< FunctionSpaceImp > >
{
public:
  typedef FunctionSpaceImp FunctionSpaceType;

private:
  typedef ExactSolution< FunctionSpaceType >                 ThisType;
  typedef Dune::Fem::Function< FunctionSpaceType, ThisType > BaseType;

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

  // evaluate 'u(x)'
  inline void evaluate(const DomainType& x,
                       RangeType& y) const {
    // approximation obtained by homogenized solution + first corrector

    // coarse part
    y = sin(2.0 * M_PI * x[0]) * sin(2.0 * M_PI * x[1]);

    // fine part // || u_fine_part ||_L2 = 0.00883883 (for eps = 0.05 )
    y += 0.5 * constants().epsilon * ( cos(2.0 * M_PI * x[0]) * sin(2.0 * M_PI * x[1]) * sin( 2.0 * M_PI * (x[0] / constants().epsilon) ) );
  } // evaluate

  // evaluate 'grad u(x)'
  inline void evaluateJacobian(const DomainType& x, JacobianRangeType& grad_u) const {
    grad_u[0][0] = 2.0* M_PI* cos(2.0 * M_PI * x[0]) * sin(2.0 * M_PI * x[1]);
    grad_u[0][1] = 2.0* M_PI* sin(2.0 * M_PI * x[0]) * cos(2.0 * M_PI * x[1]);

    grad_u[0][0] += (-1.0) * constants().epsilon * M_PI
                    * ( sin(2.0 * M_PI * x[0]) * sin(2.0 * M_PI * x[1]) * sin( 2.0 * M_PI * (x[0] / constants().epsilon) ) );
    grad_u[0][0] += M_PI * ( cos(2.0 * M_PI * x[0]) * sin(2.0 * M_PI * x[1]) * cos( 2.0 * M_PI * (x[0] / constants().epsilon) ) );

    grad_u[0][1] += constants().epsilon * M_PI
                    * ( cos(2.0 * M_PI * x[0]) * cos(2.0 * M_PI * x[1]) * sin( 2.0 * M_PI * (x[0] / constants().epsilon) ) );
  } // evaluateJacobian

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

} // namespace Nine {
}



#endif // ifndef DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_NINE
