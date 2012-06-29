#ifndef DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH
#define DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH

#include <dune/fem/function/common/function.hh>
#include <dune/multiscale/problems/constants.hh>
#include <dune/multiscale/problems/base.hh>

// ! is an exact solution available?
// this information should be provided by the 'problem specification file'
// there we define or don't define the macro EXACTSOLUTION_AVAILABLE
#define EXACTSOLUTION_AVAILABLE

// if the diffusion matrix is symmetric, we can use a CG solver, if not, default to BiCGStab.
#define SYMMETRIC_DIFFUSION_MATRIX

// ! USE THIS ONE FOR MSFEM TESTS! PURELY LINEAR ELLIPTIC!

// ! For more further details about the implementation of the following classes, see the end of the file

// in general we regard problems of the following type:

// - div ( A^{\epsilon} (x,\nabla u^{\epsilon}) ) + m^{\epsilon} u^{\epsilon} = f - div G

// Here we have:
// u^{\epsilon} = exact solution of the problem
// A^{\epsilon} = diffusion matrix (or tensor) with the structure A^{\epsilon}(x) = A(x,\frac{x}{\epsilon})
// m^{\epsilon} = a mass term (or reaction term) with the structure m^{\epsilon}(x) = m(x,\frac{x}{\epsilon})
// f = first source term with the structure f = f(x) (=> no micro-scale dependency)
// G = second source term with the structure G = G(x) (=> no micro-scale dependency).
// Note that 'G' is directly implemented! ( we do not implement '- div G'!)

// Note, that A^{\epsilon} is a monotone operator

// !############################## Elliptic Problem 9 ###################################

// ! we define:

// The entries of the operator A^{\epsilon} by
// ! a^{\epsilon}_{1}(x_1,x_2) := (**y_1**,**y_2**)
// ! a^{\epsilon}_{2}(x_1,x_2) := (**y_1**,**y_2**)

// The mass (or reaction) term m^{\epsilon} is given by:
// ! m^{\epsilon} := \epsilon
// Since \epsilon tends to zero, we may say that we do not have a real mass term for our
// problem. It is a simple condition to fix the solution which is only unique up to a constant.
// In fact we still approximate the solution of the problem without mass.

// The first source term f is given by:
// ! f(x) := ****
// since the second source is zero, f will form the right hand side (RHS) of our discrete problem

// The second source term G is constantly zero:
// ! G(x) := 0

// !FirstSource defines the right hand side (RHS) of the governing problem (i.e. it defines 'f').
// The value of the right hand side (i.e. the value of 'f') at 'x' is accessed by the method 'evaluate'.
// That means 'y := f(x)' and 'y' is returned. It is only important that 'RHSFunction' knows the function
// space ('FuncSpace') that it is part from. (f \in FunctionSpace)

// eps = 0.001 => H Ref = 16
// eps = 0.002 => H Ref = 14
// eps = 0.004 => H Ref = 12
// eps = 0.008 => H Ref = 10
// eps = 0.016 => H Ref = 8
// eps = 0.032 => H Ref = 6
// eps = 0.064 => H Ref = 4
// eps = 0.128 => H Ref = 2
// eps = 0.256 => H Ref = 0

// NOTE that (delta/epsilon_est) needs to be a positive integer!

// is an exact solution available?
// #define EXACTSOLUTION_AVAILABLE

// Note that in the following, 'Imp' abbreviates 'Implementation'
namespace Problem {
// description see below 0.05
static const double EPSILON = 0.05;
static const double EPSILON_EST = 0.1;
static const double DELTA = 0.1;

// model problem information
struct ModelProblemData
  : public IModelProblemData
{
  ModelProblemData(const std::string filename = "no_name")
    : IModelProblemData(Constants(0.05, 0.0, 0.0), filename) {
  }

  inline int get_Number_of_Model_Problem() const {
    return 9;
  }

  inline void getMacroGridFile(std::string& macroGridName) const {
    // name and location of the grid file that describes the macro-grid:
    macroGridName = ("../dune/multiscale/grids/macro_grids/elliptic/cube_three.dgf");      // _strange_grid
  }

  // get the (starting) grid refinement level for solving the reference problem
  // in genereal, this is the smallest integer (level), so that solving the reference problem on this level,
  // yields a higly accurate approximation of the exact solution
  // ( here we have heterogenious reference problem, therefore we need a high refinement level )
  inline int getRefinementLevelReferenceProblem() const {
    // required refinement level for a fine scale reference solution
    // (a saved/precomputed solution is either already available for this level or it must be computed with the
    // following refinement level)
    return 10;
  }
};

// !FirstSource defines the right hand side (RHS) of the governing problem (i.e. it defines 'f').
// The value of the right hand side (i.e. the value of 'f') at 'x' is accessed by the method 'evaluate'.
// That means 'y := f(x)' and 'y' is returned. It is only important that 'RHSFunction' knows the function
// space ('FuncSpace') that it is part from. (f \in FunctionSpace)

template< class FunctionSpaceImp >
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
  inline void evaluate(const DomainType& x,
                       RangeType& y) const {
    #if 1
    double coefficient_0 = 2.0 * ( 1.0 / (8.0 * M_PI * M_PI) ) * ( 1.0 / ( 2.0 + cos( 2.0 * M_PI * (x[0] / EPSILON) ) ) );
    double coefficient_1 = ( 1.0 / (8.0 * M_PI * M_PI) ) * ( 1.0 + ( 0.5 * cos( 2.0 * M_PI * (x[0] / EPSILON) ) ) );

    double d_x0_coefficient_0
      = pow(2.0 + cos( 2.0 * M_PI * (x[0] / EPSILON) ), -2.0) * ( 1.0 / (2.0 * M_PI) ) * (1.0 / EPSILON) * sin(
      2.0 * M_PI * (x[0] / EPSILON) );

    JacobianRangeType grad_u;
    grad_u[0][0] = 2.0* M_PI* cos(2.0 * M_PI * x[0]) * sin(2.0 * M_PI * x[1]);
    grad_u[0][1] = 2.0* M_PI* sin(2.0 * M_PI * x[0]) * cos(2.0 * M_PI * x[1]);

    grad_u[0][0] += (-1.0) * EPSILON * M_PI
                    * ( sin(2.0 * M_PI * x[0]) * sin(2.0 * M_PI * x[1]) * sin( 2.0 * M_PI * (x[0] / EPSILON) ) );
    grad_u[0][0] += M_PI * ( cos(2.0 * M_PI * x[0]) * sin(2.0 * M_PI * x[1]) * cos( 2.0 * M_PI * (x[0] / EPSILON) ) );

    grad_u[0][1] += EPSILON * M_PI
                    * ( cos(2.0 * M_PI * x[0]) * cos(2.0 * M_PI * x[1]) * sin( 2.0 * M_PI * (x[0] / EPSILON) ) );

    RangeType d_x0_x0_u(0.0);
    d_x0_x0_u -= 4.0 * pow(M_PI, 2.0) * sin(2.0 * M_PI * x[0]) * sin(2.0 * M_PI * x[1]);
    d_x0_x0_u -= 2.0
                 * pow(M_PI,
                       2.0) * ( EPSILON + (1.0 / EPSILON) ) * cos(2.0 * M_PI * x[0]) * sin(2.0 * M_PI * x[1]) * sin(
      2.0 * M_PI * (x[0] / EPSILON) );
    d_x0_x0_u -= 4.0
                 * pow(M_PI, 2.0) * sin(2.0 * M_PI * x[0]) * sin(2.0 * M_PI * x[1]) * cos( 2.0 * M_PI * (x[0] / EPSILON) );

    RangeType d_x1_x1_u(0.0);
    d_x1_x1_u -= 4.0 * pow(M_PI, 2.0) * sin(2.0 * M_PI * x[0]) * sin(2.0 * M_PI * x[1]);
    d_x1_x1_u -= 2.0
                 * pow(M_PI,
                       2.0) * EPSILON
                 * cos(2.0 * M_PI * x[0]) * sin(2.0 * M_PI * x[1]) * sin( 2.0 * M_PI * (x[0] / EPSILON) );

    y = 0.0;
    y -= d_x0_coefficient_0 * grad_u[0][0];
    y -= coefficient_0 * d_x0_x0_u;
    y -= coefficient_1 * d_x1_x1_u;
    #endif // if 1

    // y = sin( 2.0 * M_PI * x[0] ) * sin( 2.0 * M_PI * x[1] );
  } // evaluate

  inline void evaluate(const DomainType& x,
                       const TimeType& time,
                       RangeType& y) const {
    evaluate(x, y);
  }
};

// ! default class for the second source term G.

// Realization: set G(x) = 0:
template< class FunctionSpaceImp >
class SecondSource
  : public Dune::Fem::Function< FunctionSpaceImp, SecondSource< FunctionSpaceImp > >
{
public:
  typedef FunctionSpaceImp FunctionSpaceType;

private:
  typedef SecondSource< FunctionSpaceType >                  ThisType;
  typedef Dune::Fem::Function< FunctionSpaceType, ThisType > BaseType;

public:
  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType  RangeType;

  static const int dimDomain = DomainType::dimension;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType  RangeFieldType;

public:
  inline void evaluate(const DomainType& x,
                       RangeType& y) const {
    y[0] = 0;
  }

  inline void evaluate(const int i, const DomainType& x,
                       RangeType& y) const {
    y[0] = 0;
  }
};

// the (non-linear) diffusion operator A^{\epsilon}(x,\xi)
// A^{\epsilon} : R^d -> R^d

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
  // in the linear setting, use the structure
  // A^{\epsilon}_i(x,\xi) = A^{\epsilon}_{i1}(x) \xi_1 + A^{\epsilon}_{i2}(x) \xi_2
  // the usage of an evaluate method with "evaluate ( i, j, x, y, z)" should be avoided
  // use "evaluate ( i, x, y, z)" instead and return RangeType-vector.

  // instantiate all possible cases of the evaluate-method:

  // (diffusive) flux = A^{\epsilon}( x , gradient_of_a_function )
  void diffusiveFlux(const DomainType& x,
                     const JacobianRangeType& gradient,
                     JacobianRangeType& flux) const {
    double coefficient_0 = 2.0 * ( 1.0 / (8.0 * M_PI * M_PI) ) * ( 1.0 / ( 2.0 + cos( 2.0 * M_PI * (x[0] / EPSILON) ) ) );
    double coefficient_1 = ( 1.0 / (8.0 * M_PI * M_PI) ) * ( 1.0 + ( 0.5 * cos( 2.0 * M_PI * (x[0] / EPSILON) ) ) );

    double stab = 0.0;

    flux[0][0] = coefficient_0 * gradient[0][0] + stab * gradient[0][1];
    flux[0][1] = coefficient_1 * gradient[0][1] + stab * gradient[0][0];
  } // diffusiveFlux

  // the jacobian matrix (JA^{\epsilon}) of the diffusion operator A^{\epsilon} at the position "\nabla v" in direction
  // "nabla w", i.e.
  // jacobian diffusiv flux = JA^{\epsilon}(\nabla v) nabla w:

  // jacobianDiffusiveFlux = A^{\epsilon}( x , position_gradient ) direction_gradient
  void jacobianDiffusiveFlux(const DomainType& x,
                             const JacobianRangeType& position_gradient,
                             const JacobianRangeType& direction_gradient,
                             JacobianRangeType& flux) const {
    double coefficient_0 = 2.0 * ( 1.0 / (8.0 * M_PI * M_PI) ) * ( 1.0 / ( 2.0 + cos( 2.0 * M_PI * (x[0] / EPSILON) ) ) );
    double coefficient_1 = ( 1.0 / (8.0 * M_PI * M_PI) ) * ( 1.0 + ( 0.5 * cos( 2.0 * M_PI * (x[0] / EPSILON) ) ) );

    flux[0][0] = coefficient_0 * direction_gradient[0][0];
    flux[0][1] = coefficient_1 * direction_gradient[0][1];

    // std :: cout << "Do not use this evaluate method." << std :: endl;
    // abort();
  } // jacobianDiffusiveFlux

// deprecated
// ! defaults (not to be used):
// z_i = A^{\epsilon}_i(x,vec)
// instantiate all possible cases of the evaluate-method:
  inline void evaluate(const int i,
                       const int j,
                       const DomainType& x,
                       RangeType& z) const {
    if (i == j)
    {
      if (i == 0)
      { z = 2.0 * ( 1.0 / (8.0 * M_PI * M_PI) ) * ( 1.0 / ( 2.0 + cos( 2.0 * M_PI * (x[0] / EPSILON) ) ) ); } else
      { z = ( 1.0 / (8.0 * M_PI * M_PI) ) * ( 1.0 + ( 0.5 * cos( 2.0 * M_PI * (x[0] / EPSILON) ) ) ); }
    } else
    { z = 0.0; }

    std::cout << "Do not use this evaluate method." << std::endl;
    abort();
  } // evaluate

  inline void evaluate(const int i,
                       const int j,
                       const DomainType& x,
                       const DomainType& y,
                       RangeType& z) const {
    std::cout
    <<
    "WARNING! Inadmissible call for 'evaluate' method of the Diffusion class! See 'problem_specification.hh' for details."
    << std::endl;

    std::abort();

    z = 0.0;
  } // evaluate

  inline void evaluate(const int i,
                       const int j,
                       const DomainType& x,
                       const TimeType& time,
                       RangeType& z) const {
    std::cout
    << "WARNING! Call for 'evaluate' method of the Diffusion class with time variable! Skip to standard evaluation."
    << std::endl;

    std::abort();

    return evaluate(i, j, x, z);
  } // evaluate

  // dummy implementation
  inline void evaluate(const DomainType& x,
                       RangeType& y) const {
    std::cout
    <<
    "WARNING! Wrong call for 'evaluate' method of the Diffusion class (evaluate(x,y)). This is just a dummy method. Use 'diffusiveFlux(...)' instead."
    << std::endl;
    std::abort();
  } // evaluate
};

template< class FunctionSpaceImp, class FieldMatrixImp >
class HomDiffusion
  : public Dune::Fem::Function< FunctionSpaceImp, HomDiffusion< FunctionSpaceImp, FieldMatrixImp > >
{
public:
  typedef FunctionSpaceImp FunctionSpaceType;
  typedef FieldMatrixImp   FieldMatrixType;

private:
  typedef HomDiffusion< FunctionSpaceType, FieldMatrixType > ThisType;
  typedef Dune::Fem::Function< FunctionSpaceType, ThisType > BaseType;

public:
  typedef typename FunctionSpaceType::DomainType        DomainType;
  typedef typename FunctionSpaceType::RangeType         RangeType;
  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType  RangeFieldType;

  typedef DomainFieldType TimeType;

public:
  FieldMatrixType* A_hom_;

  #if 1

public:
  inline explicit HomDiffusion(FieldMatrixType& A_hom)
    : A_hom_(&A_hom)
  {}

  #endif // if 1

  // in the linear setting, use the structure
  // A^{\epsilon}_i(x,\xi) = A^{\epsilon}_{i1}(x) \xi_1 + A^{\epsilon}_{i2}(x) \xi_2
  // the usage of an evaluate method with "evaluate ( i, j, x, y, z)" should be avoided
  // use "evaluate ( i, x, y, z)" instead and return RangeType-vector.

  // instantiate all possible cases of the evaluate-method:

  // (diffusive) flux = A^{\epsilon}( x , gradient_of_a_function )
  void diffusiveFlux(const DomainType& x,
                     const JacobianRangeType& gradient,
                     JacobianRangeType& flux) const {
    flux[0][0] = (*A_hom_)[0][0] * gradient[0][0] + (*A_hom_)[0][1] * gradient[0][1];
    flux[0][1] = (*A_hom_)[1][0] * gradient[0][0] + (*A_hom_)[1][1] * gradient[0][1];
  }

  // the jacobian matrix (JA^{\epsilon}) of the diffusion operator A^{\epsilon} at the position "\nabla v" in direction
  // "nabla w", i.e.
  // jacobian diffusiv flux = JA^{\epsilon}(\nabla v) nabla w:

  // jacobianDiffusiveFlux = A^{\epsilon}( x , position_gradient ) direction_gradient
  void jacobianDiffusiveFlux(const DomainType& x,
                             const JacobianRangeType& position_gradient,
                             const JacobianRangeType& direction_gradient,
                             JacobianRangeType& flux) const {
    #ifdef LINEAR_PROBLEM
    std::cout << "Not yet implemented." << std::endl;
    std::abort();
    #else // ifdef LINEAR_PROBLEM
    std::cout << "Nonlinear example not yet implemented." << std::endl;
    std::abort();
    #endif // ifdef LINEAR_PROBLEM
  } // jacobianDiffusiveFlux

  inline void evaluate(const int i,
                       const int j,
                       const DomainType& x,
                       const DomainType& y,
                       RangeType& z) const {
    std::cout
    <<
    "WARNING! Inadmissible call for 'evaluate' method of the Diffusion class! See 'problem_specification.hh' for details."
    << std::endl;

    std::abort();

    z = 0.0;
  } // evaluate

// deprecated
// ! defaults (not to be used):
// z_i = A^{\epsilon}_i(x,vec)
// instantiate all possible cases of the evaluate-method:
  inline void evaluate(const int i,
                       const int j,
                       const DomainType& x,
                       RangeType& z) const {
    std::cout
    <<
    "WARNING! Inadmissible call for 'evaluate' method of the Diffusion class! See 'problem_specification.hh' for details."
    << std::endl;

    std::abort();
  } // evaluate

  inline void evaluate(const int i,
                       const int j,
                       const DomainType& x,
                       const TimeType& time,
                       RangeType& z) const {
    std::cout
    << "WARNING! Call for 'evaluate' method of the Diffusion class with time variable! Skip to standard evaluation."
    << std::endl;

    std::abort();

    return evaluate(i, j, x, z);
  } // evaluate

  // dummy implementation
  inline void evaluate(const DomainType& x,
                       RangeType& y) const {
    std::cout
    <<
    "WARNING! Wrong call for 'evaluate' method of the Diffusion class (evaluate(x,y)). This is just a dummy method. Use 'diffusiveFlux(...)' instead."
    << std::endl;
    std::abort();
  } // evaluate
};

// define the mass term:
template< class FunctionSpaceImp >
class MassTerm
  : public Dune::Fem::Function< FunctionSpaceImp, MassTerm< FunctionSpaceImp > >
{
public:
  typedef FunctionSpaceImp FunctionSpaceType;

private:
  typedef MassTerm< FunctionSpaceType >                      ThisType;
  typedef Dune::Fem::Function< FunctionSpaceType, ThisType > BaseType;

public:
  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType  RangeType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType  RangeFieldType;

  typedef DomainFieldType TimeType;

public:
  inline void evaluate(const DomainType& x,
                       RangeType& y) const {
    y[0] = 0.00001;
  }

  // dummy implementation
  inline void evaluate(const DomainType& x,
                       const TimeType time,
                       RangeType& y) const {
    std::cout << "WARNING! Wrong call for 'evaluate' method of the MassTerm class (evaluate(x,t,y)). Return 0.0."
              << std::endl;
    return evaluate(x, y);
  }
};

// a dummy function class for functions, vectors and matrices
template< class FunctionSpaceImp >
class DefaultDummyFunction
  : public Dune::Fem::Function< FunctionSpaceImp, DefaultDummyFunction< FunctionSpaceImp > >
{
public:
  typedef FunctionSpaceImp FunctionSpaceType;

private:
  typedef DefaultDummyFunction< FunctionSpaceType >          ThisType;
  typedef Dune::Fem::Function< FunctionSpaceType, ThisType > BaseType;

public:
  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType  RangeType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType  RangeFieldType;

  typedef DomainFieldType TimeType;

public:
  inline void evaluate(const int i,
                       const int j,
                       const DomainType& x,
                       const DomainType& y,
                       RangeType& z) const {
    z = 0;
  }

  inline void evaluate(const int i,
                       const DomainType& x,
                       const DomainType& y,
                       RangeType& z) const {
    z = 0;
  }

  inline void evaluate(const int i,
                       const int j,
                       const DomainType& x,
                       RangeType& y) const {
    y = 0;
  }

  inline void evaluate(const int i,
                       const DomainType& x,
                       RangeType& y) const {
    y = 0;
  }

  inline void evaluate(const int i,
                       const DomainType& x,
                       const TimeType& t,
                       RangeType& y) const {
    y = 0;
  }

  // dummy implementation
  inline void evaluate(const DomainType& x,
                       RangeType& y) const {
    y = 0;
  }

  // dummy implementation
  inline void evaluate(const DomainType& x,
                       const TimeType time,
                       RangeType& y) const {
    y = 0;
  }

  inline void getTimeStepSize(double& time_step_size) const {
    time_step_size = 0;
  }
};

// ! Exact solution (typically it is unknown)
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
  // in case 'u' has NO time-dependency use the following method:
  inline void evaluate(const DomainType& x,
                       RangeType& y) const {
    // approximation obtained by homogenized solution + first corrector
    y = sin(2.0 * M_PI * x[0]) * sin(2.0 * M_PI * x[1]);         // coarse part
    y += 0.5 * EPSILON * ( cos(2.0 * M_PI * x[0]) * sin(2.0 * M_PI * x[1]) * sin( 2.0 * M_PI * (x[0] / EPSILON) ) );          //
                                                                                                                              //
                                                                                                                              //
                                                                                                                              //
                                                                                                                              //
                                                                                                                              //fine
                                                                                                                              //
                                                                                                                              //
                                                                                                                              //
                                                                                                                              //
                                                                                                                              //part
    // || u_fine_part ||_L2 = 0.00883883 (fuer eps = 0.05 )
  } // evaluate

  inline void evaluateJacobian(const DomainType& x, JacobianRangeType& grad_u) const {
    grad_u[0][0] = 2.0* M_PI* cos(2.0 * M_PI * x[0]) * sin(2.0 * M_PI * x[1]);
    grad_u[0][1] = 2.0* M_PI* sin(2.0 * M_PI * x[0]) * cos(2.0 * M_PI * x[1]);

    grad_u[0][0] += (-1.0) * EPSILON * M_PI
                    * ( sin(2.0 * M_PI * x[0]) * sin(2.0 * M_PI * x[1]) * sin( 2.0 * M_PI * (x[0] / EPSILON) ) );
    grad_u[0][0] += M_PI * ( cos(2.0 * M_PI * x[0]) * sin(2.0 * M_PI * x[1]) * cos( 2.0 * M_PI * (x[0] / EPSILON) ) );

    grad_u[0][1] += EPSILON * M_PI
                    * ( cos(2.0 * M_PI * x[0]) * cos(2.0 * M_PI * x[1]) * sin( 2.0 * M_PI * (x[0] / EPSILON) ) );
  } // evaluateJacobian

  // in case 'u' HAS a time-dependency use the following method:
  // unfortunately GRAPE requires both cases of the method 'evaluate' to be
  // instantiated
  inline void evaluate(const DomainType& x,
                       const TimeType& timedummy,
                       RangeType& y) const {
    evaluate(x, y);
  }
};
}

// we need to know the term 'abstract class'.

// In short: An abstract class is only created to be a 'base class' for a set of other classes (the so called 'derived
// classes').
// These derived classes typically share a number of methods that can be subsumed by the abstract class.
// Abstract classes contain at least one virtual method, that means a method of the kind: 'virtual "datatype"
// methodname() = 0;'.
// Examples for virtual methods are: 'virtual void evaluate() = 0;' or 'virtual integer sum() = 0;'.
// Virtual methods can not be used, they will lead to error prompts! Therefor it is impossible to create objects of an
// abstract class.
// To use such a method nevertheless, the virtual method must be inherited and overwhrighten by an equally named method
// of a derived class.

#endif // ifndef DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH
