#ifndef DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_NINE
#define DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_NINE

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
namespace Nine {
// description see below 0.05
CONSTANTSFUNCTION(0.05, 0.1, 0.1)
// model problem information
struct ModelProblemData
  : public IModelProblemData
{
  ModelProblemData(const std::string filename = "no_name")
    : IModelProblemData(constants(), filename) {
  }

  inline int get_Number_of_Model_Problem() const {
    return 9;
  }

  //! \copydoc IModelProblemData::getMacroGridFile()
  inline std::string getMacroGridFile() const {
    return("../dune/multiscale/grids/macro_grids/elliptic/cube_three.dgf");      // _strange_grid
  }

  //! \copydoc IModelProblemData::getRefinementLevelReferenceProblem()
  inline int getRefinementLevelReferenceProblem() const {
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
  FirstSource(){}

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

  /** \brief default class for the second source term G.
   * Realization: set G(x) = 0: **/
  NULLFUNCTION(SecondSource)

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
    Diffusion(){}

  // in the linear setting, use the structure
  // A^{\epsilon}_i(x,\xi) = A^{\epsilon}_{i1}(x) \xi_1 + A^{\epsilon}_{i2}(x) \xi_2
  // the usage of an evaluate method with "evaluate ( i, j, x, y, z)" should be avoided
  // use "evaluate ( i, x, y, z)" instead and return RangeType-vector.

  // instantiate all possible cases of the evaluate-method:

  // (diffusive) flux = A^{\epsilon}( x , gradient_of_a_function )
  void diffusiveFlux(const DomainType& x,
                     const JacobianRangeType& gradient,
                     JacobianRangeType& flux) const {
    double coefficient_0 = 2.0 * ( 1.0 / (8.0 * M_PI * M_PI) ) * ( 1.0 / ( 2.0 + cos( 2.0 * M_PI * (x[0] / constants().epsilon) ) ) );
    double coefficient_1 = ( 1.0 / (8.0 * M_PI * M_PI) ) * ( 1.0 + ( 0.5 * cos( 2.0 * M_PI * (x[0] / constants().epsilon) ) ) );

    double stab = 0.0;

    flux[0][0] = coefficient_0 * gradient[0][0] + stab * gradient[0][1];
    flux[0][1] = coefficient_1 * gradient[0][1] + stab * gradient[0][0];
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

public:
  inline explicit HomDiffusion(FieldMatrixType& A_hom)
    : A_hom_(&A_hom)
  {}

  // in the linear setting, use the structure
  // A^{\epsilon}_i(x,\xi) = A^{\epsilon}_{i1}(x) \xi_1 + A^{\epsilon}_{i2}(x) \xi_2
  // the usage of an evaluate method with "evaluate ( i, j, x, y, z)" should be avoided
  // use "evaluate ( i, x, y, z)" instead and return RangeType-vector.

  // instantiate all possible cases of the evaluate-method:

  // (diffusive) flux = A^{\epsilon}( x , gradient_of_a_function )
  void diffusiveFlux(const DomainType& /*x*/,
                     const JacobianRangeType& gradient,
                     JacobianRangeType& flux) const {
    flux[0][0] = (*A_hom_)[0][0] * gradient[0][0] + (*A_hom_)[0][1] * gradient[0][1];
    flux[0][1] = (*A_hom_)[1][0] * gradient[0][0] + (*A_hom_)[1][1] * gradient[0][1];
  }

  /** the jacobian matrix (JA^{\epsilon}) of the diffusion operator A^{\epsilon} at the position "\nabla v" in direction
   * "nabla w", i.e.
   * jacobian diffusiv flux = JA^{\epsilon}(\nabla v) nabla w:
   * jacobianDiffusiveFlux = A^{\epsilon}( x , position_gradient ) direction_gradient 
  **/
  void jacobianDiffusiveFlux(const DomainType&,
                             const JacobianRangeType&,
                             const JacobianRangeType&,
                             JacobianRangeType&) const {
    DUNE_THROW(Dune::NotImplemented,"");
  } // jacobianDiffusiveFlux

  template < class... Args >
  void evaluate( Args... ) const
  {
    DUNE_THROW(Dune::NotImplemented, "Inadmissible call for 'evaluate'");
  }
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
  MassTerm(){}

  inline void evaluate(const DomainType& /*x*/,
                       RangeType& y) const {
    y[0] = 0.00001;
  }

  // dummy implementation
  inline void evaluate(const DomainType& x,
                       const TimeType /*time*/,
                       RangeType& y) const {
    DSC_LOG_ERROR << "WARNING! Wrong call for 'evaluate' method of the MassTerm class (evaluate(x,t,y)). Return 0.0."
              << std::endl;
    return evaluate(x, y);
  }
};

//! a dummy function class for functions, vectors and matrices
NULLFUNCTION(DefaultDummyFunction)
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
  ExactSolution(){}

  // in case 'u' has NO time-dependency use the following method:
  inline void evaluate(const DomainType& x,
                       RangeType& y) const {
    // approximation obtained by homogenized solution + first corrector
    y = sin(2.0 * M_PI * x[0]) * sin(2.0 * M_PI * x[1]);         // coarse part
    y += 0.5 * constants().epsilon * ( cos(2.0 * M_PI * x[0]) * sin(2.0 * M_PI * x[1]) * sin( 2.0 * M_PI * (x[0] / constants().epsilon) ) );          //
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
} // namespace Nine {
}



#endif // ifndef DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_NINE
