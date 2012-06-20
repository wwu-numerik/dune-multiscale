#ifndef DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH
#define DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH

#include <dune/fem/function/common/function.hh>

// ! USE THIS ONE FOR MSFEM TESTS! PURELY LINEAR ELLIPTIC!

// ! NOTE: MODEL PROBLEM 0!!!!!

// ! For more further details about the implementation of the following classes, see the end of the file

// in general we regard problems of the following type:

// - div ( A^{\epsilon} (x,\nabla u^{\epsilon}) ) + m^{\epsilon} u^{\epsilon} = f - div G

// Here we have:
// u^{\epsilon} = exact solution of the problem
// A^{\epsilon} = diffusion matrix (or tensor) with the structure A^{\epsilon}(x) = A(x,\frac{x}{\epsilon})
// m^{\epsilon} = a mass term (or reaction term) with the structure m^{\epsilon}(x) = m(x,\frac{x}{\epsilon})
// f = first source term with the structure f = f(x) (=> no micro-scale dependency)
// G = second source term with the structure G = G(x) (=> no micro-scale dependency). Note that 'G' is directly
// implemented! ( we do not implement '- div G'!)

// Note, that A^{\epsilon} is a monotone operator

// !############################## Elliptic Problem -1 ###################################

// ! we define:

// The entries of the operator A^{\epsilon} by
// ! a^{\epsilon}_{1}(x_1,x_2) := (**y_1**,**y_2**)
// ! a^{\epsilon}_{2}(x_1,x_2) := (**y_1**,**y_2**)

// The mass (or reaction) term m^{\epsilon} is given by:
// ! m^{\epsilon} := \epsilon
// Since \epsilon tends to zero, we may say that we do not have a real mass term for our problem. It is a simple
// condition to fix the solution which is only unique up to a constant. In fact we still approximate the solution of the
// problem without mass.

// The first source term f is given by:
// ! f(x) := ****
// since the second source is zero, f will form the right hand side (RHS) of our discrete problem

// The second source term G is constantly zero:
// ! G(x) := 0

// !FirstSource defines the right hand side (RHS) of the governing problem (i.e. it defines 'f').
// The value of the right hand side (i.e. the value of 'f') at 'x' is accessed by the method 'evaluate'. That means 'y
// := f(x)' and 'y' is returned. It is only important that 'RHSFunction' knows the function space ('FuncSpace') that it
// is part from. (f \in FunctionSpace)

// NOTE that (delta/epsilon_est) needs to be a positive integer!

// is an exact solution available?
// #define EXACTSOLUTION_AVAILABLE

// Note that in the following, 'Imp' abbreviates 'Implementation'
namespace Problem {
// description see below 0.05
static const double EPSILON = 1.0;
static const double EPSILON_EST = 1.0;
static const double DELTA = 1.0;

// model problem information
class ModelProblemData
{
protected:
  // name of the file where data is saved
  const std::string file_name_;

  int current_number_of_cell_problem_;

public:
  // Constructor for ModelProblemData
  inline explicit ModelProblemData(const std::string& file_name)
    : file_name_(file_name)
      , current_number_of_cell_problem_(-1)
  {}

  inline explicit ModelProblemData()
    : file_name_("no_name")
      , current_number_of_cell_problem_(-1)
  {}

public:
  inline int get_Number_of_Model_Problem() const {
    return -1;
  }

  // epsilon (the smaller epsilon, the finer the micro-structure)
  // in the periodic setting, epsilon denotes the periode of the fine-scale oscillations
  // in the non-periodic setting, can be seen as a representative size for the fine-scale behaviour
  inline double getEpsilon() const {
    const double epsilon = EPSILON;

    return epsilon;
  }

  // epsilon (the smaller epsilon, the finer the micro-structure)
  inline double getEpsilonEstimated() const {
    const double epsilon_est = EPSILON_EST;

    return epsilon_est;
  }

  // edge length of a cell (where we solve the cell problems)
  // we need delta >= epsilon
  inline double getDelta() const {
    const double delta = DELTA;

    return delta;
    // NOTE that (delta/epsilon_est) needs to be a positive integer!
  }

  inline void getMacroGridFile(std::string& macroGridName) const {
    // name and location of the grid file that describes the macro-grid:
    std::string macro_grid_location("../dune/multiscale/grids/macro_grids/elliptic/cube_three.dgf");

    macroGridName = macro_grid_location;
  }

  // get an information on whether we use the solutions of cell problems that are already computed and saved in a file
  // with the name 'name_'
  inline std::string getName_and_getBool(bool& use_saved) const {
    if (file_name_ == "no_name")
    { use_saved = false; } else
    { use_saved = true; }

    return file_name_;
  } // getName_and_getBool

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

  inline void set_current_number_of_cell_problem(int number) {
    current_number_of_cell_problem_ = number;
  }

  inline int get_current_number_of_cell_problem() {
    return current_number_of_cell_problem_;
  }
};

// !FirstSource defines the right hand side (RHS) of the governing problem (i.e. it defines 'f').
// The value of the right hand side (i.e. the value of 'f') at 'x' is accessed by the method 'evaluate'. That means 'y
// := f(x)' and 'y' is returned. It is only important that 'RHSFunction' knows the function space ('FuncSpace') that it
// is part from. (f \in FunctionSpace)

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
  typedef typename FunctionSpaceType::DomainType        DomainType;
  typedef typename FunctionSpaceType::RangeType         RangeType;
  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

  static const int dimDomain = DomainType::dimension;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType  RangeFieldType;

  typedef DomainFieldType TimeType;

public:
  inline void evaluate(const DomainType& x,
                       RangeType& y) const {
    double a_0_x_0 = 1.0 + pow(x[0], 2.0);
    double a_1_x_1 = 1.0 + pow(x[0], 2.0);

    double grad_a_0_x_0 = 2.0 * x[0];
    double grad_a_1_x_1 = 0.0;

    JacobianRangeType grad_u(0.0);

    grad_u[0][0] = (1.0 - x[0]) * (1.0 - x[1]) * x[1] - x[0] * (1.0 - x[1]) * x[1];
    grad_u[0][1] = x[0] * (1.0 - x[0]) * (1.0 - x[1]) - x[0] * (1.0 - x[0]) * x[1];

    RangeType u = x[0] * (1.0 - x[0]) * (1.0 - x[1]) * x[1];
    RangeType d_xx_u = (-2.0) * (1.0 - x[1]) * x[1];
    RangeType d_yy_u = (-2.0) * (1.0 - x[0]) * x[0];

    y = 0.0;
    y -= grad_a_0_x_0 * grad_u[0][0];
    y -= a_0_x_0 * d_xx_u;

    y -= grad_a_1_x_1 * grad_u[0][1];
    y -= a_1_x_1 * d_yy_u;
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
    double a_0 = 1.0 + pow(x[0], 2.0);

    flux[0][0] = a_0 * gradient[0][0];
    flux[0][1] = a_0 * gradient[0][1];
  }

// deprecated
// ! defaults (not to be used):
// z_i = A^{\epsilon}_i(x,vec)
// instantiate all possible cases of the evaluate-method:
  inline void evaluate(const int i,
                       const int j,
                       const DomainType& x,
                       RangeType& z) const {
    if (i == j)
    { z = 0.0; } else
    { z = 0.0; }
    abort();
  }

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
  // in case 'u' has NO time-dependency use the following method:
  inline void evaluate(const DomainType& x,
                       RangeType& y) const {
    y = x[0] * (1.0 - x[0]) * (1.0 - x[1]) * x[1];
  }

  // in case 'u' HAS a time-dependency use the following method:
  // unfortunately GRAPE requires both cases of the method 'evaluate' to be
  // instantiated
  inline void evaluate(const DomainType& x,
                       const TimeType& timedummy,
                       RangeType& y) const {
    evaluate(x, y);
  }

  inline void evaluateJacobian(const DomainType& x, JacobianRangeType& grad_u) const {
    grad_u[0][0] = (1.0 - x[0]) * (1.0 - x[1]) * x[1] - x[0] * (1.0 - x[1]) * x[1];
    grad_u[0][1] = x[0] * (1.0 - x[0]) * (1.0 - x[1]) - x[0] * (1.0 - x[0]) * x[1];
  }
};
}

// we need to know the term 'abstract class'.

// In short: An abstract class is only created to be a 'base class' for a set of other classes (the so called 'derived
// classes'). These derived classes typically share a number of methods that can be subsumed by the abstract class.
// Abstract classes contain at least one virtual method, that means a method of the kind: 'virtual "datatype"
// methodname() = 0;'.
// Examples for virtual methods are: 'virtual void evaluate() = 0;' or 'virtual integer sum() = 0;'.
// Virtual methods can not be used, they will lead to error prompts! Therefor it is impossible to create objects of an
// abstract class. To use such a method nevertheless, the virtual method must be inherited and overwhrighten by an
// equally named method of a derived class.

// Example: We construct a class 'FirstSource' that inherits the features/methods of the abstract class 'Function'. In
// order to fill for instance the virtual method 'evaluate' (of //the abstract class 'Function') with life, we need to
// define an equally named method 'evaluate' in 'FirstSource'.
// Unfortunately, if used very often, virtual methods cause long running times. Therefore, to avoid virtual methods for
// the sake of efficiency, Dune-programmers often make use of the so called Barton-Nackman trick: the derived class
// 'FirstSource' inherits the (virtual) methods of the class 'Function< FuncSpace, FirstSource >'. But the class
// template
// 'Function' recieves the template parameter 'FirstSource', which itself is derived of the class 'Function< FuncSpace,
// FirstSource >'. The snake seems to bite its tail! But it won't cause any errors. Simplified, the following happens:
// if
// the program wants to execute the method 'evaluate(...)' it is normaly searching in 'Function', there it finds out
// that
// the method is virtual and it has to search in the derived class 'FirstSource'. Such a process is not very efficient.
// By using the Barton-Nackman trick this process of 'searching' is reversed, the base class 'turns into' the derived
// class and the derived class 'turns into' the base class. Since the program always searches at first in the base class
// (now: 'FirstSource') it will immediately find a non-virtual method 'evaluate()' that can be executed.
// To realize the Barton-Nackman trick we essentially need the order 'asimp()'. This order must be somewhere in the
// orginal base class ('Function'), before the virtual methods! It is responsible for the reversion of base class and
// derived class.

#endif // ifndef DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH
