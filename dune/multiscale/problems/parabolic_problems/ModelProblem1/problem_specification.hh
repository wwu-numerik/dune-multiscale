#ifndef DUNE_PARABOLIC_MODEL_PROBLEM_SPECIFICATION_HH
#define DUNE_PARABOLIC_MODEL_PROBLEM_SPECIFICATION_HH

#include <dune/fem/function/common/function.hh>
#include <vector>

// ! This is a heterogeneous test problem. (heterogeneous diffusion matrix)
// No heterogenity in the advective part so that we can use that b is divergence-free and with zero meanvalue. This
// means that we can make a very fine scale
// calculation to see what happens. The heterogenity will be in the diffusive part.

// !#####################################################################################################
// !######################### HM-FEM Parabolic Problem 6 - DIM = 2 ######################################
// !#####################################################################################################

// ! For more further details about the implementation of the following classes, see the end of the file

// in general we regard problems of the following type:

// \dt u^{\eps} - div ( A^{\epsilon} \nabla u^{\epsilon} ) + \frac{1}{\epsilon} b^{\epsilon} \nabla u^{\epsilon}
// + \frac{1}{\epsilon^2} c^{\epsilon} u^{\epsilon} = f - div G

// Here we have:
// u^{\epsilon}(t,x) = exact solution of the problem
// A^{\epsilon}(t,x) = diffusion matrix (or tensor) with the structure A^{\epsilon}(t,x) = A(t,\frac{x}{\epsilon})
// b^{\epsilon}(t,x) = advection term with the structure b^{\epsilon}(t,x) = b(t,\frac{x}{\epsilon})
// c^{\epsilon}(t,x) = a mass term (or reaction term) with the structure c^{\epsilon}(x) = c(x,\frac{x}{\epsilon})
// f(t,x) = first source term with the structure f = f(t,x) (=> no micro-scale dependency)
// G(t,x) = second source term with the structure G = G(t,x) (=> no micro-scale dependency). Note that 'G' is directly
// implemented! ( we do not implement '- div G'!)

// in all of the model problems, we ignore the reaction and all the sources, i.e.: c = f = G_i = 0

// ! we define:

// exact solution:
// ! Unknown! An approximation is determined by homogenizing the problem.

// The entries of the matrix A^{\epsilon}(x) are given by
// ! a_{1,1}(t,y_1,y_2) :=
// ! a_{2,2}(t,y_1,y_2) :=
// ! a_{1,2}(t,y_1,y_2) := a_{2,1}(y_1,y_2) = 0

// The mass term m^{\epsilon} ist given by:
// ! m^{\epsilon}(t,y_1,y_2) := 0
// Since \epsilon tends to zero, we may say that we do not have a real mass term for our problem. It is a simple
// condition to fix the solution which is only unique up to a constant. In fact we still approximate the solution of the
// problem without mass.

// The first source term f is given by:
// ! f(t,x) := 0
// since the second source is zero, f will form the right hand side (RHS) of our discrete problem

// The second source term G is constantly zero:
// ! G(t,x) := 0

// The advection term b ( where b^{\epsilon}(t,x)= \frac{1}{\epsilon} b(t,\frac{x}{\epsilon}) ) is given by:
// ! b_1(t,y_1,y_2) :=
// ! b_1(t,y_1,y_2) :=
// Note that b is divergence-free and that it has zero average. This implies that that there is no drift occuring.

// The initial value v_0 is given by:
// ! v_0(x_1,x_2) =

// Note that in the following, 'Imp' abbreviates 'Implementation'

namespace Dune {
// Diese Klasse kann man natürlich auch komplett über Makros definieren (ich persönlich finde es so übersichtlicher)
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
      current_number_of_cell_problem_(-1)
  {}

public:
  inline int get_Number_of_Model_Problem() const {
    return 7;
  }

  inline bool is_Model_Problem_homogenious() const {
    return false;
  }

  inline double getEpsilon() const {
    const double epsilon = 0.01;

    return epsilon;
  }

  inline double getEpsilonEstimated() const {
    const double epsilon_est = 0.01;

    return epsilon_est;
  }

  inline double getDelta() const {
    const double delta = 0.02;

    return delta;
  }

  inline double getMaxTime() const {
    const double max_T = 3.0;

    return max_T;
  }

  inline std::string getMacroGridFile() const {
    // name and location of the grid file that describes the macro-grid:
    macroGridName = ("grids/macro_grids/parabolic/huge_square2D.dgf");

    
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
  // ( here we have heterogenious reference problem, therefor we need a high uniform refinement level )
  inline int getRefinementLevelReferenceProblem() const {
    return 14;
  }

  inline void set_current_number_of_cell_problem(int number) {
    current_number_of_cell_problem_ = number;
  }

  inline int get_current_number_of_cell_problem() {
    return current_number_of_cell_problem_;
  }
};

// v_0 denotes the initial value.
template< class FunctionSpaceImp >
class InitialValue
  : public Function< FunctionSpaceImp, InitialValue< FunctionSpaceImp > >

    // we derive a class 'FirstSource' of the class 'Function< FunctionSpaceImp, FirstSource< FunctionSpaceImp > >'
    // Note: 'Function' is only a class-template, 'Function< FunctionSpaceImp, FirstSource< FunctionSpaceImp > >' is a
    // class!

    // template parameters of 'Function': - 'FunctionSpaceImp' the function space where the Function ('FirstSource')
    // shall be in, and - 'FirstSource': see Barton-Nackman trick.
{
public:
  typedef FunctionSpaceImp FunctionSpaceType;

private:
  typedef InitialValue< FunctionSpaceType >       ThisType;
  typedef Function< FunctionSpaceType, ThisType > BaseType;

public:
  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType  RangeType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType  RangeFieldType;

public:
  inline explicit InitialValue(const FunctionSpaceType& fSpace)
    : BaseType(fSpace)
  {}

  inline void evaluate(const DomainType& x, RangeType& y) const
  // x = (x_1,x_2), y = f(x_1,x_2)
  // Note: in C++, arrays are non-permissible return values. ('array evaluate(...)' causes error prompts)
  // Therefor we need to declare 'evaluate' like above to avoid such problems.
  {
    if ( (x[0] <= 0.5) && (x[0] >= -0.5) && (x[1] <= 0.5) && (x[1] >= -0.5) )
    {
      y = 2.0 * sin(2.0 * M_PI * x[0]) * cos(M_PI * x[1]);
    } else {
      y[0] = 0.0;
      y[1] = 0.0;
    }
  } // evaluate

  inline void jacobian(const DomainType& x, DomainType& y) const {
    if ( (x[0] <= 0.5) && (x[0] >= -0.5) && (x[1] <= 0.5) && (x[1] >= -0.5) )
    {
      y[0] = 4.0* M_PI* cos(2.0 * M_PI * x[0]) * cos(M_PI * x[1]);
      y[1] = (-2.0) * M_PI * sin(2.0 * M_PI * x[0]) * sin(M_PI * x[1]);
    } else {
      y[0] = 0.0;
      y[1] = 0.0;
    }
  } // jacobian
};

// !FirstSource defines the right hand side (RHS) of the governing problem (i.e. it defines 'f').
// The value of the right hand side (i.e. the value of 'f') at 'x' is accessed by the method 'evaluate'. That means 'y
// := f(x)' and 'y' is returned. It is only important that 'RHSFunction' knows the function space ('FuncSpace') that it
// is part from. (f \in FunctionSpace)

// ! default class for the second source term. Set f(x) = 0:
template< class FunctionSpaceImp >
class FirstSource
  : public Function< FunctionSpaceImp, FirstSource< FunctionSpaceImp > >

    // we derive a class 'FirstSource' of the class 'Function< FunctionSpaceImp, FirstSource< FunctionSpaceImp > >'
    // Note: 'Function' is only a class-template, 'Function< FunctionSpaceImp, FirstSource< FunctionSpaceImp > >' is a
    // class!

    // template parameters of 'Function': - 'FunctionSpaceImp' the function space where the Function ('FirstSource')
    // shall be in, and - 'FirstSource': see Barton-Nackman trick.
{
public:
  typedef FunctionSpaceImp FunctionSpaceType;

private:
  typedef FirstSource< FunctionSpaceType >        ThisType;
  typedef Function< FunctionSpaceType, ThisType > BaseType;

public:
  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType  RangeType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType  RangeFieldType;

  typedef DomainFieldType TimeType;

  // theta is a parameter that is only required for the cell problems. In our case it is completly ignored by defaulting
  // to 'theta = 1'
  const double theta_;
  const TimeType t_;

public:
  inline explicit FirstSource(const FunctionSpaceType& fSpace, const TimeType& t)
    : BaseType(fSpace)
      , theta_(1)
      , t_(t)
  {}

  inline explicit FirstSource(const FunctionSpaceType& fSpace)
    : BaseType(fSpace)
      , theta_(1)
      , t_(1)
  {}

  inline void evaluate(const DomainType& x, RangeType& y) const
  // x = (x_1,x_2), y = f(x_1,x_2)
  // Note: in C++, arrays are non-permissible return values. ('array evaluate(...)' causes error prompts)
  // Therefor we need to declare 'evaluate' like above to avoid such problems.
  {
    y[0] = 0;
  }

  inline void evaluate(const DomainType& x, const TimeType& t, RangeType& y) const
  // x = (x_1,x_2), y = f(x_1,x_2,t)
  {
    y[0] = 0;
  }
};

// ! default class for the second source term. Set G(x) = 0:
template< class FunctionSpaceImp >
class SecondSource
  : public Function< FunctionSpaceImp, SecondSource< FunctionSpaceImp > >
{
public:
  typedef FunctionSpaceImp FunctionSpaceType;

private:
  typedef SecondSource< FunctionSpaceType >       ThisType;
  typedef Function< FunctionSpaceType, ThisType > BaseType;

public:
  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType  RangeType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType  RangeFieldType;

public:
  inline explicit SecondSource(const FunctionSpaceType& fSpace)
    : BaseType(fSpace)
  {}

  inline void evaluate(const DomainType& x,
                       RangeType& y) const {
    y[0] = 0;
  }

  inline void evaluate(const int i, const DomainType& x,
                       RangeType& y) const {
    y[0] = 0;
  }
};

// ! The next class could be used to define the exact solution of the problem above. But since it very problematic to
// find this exact solution, the following problem is only a dummy. To compare (later on) the numerical solution with
// the
// exact solution in order to calculate the EOC (Estimated Order of Convergence), we will need to homogenize the problem
// to get a numerical approximation of the exact solution.

// ! Exact solution unknown
template< class FunctionSpaceImp >
class ExactSolution
  : public Function< FunctionSpaceImp, ExactSolution< FunctionSpaceImp > >
{
public:
  typedef FunctionSpaceImp FunctionSpaceType;

private:
  typedef ExactSolution< FunctionSpaceType >      ThisType;
  typedef Function< FunctionSpaceType, ThisType > BaseType;

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

  const TimeType t_;

public:
  inline explicit ExactSolution(const FunctionSpaceType& fSpace, const TimeType& t)
    : BaseType(fSpace)
      , t_(t)
  {}

  inline explicit ExactSolution(const FunctionSpaceType& fSpace)
    : BaseType(fSpace)
      , t_(1)
  {}

  // in case 'u' has NO time-dependency use the following method:
  inline void evaluate(const DomainType& x,
                       RangeType& y) const {
    y = 0;
  }

  // in case 'u' HAS a time-dependency use the following method:
  // unfortunately GRAPE requires both cases of the method 'evaluate' to be
  // instantiated
  inline void evaluate(const DomainType& x,
                       const TimeType& timedummy,
                       RangeType& y) const {
    evaluate(x, y);
  }
};

// In the following we want to implement the diffusion tensor A^{\epsilon} //Not the direct tensor A!.
template< class FunctionSpaceImp >
class Diffusion
  : public Function< FunctionSpaceImp, Diffusion< FunctionSpaceImp > >
{
public:
  typedef FunctionSpaceImp FunctionSpaceType;

private:
  typedef Diffusion< FunctionSpaceType >          ThisType;
  typedef Function< FunctionSpaceType, ThisType > BaseType;

public:
  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType  RangeType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType  RangeFieldType;

  typedef DomainFieldType TimeType;

  const FunctionSpaceType* functionSpace_;

protected:
  const double epsilon_;
  const double epsilon_est_;   // estimated_epsilon
  const double delta_;   // the edge length of the coputational cell

  const double time_step_size_;   // the time step size

  // name of the diffusion matrix
  const std::string name_;

public:
  // Constructor for Diffusion
  inline explicit Diffusion(const FunctionSpaceType& functionSpace, const double& epsilon)
    : BaseType(functionSpace)
      , functionSpace_(&functionSpace)
      , epsilon_(epsilon)
      , epsilon_est_(epsilon)
      , delta_(epsilon)
      , time_step_size_(1)
      , name_("no_name")
  {}

  inline explicit Diffusion(const FunctionSpaceType& functionSpace, const double& epsilon, const double& time_step_size)
    : BaseType(functionSpace)
      , functionSpace_(&functionSpace)
      , epsilon_(epsilon)
      , epsilon_est_(epsilon)
      , delta_(epsilon)
      , time_step_size_(time_step_size)
      , name_("no_name")
  {}

  inline explicit Diffusion(const FunctionSpaceType& functionSpace,
                            const double& epsilon,
                            const double& time_step_size,
                            const std::string& name)
    : BaseType(functionSpace)
      , functionSpace_(&functionSpace)
      , epsilon_(epsilon)
      , epsilon_est_(epsilon)
      , delta_(epsilon)
      , time_step_size_(time_step_size)
      , name_(name)
  {}

  inline explicit Diffusion(const FunctionSpaceType& functionSpace,
                            const double& epsilon,
                            const double& epsilon_est,     // estimated epsilon
                            const double& time_step_size)
    : BaseType(functionSpace)
      , functionSpace_(&functionSpace)
      , epsilon_(epsilon)
      , epsilon_est_(epsilon_est)
      , delta_(epsilon_est)
      , time_step_size_(time_step_size)
      , name_("no_name")
  {}

  inline explicit Diffusion(const FunctionSpaceType& functionSpace,
                            const double& epsilon,
                            const double& epsilon_est,     // estimated epsilon
                            const double& time_step_size,
                            const std::string& name)
    : BaseType(functionSpace)
      , functionSpace_(&functionSpace)
      , epsilon_(epsilon)
      , epsilon_est_(epsilon_est)
      , delta_(epsilon_est)
      , time_step_size_(time_step_size)
      , name_(name)
  {}

  inline explicit Diffusion(const FunctionSpaceType& functionSpace,
                            const double& epsilon,
                            const double& epsilon_est,     // estimated epsilon
                            const double& delta,
                            const double& time_step_size,
                            const std::string& name)
    : BaseType(functionSpace)
      , functionSpace_(&functionSpace)
      , epsilon_(epsilon)
      , epsilon_est_(epsilon_est)
      , delta_(delta)
      , time_step_size_(time_step_size)
      , name_(name)
  {}

  // t is a fixed point in time for all evaluations

  // evaluate A^{\epsilon}(t,x) = y
  inline void evaluate(const int i,
                       const int j,
                       const DomainType& x,
                       const TimeType& t,
                       RangeType& y) const {
    DomainType x_new;

    x_new[0] = log(fabs(x[0]) + 1) / epsilon_;
    x_new[1] = log(fabs(x[1]) + 1) / epsilon_;

    if (i == j)
    {
      if (i == 0)
      { y = log( 4.0 + ( 2.0 * sin(2.0 * M_PI * x_new[0]) * sin(2.0 * M_PI * x_new[0]) ) ); }
      if (i == 1)
      { y = log( 4.0 + ( 2.0 * cos(2.0 * M_PI * x_new[1]) * cos(2.0 * M_PI * x_new[1]) ) ); }

      y *= 0.01;
    } else {
      y = 0.0;
    }

  } // evaluate

  inline void evaluate(const int i,
                       const int j,
                       const DomainType& x,
                       RangeType& y) const {
    evaluate(i, j, x, 0.0, y);
    std::cout << "WARNING! Call for 'evaluate' method of the Diffusion class without time! Set t=0." << std::endl;
  }

  // dummy implementation
  inline void evaluate(const DomainType& x,
                       RangeType& y) const {
    y = 1.0;
  } // evaluate

  // dummy implementation
  inline void evaluate(const DomainType& x,
                       const TimeType time,
                       RangeType& y) const {
    return evaluate(x, y);
  }

  inline void getTimeStepSize(double& time_step_size) const {
    time_step_size = time_step_size_;
  }

  inline void getEpsilon(double& epsilon) const {
    epsilon = epsilon_;
  }

  inline void getEpsilonEstimated(double& epsilon_est) const {
    epsilon_est = epsilon_est_;
  }

  inline void getDelta(double& delta) const {
    delta = delta_;
  }

  // get an information on whether we use the solutions of cell problems that are already computed and saved in a file
  // with the name 'name_'
  inline std::string getName_and_getBool(bool& use_saved) const {
    if (name_ == "no_name")
    { use_saved = false; } else
    { use_saved = true; }

    return name_;
  } // getName_and_getBool
};

// In the following we want to implement the diffusion tensor A^{\epsilon} //Not the direct tensor A!.
template< class FunctionSpaceImp >
class DriftedDiffusion
  : public Function< FunctionSpaceImp, DriftedDiffusion< FunctionSpaceImp > >
{
public:
  typedef FunctionSpaceImp FunctionSpaceType;

private:
  typedef DriftedDiffusion< FunctionSpaceType >   ThisType;
  typedef Function< FunctionSpaceType, ThisType > BaseType;

public:
  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType  RangeType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType  RangeFieldType;

  typedef DomainFieldType                TimeType;
  typedef Diffusion< FunctionSpaceType > DiffusionType;

  const FunctionSpaceType* functionSpace_;
  const DiffusionType& diffusion_;
  const DomainType& drift_velocity_;

public:
  // Constructor for Diffusion
  inline explicit DriftedDiffusion(const FunctionSpaceType& functionSpace,
                                   const DiffusionType& diffusion,
                                   const DomainType& drift_velocity)
    : BaseType(functionSpace)
      , functionSpace_(&functionSpace)
      , diffusion_(diffusion)
      , drift_velocity_(drift_velocity)
  {}

  // t is a fixed point in time for all evaluations

  // evaluate A^{\epsilon}(t,x) = A(t,\frac{x}{\epsilon}) = y
  inline void evaluate(const int i,
                       const int j,
                       const DomainType& x,
                       const TimeType& t,
                       RangeType& y) const {
    double epsilon = 0.0;

    diffusion_.getEpsilon(epsilon);

    DomainType x_drift(0.0);
    x_drift[0] = x[0] + ( (1.0 / epsilon) * drift_velocity_[0] * t );
    x_drift[1] = x[1] + ( (1.0 / epsilon) * drift_velocity_[1] * t );

    diffusion_.evaluate(i, j, x_drift, t, y);
  } // evaluate

  inline void evaluate(const int i,
                       const int j,
                       const DomainType& x,
                       RangeType& y) const {
    evaluate(i, j, x, 0.0, y);
    std::cout << "WARNING! Call for 'evaluate' method of the Diffusion class without time! Set t=0." << std::endl;
  }

  // dummy implementation
  inline void evaluate(const DomainType& x,
                       RangeType& y) const {
    y = 1.0;
  }

  // dummy implementation
  inline void evaluate(const DomainType& x,
                       const TimeType time,
                       RangeType& y) const {
    return evaluate(x, y);
  }

  inline void getTimeStepSize(double& time_step_size) const {
    diffusion_.getTimeStepSize(time_step_size);
  }

  inline void getEpsilon(double& epsilon) const {
    diffusion_.getEpsilon(epsilon);
  }

  inline void getEpsilonEstimated(double& epsilon_est) const {
    diffusion_.getEpsilonEstimated(epsilon_est);
  }

  // get an information on whether we use the solutions of cell problems that are already computed and saved in a file
  // with the name 'name_'
  inline std::string getName_and_getBool(bool& use_saved) const {
    std::string name = diffusion_.getName_and_getBool(use_saved);

    return name;
  }
};

// since the coefficients are heterogneous, this class is a dummy
template< class FunctionSpaceImp, class DiffusionImp >
class DirectDiffusion
  : public Function< FunctionSpaceImp, DirectDiffusion< FunctionSpaceImp, DiffusionImp > >
{
public:
  typedef FunctionSpaceImp FunctionSpaceType;
  typedef DiffusionImp     DiffusionType;

private:
  typedef DirectDiffusion< FunctionSpaceType, DiffusionType > ThisType;
  typedef Function< FunctionSpaceType, ThisType >             BaseType;

public:
  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType  RangeType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType  RangeFieldType;

  typedef DomainFieldType TimeType;

  const FunctionSpaceType* functionSpace_;

protected:
  const DiffusionType& diffusion_;

public:
  // Constructor for DirectDiffusion
  inline explicit DirectDiffusion(const FunctionSpaceType& functionSpace, const DiffusionType& diffusion)
    : BaseType(functionSpace)
      , functionSpace_(&functionSpace)
      , diffusion_(diffusion)
  {}

  // evaluate A(t,x) = A^{\epsilon}(t,\epsilon} x) = A(t,\epsilon} \frac{x}{\epsilon}) = y
  inline void evaluate(const int i,
                       const int j,
                       const DomainType& x,
                       const TimeType& t,
                       RangeType& y) const {
    double epsilon;
    diffusion_.getEpsilon(epsilon);
    DomainType x_new;
    x_new[0] = epsilon * x[0];
    x_new[1] = epsilon * x[1];

    // no perturbation ( perturbation = 0.0 )
    diffusion_.evaluate(i, j, x_new, t, y);
  } // evaluate

  inline void evaluate(const int i,
                       const int j,
                       const DomainType& x,
                       RangeType& y) const {
    evaluate(i, j, x, 0.0, y);
    std::cout << "WARNING! Call for 'evaluate' method of the DirectDiffusion class without time! Set t=0." << std::endl;
  }

  // dummy implementation
  inline void evaluate(const DomainType& x,
                       RangeType& y) const {
    y = 1.0;
  }

  // dummy implementation
  inline void evaluate(const DomainType& x,
                       const TimeType& t,
                       RangeType& y) const {
    return evaluate(x, t, y);
  }
};

// In the following we want to implement the diffusion tensor A (and therefore A^{\epsilon}.
// Two different kinds of evaluate-methods will be given. One to evaluate A and one to evaluate
// A^{\epsilon}=A(\frac{\cdot}{\epsilon}).
template< class FunctionSpaceImp >
class Advection
  : public Function< FunctionSpaceImp, Advection< FunctionSpaceImp > >
{
public:
  typedef FunctionSpaceImp FunctionSpaceType;

private:
  typedef Advection< FunctionSpaceType >          ThisType;
  typedef Function< FunctionSpaceType, ThisType > BaseType;

public:
  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType  RangeType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType  RangeFieldType;

  typedef DomainFieldType                TimeType;
  typedef Diffusion< FunctionSpaceType > DiffusionType;

  const FunctionSpaceType* functionSpace_;
  const DiffusionType& diffusion_;

public:
  // Constructor for Advection
  inline explicit Advection(const FunctionSpaceType& functionSpace, const DiffusionType& diffusion)
    : BaseType(functionSpace)
      , functionSpace_(&functionSpace)
      , diffusion_(diffusion)
  {}

  inline void getTimeStepSize(double& time_step_size) const {
    diffusion_.getTimeStepSize(time_step_size);
  }

  inline void getEpsilon(double& epsilon) const {
    diffusion_.getEpsilon(epsilon);
  }

  inline void getEpsilonEstimated(double& epsilon_est) const {
    diffusion_.getEpsilonEstimated(epsilon_est);
  }

  inline void getDelta(double& delta) const {
    diffusion_.getDelta(delta);
  }

  inline void getAverage(DomainType& bar_b) const {
    #if true
    DomainType zero(0.0);
    bar_b = zero;
    #endif // if true

    // ! Note: restriction: at the moment, only constant b_average works!
    // for generalization change for instance (in the classes DriftedDiffusion and DriftedAdvection) "bar_b * t" to some
    // "B(t)"!
  } // getAverage

  // evaluate b(t,y)=z
  inline void evaluate(const int i,
                       const DomainType& y,
                       const TimeType& t,
                       RangeType& z) const {
    double epsilon;

    getEpsilon(epsilon);

    if (i == 0)
    { z = 1.5 * sin( 2.0 * M_PI * (y[0] / epsilon) ) * sin( 2.0 * M_PI * (y[1] / epsilon) ); } else
    { z = 1.5 * cos( 2.0 * M_PI * (y[0] / epsilon) ) * cos( 2.0 * M_PI * (y[1] / epsilon) ); }

    z = (0.01 * z) / epsilon;
    // z = 0.0; //! LOESCHEN!!!!!!!!!!
  } // evaluate

  inline void evaluate(const int i,
                       const DomainType& y,
                       RangeType& z) const {
    evaluate(i, y, 0.0, z);
    std::cout << "WARNING! Call for 'evaluate' method of the Advection class without time! Set t=0." << std::endl;
  }

  // dummy implementation
  inline void evaluate(const DomainType& x,
                       RangeType& y) const {
    y = 1.0;
  }

  // dummy implementation
  inline void evaluate(const DomainType& x,
                       const TimeType time,
                       RangeType& y) const {
    return evaluate(x, y);
  }
};


// In the following we want to implement the diffusion tensor A (and therefore A^{\epsilon}.
// Two different kinds of evaluate-methods will be given. One to evaluate A and one to evaluate
// A^{\epsilon}=A(\frac{\cdot}{\epsilon}).
template< class FunctionSpaceImp >
class DriftedAdvection
  : public Function< FunctionSpaceImp, DriftedAdvection< FunctionSpaceImp > >
{
public:
  typedef FunctionSpaceImp FunctionSpaceType;

private:
  typedef DriftedAdvection< FunctionSpaceType >   ThisType;
  typedef Function< FunctionSpaceType, ThisType > BaseType;

public:
  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType  RangeType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType  RangeFieldType;

  typedef DomainFieldType                TimeType;
  typedef Advection< FunctionSpaceType > AdvectionType;

  const FunctionSpaceType* functionSpace_;
  const AdvectionType& advection_;
  const DomainType& drift_velocity_;

public:
  // Constructor for Diffusion
  inline explicit DriftedAdvection(const FunctionSpaceType& functionSpace,
                                   const AdvectionType& advection,
                                   const DomainType& drift_velocity)
    : BaseType(functionSpace)
      , functionSpace_(&functionSpace)
      , advection_(advection)
      , drift_velocity_(drift_velocity)
  {}

  inline void getTimeStepSize(double& time_step_size) const {
    advection_.getTimeStepSize(time_step_size);
  }

  inline void getEpsilon(double& epsilon) const {
    advection_.getEpsilon(epsilon);
  }

  inline void getEpsilonEstimated(double& epsilon_est) const {
    advection_.getEpsilonEstimated(epsilon_est);
  }

  inline void getAverage(DomainType& bar_b) const {
    advection_.getAverage(bar_b);
  }

  // evaluate b(t,y)=z
  inline void evaluate(const int i,
                       const DomainType& y,
                       const TimeType& t,
                       RangeType& z) const {
    double epsilon = 0.0;

    advection_.getEpsilon(epsilon);

    DomainType x_drift(0.0);
    x_drift[0] = y[0] + ( (1.0 / epsilon) * drift_velocity_[0] * t );
    x_drift[1] = y[1] + ( (1.0 / epsilon) * drift_velocity_[1] * t );

    RangeType eps_b_eps_i = 0.0;
    advection_.evaluate(i, x_drift, t, eps_b_eps_i);
    // note: eps_b_eps_i = (1/eps) b^eps[i]

    z = eps_b_eps_i - ( (1.0 / epsilon) * drift_velocity_[i] );
  } // evaluate

  inline void evaluate(const int i,
                       const DomainType& y,
                       RangeType& z) const {
    evaluate(i, y, 0.0, z);
    std::cout << "WARNING! Call for 'evaluate' method of the Advection class without time! Set t=0." << std::endl;
  }

  // dummy implementation
  inline void evaluate(const DomainType& x,
                       RangeType& y) const {
    y = 1.0;
  }

  // dummy implementation
  inline void evaluate(const DomainType& x,
                       const TimeType time,
                       RangeType& y) const {
    return evaluate(x, y);
  }
};

// since the coefficients are heterogneous, this class is a dummy
template< class FunctionSpaceImp, class AdvectionImp >
class DirectAdvection
  : public Function< FunctionSpaceImp, DirectAdvection< FunctionSpaceImp, AdvectionImp > >
{
public:
  typedef FunctionSpaceImp FunctionSpaceType;
  typedef AdvectionImp     AdvectionType;

private:
  typedef DirectAdvection< FunctionSpaceType, AdvectionType > ThisType;
  typedef Function< FunctionSpaceType, ThisType >             BaseType;

public:
  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType  RangeType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType  RangeFieldType;

  typedef DomainFieldType TimeType;

  const FunctionSpaceType* functionSpace_;

protected:
  const AdvectionType& advection_;

public:
  // Constructor for Advection
  inline explicit DirectAdvection(const FunctionSpaceType& functionSpace, const AdvectionType& advection)
    : BaseType(functionSpace)
      , functionSpace_(&functionSpace)
      , advection_(advection)
  {}

  // evaluate b(t,x) = \eps b^{\epsilon}(t,\epsilon} x) = b(t,\epsilon} \frac{x}{\epsilon}) = y
  inline void evaluate(const int i,
                       const DomainType& x,
                       const TimeType& t,
                       RangeType& y) const {
    double epsilon;

    advection_.getEpsilon(epsilon);

    DomainType x_new;
    x_new[0] = epsilon * x[0];
    x_new[1] = epsilon * x[1];

    advection_.evaluate(i, x_new, t, y);

    y = y * epsilon;
  } // evaluate

  inline void evaluate(const int i,
                       const DomainType& y,
                       RangeType& z) const {
    evaluate(i, y, 0.0, z);
    std::cout << "WARNING! Call for 'evaluate' method of the DirectAdvection class without time! Set t=0." << std::endl;
  }

  // dummy implementation
  inline void evaluate(const DomainType& x,
                       RangeType& y) const {
    y = 1.0;
  }

  // dummy implementation
  inline void evaluate(const DomainType& x,
                       const TimeType time,
                       RangeType& y) const {
    return evaluate(x, y);
  }

  inline void getAverage(DomainType& bar_b) const {
    advection_.getAverage(bar_b);
  }
};

// define the mass term:
template< class FunctionSpaceImp >
class MassTerm
  : public Function< FunctionSpaceImp, MassTerm< FunctionSpaceImp > >
{
public:
  typedef FunctionSpaceImp FunctionSpaceType;

private:
  typedef MassTerm< FunctionSpaceType >           ThisType;
  typedef Function< FunctionSpaceType, ThisType > BaseType;

public:
  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType  RangeType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType  RangeFieldType;

  typedef DomainFieldType TimeType;

protected:
  const double epsilon_;
  const double time_step_size_;
  const TimeType t_;

public:
  // Constructor
  inline explicit MassTerm(const FunctionSpaceType& functionSpace)
    : BaseType(functionSpace)
      , epsilon_(1)
      , time_step_size_(1)
      , t_(0.0)
  {}

  inline explicit MassTerm(const FunctionSpaceType& functionSpace,
                           const RangeFieldType epsilon)
    : BaseType(functionSpace)
      , epsilon_(epsilon)
      , time_step_size_(1)
      , t_(0.0)
  {}

  inline void evaluate(const DomainType& x,
                       RangeType& y) const {
    y[0] = 0;
  }

  // dummy implementation
  inline void evaluate(const DomainType& x,
                       const TimeType time,
                       RangeType& y) const {
    return evaluate(x, y);
  }

  inline void getTimeStepSize(double& time_step_size) const {
    time_step_size = time_step_size_;
  }

  inline void getEpsilon(double& epsilon) const {
    epsilon = epsilon_;
  }
};

// ! Implementation of a simple tensor that can be filled with constant entries.
// This implementation is required to construct a tensor with the values that have been determined by the homogenizer
template< class FunctionSpaceImp >
class FillTensor2D
  : public Function< FunctionSpaceImp, FillTensor2D< FunctionSpaceImp > >
{
public:
  typedef FunctionSpaceImp FunctionSpaceType;

private:
  typedef FillTensor2D< FunctionSpaceType >       ThisType;
  typedef Function< FunctionSpaceType, ThisType > BaseType;

public:
  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType  RangeType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType  RangeFieldType;

  typedef DomainFieldType TimeType;

  const double& time_step_size_;

  const RangeType& a_0_0_;
  const RangeType& a_0_1_;
  const RangeType& a_1_0_;
  const RangeType& a_1_1_;

  const FunctionSpaceType& functionSpace_;

public:
  // Constructor for Tensor
  inline explicit FillTensor2D(const FunctionSpaceType& functionSpace,
                               const double& time_step_size,
                               const RangeType& a_0_0,
                               const RangeType& a_0_1,
                               const RangeType& a_1_0,
                               const RangeType& a_1_1)
    : BaseType(functionSpace)
      , time_step_size_(time_step_size)
      , a_0_0_(a_0_0)
      , a_0_1_(a_0_1)
      , a_1_0_(a_1_0)
      , a_1_1_(a_1_1)
      , functionSpace_(functionSpace)
  {}

  // instantiate all possible cases of the evaluate-method:
  inline void evaluate(const int i,
                       const int j,
                       const DomainType& x,
                       RangeType& z) const {
    if (i == j)
    {
      if (i == 0)
      { z = a_0_0_; } else
      { z = a_1_1_; }
    } else {
      if (i == 0)
      { z = a_0_1_; } else
      { z = a_1_0_; }
    }
  } // evaluate

  // the tensor will be called with time t, since it is newly created for any time step, t is just a dummy (the entries
  // change automaticaly)
  inline void evaluate(const int i,
                       const int j,
                       const DomainType& x,
                       const TimeType& t,
                       RangeType& z) const {
    evaluate(i, j, x, z);
  }

  inline void evaluate(const DomainType& x,
                       RangeType& y) const {
    y = 1;
  }

  inline void evaluate(const DomainType& x,
                       const TimeType time,
                       RangeType& y) const {
    return evaluate(x, y);
  }

  inline void getTimeStepSize(double& time_step_size) const {
    time_step_size = time_step_size_;   // 0.03 / 100;
  }
};

// a dummy function class for functions, vectors and matrices
template< class FunctionSpaceImp >
class DefaultDummyFunction
  : public Function< FunctionSpaceImp, DefaultDummyFunction< FunctionSpaceImp > >
{
public:
  typedef FunctionSpaceImp FunctionSpaceType;

private:
  typedef DefaultDummyFunction< FunctionSpaceType > ThisType;
  typedef Function< FunctionSpaceType, ThisType >   BaseType;

public:
  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType  RangeType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType  RangeFieldType;

  typedef DomainFieldType TimeType;

  const FunctionSpaceType* functionSpace_;

protected:
  const double epsilon_;
  const double time_step_size_;

public:
  // Constructor for Dummy
  inline explicit DefaultDummyFunction(const FunctionSpaceType& functionSpace)
    : BaseType(functionSpace)
      , functionSpace_(&functionSpace)
      , epsilon_(0.0)
      , time_step_size_(0.0)
  {}

  inline explicit DefaultDummyFunction(const FunctionSpaceType& functionSpace,
                                       const double& epsilon,
                                       const double& time_step_size)
    : BaseType(functionSpace)
      , functionSpace_(&functionSpace)
      , epsilon_(epsilon)
      , time_step_size_(time_step_size)
  {}

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

template< class DiscreteFunctionImp >
class DiscFuncAdapter
{
public:
  typedef DiscreteFunctionImp DiscreteFunctionType;

  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType
  DiscreteFunctionSpaceType;

  typedef typename DiscreteFunctionType::RangeType
  RangeType;

  typedef typename DiscreteFunctionSpaceType::IteratorType
  IteratorType;

  typedef typename DiscreteFunctionSpaceType::GridType
  GridType;

  typedef typename DiscreteFunctionSpaceType::GridPartType
  GridPartType;

  typedef typename GridType::template Codim< 0 >::Entity
  EntityType;

  typedef typename GridType::template Codim< 0 >::Geometry
  EnGeometryType;

  typedef typename EntityType::ctype
  coordType;

  typedef typename DiscreteFunctionType::LocalFunctionType
  LocalFunctionType;

  enum { dimension = GridType::dimension };
  enum { spacePolOrd = DiscreteFunctionSpaceType::polynomialOrder };

public:
  // discreteFunction = cof * discreteFunction
  // template <class FunctionType>
  static void multScalar(DiscreteFunctionType& discreteFunction, double& scalar) {
    typedef typename DiscreteFunctionType::DofIteratorType DofIteratorType;

    const DofIteratorType end = discreteFunction.dend();
    for (DofIteratorType it = discreteFunction.dbegin(); it != end; ++it)
      *it = (*it) * scalar;
  } // multScalar

  // discreteFunction1 = discreteFunction1 + discreteFunction2
  static void add(DiscreteFunctionType& discreteFunction1,
                  DiscreteFunctionType& discreteFunction2) {
    typedef typename DiscreteFunctionType::DofIteratorType DofIteratorType;

    double value[10000];
    for (int k = 0; k < 10000; ++k)
    {
      value[k] = 0;
    }

    int globalDofNumber = 0;

    const DofIteratorType end2 = discreteFunction2.dend();
    for (DofIteratorType it = discreteFunction2.dbegin(); it != end2; ++it)
    {
      value[globalDofNumber] = (*it);
      globalDofNumber += 1;
    }

    globalDofNumber = 0;
    const DofIteratorType end1 = discreteFunction1.dend();
    for (DofIteratorType it = discreteFunction1.dbegin(); it != end1; ++it)
    {
      *it = value[globalDofNumber] + (*it);
      globalDofNumber += 1;
    }
  } // add

  static void setZero(DiscreteFunctionType& discreteFunction) {
    typedef typename DiscreteFunctionType::DofIteratorType DofIteratorType;

    const DofIteratorType end = discreteFunction.dend();
    for (DofIteratorType it = discreteFunction.dbegin(); it != end; ++it)
    {
      *it = 0;
    }
  } // setZero
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

#endif // ifndef DUNE_PARABOLIC_MODEL_PROBLEM_SPECIFICATION_HH
