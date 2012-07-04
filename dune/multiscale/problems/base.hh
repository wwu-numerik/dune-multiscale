#ifndef DUNE_MS_PROBLEMS_BASE_HH
#define DUNE_MS_PROBLEMS_BASE_HH

#include <string>
#include <dune/multiscale/problems/constants.hh>
#include <dune/stuff/functions.hh>

namespace Problem {
/**
 *
 * in general we regard problems of the following type:

 * - div ( A^{\epsilon} (x,\nabla u^{\epsilon}) ) + m^{\epsilon} u^{\epsilon} = f - div G

 * Here we have:
 * u^{\epsilon} = exact solution of the problem
 * A^{\epsilon} = diffusion matrix (or tensor) with the structure A^{\epsilon}(x) = A(x,\frac{x}{\epsilon})
 * m^{\epsilon} = a mass term (or reaction term) with the structure m^{\epsilon}(x) = m(x,\frac{x}{\epsilon})
 * f = first source term with the structure f = f(x) (=> no micro-scale dependency)
 * G = second source term with the structure G = G(x) (=> no micro-scale dependency). Note that 'G' is directly
 * implemented! ( we do not implement '- div G'!)

 * Note, that A^{\epsilon} is a monotone operator

 * ! we define:

 * The entries of the operator A^{\epsilon} by
 * ! a^{\epsilon}_{1}(x_1,x_2) := (**y_1**,**y_2**)
 * ! a^{\epsilon}_{2}(x_1,x_2) := (**y_1**,**y_2**)

 * The mass (or reaction) term m^{\epsilon} is given by:
 * ! m^{\epsilon} := \epsilon
 * Since \epsilon tends to zero, we may say that we do not have a real mass term for our problem. It is a simple
 * condition to fix the solution which is only unique up to a constant. In fact we still approximate the solution of the
 * problem without mass.

 * The first source term f is given by:
 * ! f(x) := ****
 * since the second source is zero, f will form the right hand side (RHS) of our discrete problem

 * The second source term G is constantly zero:
 * ! G(x) := 0

 * !FirstSource defines the right hand side (RHS) of the governing problem (i.e. it defines 'f').
 * The value of the right hand side (i.e. the value of 'f') at 'x' is accessed by the method 'evaluate'. That means 'y
 * := f(x)' and 'y' is returned. It is only important that 'RHSFunction' knows the function space ('FuncSpace') that it
 * is part from. (f \in FunctionSpace)
 *
 */
class IModelProblemData
{
protected:
  //! name of the file where data is saved
  const std::string file_name_;
  const Constants constants_;
  int current_number_of_cell_problem_;

public:

  static const bool has_exact_solution = false;
  //! Constructor for ModelProblemData
  inline IModelProblemData(const Constants constants, const std::string file_name = "no_name")
    : file_name_(file_name)
      , constants_(constants)
      , current_number_of_cell_problem_(-1)
  {}
  virtual ~IModelProblemData()
  {}

  /** epsilon (the smaller epsilon, the finer the micro-structure)
   *  in the periodic setting, epsilon denotes the periode of the fine-scale oscillations
   *  in the non-periodic setting, can be seen as a representative size for the fine-scale behaviour
   **/
  inline double getEpsilon() const {
    return constants_.epsilon;
  }

  //! epsilon (the smaller epsilon, the finer the micro-structure)
  inline double getEpsilonEstimated() const {
    return constants_.epsilon_est;
  }

  /** edge length of a cell (where we solve the cell problems)
   *  we need delta >= epsilon
   **/
  inline double getDelta() const {
    return constants_.delta;
    // NOTE that (delta/epsilon_est) needs to be a positive integer!
  }

  /** get an information on whether we use the solutions of cell problems that are already computed and saved in a file
   * with the name 'name_'
   **/
  inline std::string getName_and_getBool(bool& use_saved) const {
    use_saved = (file_name_ != "no_name");
    return file_name_;
  }

  inline void set_current_number_of_cell_problem(int number) {
    current_number_of_cell_problem_ = number;
  }

  inline int get_current_number_of_cell_problem() const {
    return current_number_of_cell_problem_;
  }

  /**
   * get the (starting) grid refinement level for solving the reference problem
   * in genereal, this is the smallest integer (level), so that solving the reference problem on this level,
   * yields a higly accurate approximation of the exact solution
   * ( here we have heterogenious reference problem, therefore we need a high refinement level )
    * required refinement level for a fine scale reference solution
    * (a saved/precomputed solution is either already available for this level or it must be computed with the
    * following refinement level)
   */
  virtual int  getRefinementLevelReferenceProblem() const = 0;
  /**
   * @brief getMacroGridFile returns a path to a Dune::Grid loadable file (dgf)
   * @param macroGridName is set to said path
   * \todo paths need to be relative to binary
   */
  virtual void getMacroGridFile(std::string& macroGridName) const = 0;
  //! a unique integral identifier for this problem
  virtual int  get_Number_of_Model_Problem() const = 0;
};
}

#endif // DUNE_MS_PROBLEMS_BASE_HH
