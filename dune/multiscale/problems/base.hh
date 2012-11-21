#ifndef DUNE_MS_PROBLEMS_BASE_HH
#define DUNE_MS_PROBLEMS_BASE_HH

#include <string>
#include <dune/multiscale/problems/constants.hh>
#include <dune/stuff/fem/functions/analytical.hh>

namespace Problem {
/**
 *
 * in general we regard problems of the following type:

 * - div ( A^{\epsilon} (x,\nabla u^{\epsilon}) ) + m^{\epsilon} u^{\epsilon} = f - div G

 * Here we have:
 * u^{\epsilon} = exact solution of the problem
 * A^{\epsilon} = diffusion matrix (or tensor), e.g. with the structure A^{\epsilon}(x) = A(x,\frac{x}{\epsilon})
 * m^{\epsilon} = a mass term (or reaction term), e.g. with the structure m^{\epsilon}(x) = m(x,\frac{x}{\epsilon})
 * f = first source term with the structure f = f(x) (=> no micro-scale dependency)
 * G = second source term with the structure G = G(x) (=> no micro-scale dependency).
 * (Note that 'G' is directly implemented! We do not implement '- div G'!)

 * A^{\epsilon} is can be a monotone operator (=> use HMM, the MsFEM is not implemented for nonlinear problems)


 * ! we use the following class names:


 * class ExactSolution -> describes u^{\epsilon}
 * methods:
 *   evaluate  u^{\epsilon}( x )        --> evaluate
 *   evaluate  \grad u^{\epsilon}( x )  --> evaluateJacobian


 * class Diffusion -> describes A^{\epsilon}
 * methods:
 *   evaluate A^{\epsilon}( x , direction )            --> diffusiveFlux
 *   evaluate DA^{\epsilon}( x , position ) direction  --> jacobianDiffusiveFlux


 * class MassTerm -> describes m^{\epsilon}
 * methods:
 *   evaluate m^{\epsilon}( x )         --> evaluate
 *

 * class FirstSource -> describes f
 * methods:
 *   evaluate f( x )                    --> evaluate


 * class SecondSource -> describes G
 * methods:
 *   evaluate G( x )                    --> evaluate


 ! See 'elliptic_problems/example.hh' for details


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
 *
**/
class IModelProblemData
{
protected:
  const Constants constants_;

public:

  static const bool has_exact_solution = false;
  //! Constructor for ModelProblemData
  inline IModelProblemData(const Constants constants)
    : constants_(constants)
  {}
  virtual ~IModelProblemData()
  {}

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
   * @return macroGridName is set to said path
   * \todo paths need to be relative to binary
   */
  virtual std::string getMacroGridFile() const = 0;
};
}

#endif // DUNE_MS_PROBLEMS_BASE_HH
