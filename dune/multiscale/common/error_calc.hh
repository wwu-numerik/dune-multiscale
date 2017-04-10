#ifndef ERROR_CALC_HH
#define ERROR_CALC_HH

#include <dune/xt/common/timings.hh>
#include <dune/multiscale/common/traits.hh>
#include <dune/multiscale/problems/selector.hh>
#include <dune/multiscale/msfem/localsolution_proxy.hh>
#include <iosfwd>
#include <memory>

namespace Dune {
namespace Multiscale {

class Elliptic_FEM_Solver;
class LocalsolutionProxy;

class ErrorCalculator
{

public:
  ErrorCalculator(const Problem::ProblemContainer& problem,
                  const std::unique_ptr<LocalsolutionProxy>& msfem_solution,
                  CommonTraits::ConstDiscreteFunctionType* fem_solution);

  //! this one runs cg-fem sim if mandated by config
  ErrorCalculator(const Problem::ProblemContainer& problem, const std::unique_ptr<LocalsolutionProxy>& msfem_solution);

  /** Compute and print errors between exact, msfem and fem solution
   *
   *  This computes the errors between the exact solution, the msfem solution
   *  and the fem solution (if they exist respectively). The results are printed
   *  to the given ostream and returned in a std::map.
   *
   *  \param[in] out The ostream for error printing
   *  \return Returns the computed errors as std::map with the following entries:
   *          - "msfem_exact_L2" The L2 error between the msfem solution and the
   *            exact solution
   *          - "msfem_exact_H1" The H1 error between the msfem solution and the
   *            exact solution
   *          - "fem_exact_L2" The L2 error between the fem solution and the
   *            exact solution
   *          - "fem_exact_H1" The H1 error between the fem solution and the
   *            exact solution
   *          - "msfem_fem_L2" The L2 error between the msfem solution and the
   *            fem solution
   *          Here, each entry will only be present if the corresponding error was
   *          computed.
   */
  std::map<std::string, double> print(std::ostream& out);

private:
  const DMP::ProblemContainer& problem_;
  const std::unique_ptr<LocalsolutionProxy>& msfem_solution_;
  std::unique_ptr<CommonTraits::DiscreteFunctionType> fem_solution_ptr_;
  std::unique_ptr<CommonTraits::DiscreteFunctionType> coarse_fem_solution_ptr_;
  CommonTraits::ConstDiscreteFunctionType* fem_solution_;
  std::unique_ptr<Elliptic_FEM_Solver> fem_solver_;
  Dune::XT::Common::ScopedTiming timing_;
};

} // namespace Multiscale
} // namespace Dune

#endif // ERROR_CALC_HH
