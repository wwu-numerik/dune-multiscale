#include <config.h>

#include <dune/common/exceptions.hh>
#include <dune/common/timer.hh>
#include <dune/multiscale/msfem/coarse_rhs_functional.hh>
#include <dune/multiscale/common/traits.hh>
#include <dune/multiscale/common/grid_creation.hh>
#include <dune/multiscale/msfem/coarse_scale_operator.hh>
#include <dune/xt/common/logging.hh>
#include <dune/xt/common/timings.hh>
#include <dune/xt/common/configuration.hh>

#include <limits>
#include <sstream>
#include <string>

#include "fem_solver.hh"

namespace Dune {
namespace Multiscale {

Elliptic_FEM_Solver::Elliptic_FEM_Solver(const Problem::ProblemContainer& problem, GridPtrType grid)
  : grid_(grid)
  , space_(CommonTraits::SpaceChooserType::PartViewType::create(*grid_, CommonTraits::st_gdt_grid_level))
  , solution_(space_, "fem_solution")
  , problem_(problem)
{
}

Elliptic_FEM_Solver::Elliptic_FEM_Solver(const DMP::ProblemContainer& problem)
  : Elliptic_FEM_Solver(problem, make_fine_grid(problem, nullptr, false))
{
}

CommonTraits::ConstDiscreteFunctionType& Elliptic_FEM_Solver::solve()
{
  apply(solution_);
  return solution_;
}

void Elliptic_FEM_Solver::apply(CommonTraits::DiscreteFunctionType& solution) const
{
  solution.vector() *= 0;
  CoarseScaleOperator elliptic_msfem_op(problem_, space_, nullptr);
  elliptic_msfem_op.apply_inverse(solution);
}

} // namespace Multiscale
}
