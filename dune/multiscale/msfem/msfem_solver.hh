// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef Elliptic_MSEM_Solver_HH
#define Elliptic_MSEM_Solver_HH

#include <dune/multiscale/common/traits.hh>
#include <dune/multiscale/msfem/msfem_traits.hh>

namespace Dune {
namespace Multiscale {

namespace Problem {
struct ProblemContainer;
}

class LocalGridList;
class LocalsolutionProxy;

//! \TODO needs a better name
class Elliptic_MsFEM_Solver {
private:
  //! identify fine scale part of MsFEM solution (including the projection!)
  void identify_fine_scale_part(const DMP::ProblemContainer& problem, LocalGridList& localgrid_list,
                                const CommonTraits::DiscreteFunctionType& coarse_msfem_solution,
                                const CommonTraits::SpaceType& coarse_space,
                                std::unique_ptr<LocalsolutionProxy>& msfem_solution) const;

public:
  /** - ∇ (A(x,∇u)) + b ∇u + c u = f - divG
   then:
   A --> diffusion operator ('DiffusionOperatorType')
   b --> advective part ('AdvectionTermType')
   c --> reaction part ('ReactionTermType')
   f --> 'first' source term, scalar ('SourceTermType')
   homogenous Dirchilet boundary condition!:
   **/
  void apply(Problem::ProblemContainer& problem, const CommonTraits::SpaceType& coarse_space,
             std::unique_ptr<LocalsolutionProxy>& msfem_solution, LocalGridList& localgrid_list) const;
};

} // namespace Multiscale {
} // namespace Dune {

#endif // #ifndef Elliptic_MSEM_Solver_HH
