// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DiscreteEllipticMsFEMLocalProblem_HH
#define DiscreteEllipticMsFEMLocalProblem_HH

#include <dune/multiscale/common/traits.hh>
#include <dune/multiscale/common/la_backend.hh>
#include <dune/multiscale/msfem/msfem_traits.hh>
#include <dune/stuff/common/parallel/threadmanager.hh>

namespace Dune {
template <class K, int SIZE>
class FieldVector;
} // namespace Dune

namespace Dune {
namespace Multiscale {

struct LocalFunctor;
class LocalGridList;

//! the essential local msfem problem solver class
class LocalProblemSolver {
private:
  friend struct LocalFunctor;
  typedef typename BackendChooser<MsFEMTraits::LocalSpaceType>::InverseOperatorType InverseOperatorType;

  LocalGridList& localgrid_list_;
  const DS::PerThreadValue<CommonTraits::SpaceType> coarse_space_;

public:
  typedef typename BackendChooser<MsFEMTraits::LocalSpaceType>::LinearOperatorType LinearOperatorType;

  /** \brief constructor - with diffusion operator A^{\epsilon}(x)
   * \param localgrid_list cannot be const because Dune::Fem does not provide Gridparts that can be build on a const
   *grid
   **/
  LocalProblemSolver(CommonTraits::SpaceType coarse_space, LocalGridList& localgrid_list);

  /** method for solving and saving the solutions of the local msfem problems
    * for the whole set of macro-entities and for every unit vector e_i
    * ---- method: solve and save the whole set of local msfem problems -----
    * Use the host-grid entities of Level 'computational_level' as computational domains for the subgrid computations
    * **/
  void solve_for_all_cells();

private:
  //! Solve all local MsFEM problems for one coarse entity at once.
  void solve_all_on_single_cell(const MsFEMTraits::CoarseEntityType& coarseCell,
                                MsFEMTraits::LocalSolutionVectorType& allLocalSolutions) const;
}; // end class

} // namespace Multiscale {
} // namespace Dune {

#endif // #ifndef DiscreteEllipticMsFEMLocalProblem_HH
