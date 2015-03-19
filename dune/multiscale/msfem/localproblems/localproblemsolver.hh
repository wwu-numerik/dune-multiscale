// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DiscreteEllipticMsFEMLocalProblem_HH
#define DiscreteEllipticMsFEMLocalProblem_HH

#include <dune/common/fmatrix.hh>
#include <dune/common/typetraits.hh>

#include <dune/multiscale/common/traits.hh>
#include <dune/multiscale/msfem/localproblems/localgridlist.hh>
#include <dune/multiscale/tools/misc/outputparameter.hh>
#include <cstddef>
#include <map>
#include <memory>
#include <vector>

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

//! --------------------- the essential local msfem problem solver class ---------------------------
class LocalProblemSolver {

  typedef typename MsFEMTraits::LocalSpaceType LocalSpaceType;
  typedef typename LocalSpaceType::BaseFunctionSetType::RangeType RangeType;
  typedef typename LocalSpaceType::BaseFunctionSetType::JacobianRangeType JacobianRangeType;
  typedef typename LocalSpaceType::DomainType DomainType;

  friend struct LocalFunctor;

public:
  typedef typename BackendChooser<LocalSpaceType>::LinearOperatorType LinearOperatorType;

private:
  typedef typename BackendChooser<LocalSpaceType>::InverseOperatorType InverseOperatorType;

  LocalGridList& localgrid_list_;
  const DS::PerThreadValue<CommonTraits::SpaceType> coarse_space_;

public:
  /** \brief constructor - with diffusion operator A^{\epsilon}(x)
   * \param localgrid_list cannot be const because Dune::Fem does not provide Gridparts that can be build on a const
   *grid
   **/
  LocalProblemSolver(CommonTraits::SpaceType coarse_space, LocalGridList& localgrid_list);

private:
  //! Solve all local MsFEM problems for one coarse entity at once.
  void solve_all_on_single_cell(const MsFEMTraits::CoarseEntityType& coarseCell,
                                MsFEMTraits::LocalSolutionVectorType& allLocalSolutions) const;

public:
  /** method for solving and saving the solutions of the local msfem problems
    * for the whole set of macro-entities and for every unit vector e_i
    * ---- method: solve and save the whole set of local msfem problems -----
    * Use the host-grid entities of Level 'computational_level' as computational domains for the subgrid computations
    * **/
  void solve_for_all_cells();

}; // end class

} // namespace Multiscale {
} // namespace Dune {

#endif // #ifndef DiscreteEllipticMsFEMLocalProblem_HH
