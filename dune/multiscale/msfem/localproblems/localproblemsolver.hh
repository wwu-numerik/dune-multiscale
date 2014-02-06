// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DiscreteEllipticMsFEMLocalProblem_HH
#define DiscreteEllipticMsFEMLocalProblem_HH

#include <dune/common/fmatrix.hh>
#include <dune/common/typetraits.hh>
#include <dune/fem/gridpart/common/gridpart.hh>
#include <dune/fem/operator/2order/lagrangematrixsetup.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/operator/matrix/spmatrix.hh>
#include <dune/fem/operator/matrix/spmatrix.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/solver/cginverseoperator.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/matrix.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solver.hh>
#include <dune/multiscale/common/traits.hh>
#include <dune/multiscale/msfem/localproblems/subgrid-list.hh>
#include <dune/multiscale/tools/misc/outputparameter.hh>
#include <cstddef>
#include <map>
#include <memory>
#include <vector>

#include "dune/multiscale/common/la_backend.hh"
#include "dune/multiscale/msfem/msfem_traits.hh"
#include "dune/multiscale/problems/base.hh"

namespace Dune {
template <class K, int SIZE>
class FieldVector;
} // namespace Dune

namespace Dune {
namespace Multiscale {
namespace MsFEM {

//! --------------------- the essential local msfem problem solver class ---------------------------
class LocalProblemSolver {

  typedef typename MsFEMTraits::LocalGridDiscreteFunctionSpaceType LocalGridDiscreteFunctionSpaceType;
  typedef typename LocalGridDiscreteFunctionSpaceType::RangeType RangeType;
  typedef typename LocalGridDiscreteFunctionSpaceType::JacobianRangeType JacobianRangeType;
  typedef typename LocalGridDiscreteFunctionSpaceType::DomainType DomainType;

public:
  typedef typename BackendChooser<LocalGridDiscreteFunctionSpaceType>::LinearOperatorType LinearOperatorType;

private:
  typedef typename BackendChooser<LocalGridDiscreteFunctionSpaceType>::InverseOperatorType InverseOperatorType;

  static std::unique_ptr<InverseOperatorType> make_inverse_operator(LinearOperatorType& problem_matrix);

  const CommonTraits::DiffusionType& diffusion_;
  LocalGridList& subgrid_list_;
  const CommonTraits::DiscreteFunctionSpaceType& coarse_space_;

public:
  /** \brief constructor - with diffusion operator A^{\epsilon}(x)
   * \param subgrid_list cannot be const because Dune::Fem does not provide Gridparts that can be build on a const grid
   **/
  LocalProblemSolver(const CommonTraits::DiscreteFunctionSpaceType& coarse_space, LocalGridList& subgrid_list,
                     const CommonTraits::DiffusionType& diffusion_operator);

private:
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

} // namespace MsFEM {
} // namespace Multiscale {
} // namespace Dune {

#endif // #ifndef DiscreteEllipticMsFEMLocalProblem_HH
