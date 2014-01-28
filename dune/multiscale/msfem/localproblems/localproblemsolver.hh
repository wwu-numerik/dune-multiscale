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
#include "dune/multiscale/msfem/msfem_grid_specifier.hh"
#include "dune/multiscale/msfem/msfem_traits.hh"
#include "dune/multiscale/problems/base.hh"

namespace Dune {
template <class K, int SIZE> class FieldVector;
}  // namespace Dune

namespace Dune {
namespace Multiscale {
namespace MsFEM {

/** \brief define output parameters for local problems
 *  appends "local_problems" for path
 **/
struct LocalProblemDataOutputParameters : public OutputParameters {
public:
  explicit LocalProblemDataOutputParameters();
};

//! --------------------- the essential local msfem problem solver class ---------------------------
class LocalProblemSolver {
public:
  typedef typename MsFEMTraits::LocalGridDiscreteFunctionType LocalGridDiscreteFunctionType;

private:
  typedef typename MsFEMTraits::LocalGridType LocalGridType;
  typedef typename MsFEMTraits::LocalGridPartType LocalGridPartType;
  typedef typename MsFEMTraits::LocalGridDiscreteFunctionSpaceType LocalGridDiscreteFunctionSpaceType;

  typedef CommonTraits::DiffusionType DiffusionOperatorType;
  typedef CommonTraits::NeumannBCType NeumannBoundaryType;

  //! @todo LocalGridDiscreteFunctionType should be replaced by some kind of coarse function type
  typedef std::vector<std::shared_ptr<LocalGridDiscreteFunctionType>> CoarseBasisFunctionListType;

  typedef typename LocalGridDiscreteFunctionSpaceType::FunctionSpaceType FunctionSpaceType;
  typedef typename LocalGridDiscreteFunctionSpaceType::RangeType RangeType;
  //! type of value of a gradient of a function
  typedef typename LocalGridDiscreteFunctionSpaceType::JacobianRangeType JacobianRangeType;
  typedef typename LocalGridDiscreteFunctionSpaceType::DomainType DomainType;

  typedef typename LocalGridType::Traits::LeafIndexSet LocalGridLeafIndexSet;
  typedef typename LocalGridType::Traits::GlobalIdSet::IdType IdType;
  typedef typename LocalGridDiscreteFunctionSpaceType::IteratorType LocalGridEntityIteratorType;
  typedef typename LocalGridEntityIteratorType::Entity LocalEntityType;
  typedef typename LocalEntityType::EntityPointer LocalEntityPointerType;
  typedef typename LocalGridType::Codim<0>::Geometry LocalGridEntityGeometry;
  typedef typename LocalGridDiscreteFunctionType::LocalFunctionType HostLocalFunctionType;
  typedef typename LocalGridPartType::IntersectionIteratorType HostIntersectionIterator;

  typedef MsFEMTraits::CoarseEntityType CoarseEntityType;
  typedef MsFEMTraits::LocalSolutionVectorType LocalGridDiscreteFunctionVectorType;

  typedef typename LocalGridDiscreteFunctionSpaceType::IteratorType LocalGridIteratorType;
  typedef typename LocalGridIteratorType::Entity LocalGridEntityType;
  typedef typename LocalGridEntityType::EntityPointer LocalGridEntityPointerType;
  typedef typename LocalGridDiscreteFunctionType::LocalFunctionType SubLocalFunctionType;
  typedef typename LocalGridDiscreteFunctionSpaceType::LagrangePointSetType SGLagrangePointSetType;
  typedef typename LocalGridDiscreteFunctionSpaceType::LagrangePointSetType LocalGridLagrangePointSetType;
  typedef typename LocalGridType::Codim<0>::Geometry SubGridEntityGeometry;

  static const int faceCodim = 1;
  typedef typename LocalGridLagrangePointSetType::Codim<faceCodim>::SubEntityIteratorType LocalGridFaceDofIteratorType;


public:
  typedef typename BackendChooser<LocalGridDiscreteFunctionSpaceType>::LinearOperatorType LocProbLinearOperatorTypeType;

private:
  typedef typename BackendChooser<LocalGridDiscreteFunctionSpaceType>::InverseOperatorType
  InverseLocProbLinearOperatorTypeType;

  static std::unique_ptr<InverseLocProbLinearOperatorTypeType>
  make_inverse_operator(LocProbLinearOperatorTypeType& problem_matrix);

  const DiffusionOperatorType& diffusion_;
  LocalGridList& subgrid_list_;
  const CommonTraits::DiscreteFunctionSpaceType& coarse_space_;

public:
  /** \brief constructor - with diffusion operator A^{\epsilon}(x)
   * \param subgrid_list cannot be const because Dune::Fem does not provide Gridparts that can be build on a const grid
   **/
  LocalProblemSolver(const CommonTraits::DiscreteFunctionSpaceType& coarse_space, LocalGridList& subgrid_list,
                          const DiffusionOperatorType& diffusion_operator);

private:
  void solve_all_on_single_cell(const CoarseEntityType& coarseCell,
                             LocalGridDiscreteFunctionVectorType& allLocalSolutions) const;
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
