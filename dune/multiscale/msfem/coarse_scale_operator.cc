#include <config.h>

#include "coarse_scale_operator.hh"

#include <dune/stuff/common/configuration.hh>
#include <dune/stuff/fem/functions/integrals.hh>
#include <dune/stuff/common/profiler.hh>
#include <dune/gdt/operators/projections.hh>
#include <dune/gdt/operators/prolongations.hh>
#include <dune/gdt/spaces/constraints.hh>
#include <dune/gdt/functionals/l2.hh>
#include <dune/multiscale/msfem/localproblems/localproblemsolver.hh>
#include <dune/multiscale/msfem/localproblems/localsolutionmanager.hh>
#include <dune/multiscale/msfem/localproblems/subgrid-list.hh>
#include <dune/multiscale/msfem/msfem_traits.hh>
#include <dune/multiscale/problems/base.hh>
#include <dune/multiscale/problems/selector.hh>
#include <dune/multiscale/tools/discretefunctionwriter.hh>
#include <dune/multiscale/tools/misc.hh>
#include <sstream>

namespace Dune {
namespace Multiscale {
namespace MsFEM {

Stuff::LA::SparsityPatternDefault CoarseScaleOperator::pattern(const CoarseScaleOperator::RangeSpaceType& range_space,
                                                               const CoarseScaleOperator::SourceSpaceType& source_space,
                                                               const CoarseScaleOperator::GridViewType& grid_view) {
  return range_space.compute_volume_pattern(grid_view, source_space);
}

CoarseScaleOperator::CoarseScaleOperator(const CoarseScaleOperator::SourceSpaceType& src_spc,
                                         LocalGridList& localGridList)
  : OperatorBaseType(global_matrix_, src_spc)
  , AssemblerBaseType(src_spc)
  , global_matrix_(src_spc.mapper().size(), src_spc.mapper().size(), EllipticOperatorType::pattern(src_spc))
  , local_assembler_(local_operator_, localGridList) {
  this->add_codim0_assembler(local_assembler_, this->matrix());
}

void CoarseScaleOperator::assemble() { AssemblerBaseType::assemble(); }

void CoarseScaleOperator::apply_inverse(const CoarseScaleOperator::CoarseDiscreteFunction& rhs,
                                        CoarseScaleOperator::CoarseDiscreteFunction& solution) {
  BOOST_ASSERT_MSG(rhs.dofs_valid(), "Coarse scale RHS DOFs need to be valid!");
  DSC_PROFILER.startTiming("msfem.solveCoarse");
  const typename BackendChooser<CoarseDiscreteFunctionSpace>::InverseOperatorType inverse(
      global_matrix_,
      rhs.space().communicator()); /*, 1e-8, 1e-8, DSC_CONFIG_GET("msfem.solver.iterations", rhs.size()),
DSC_CONFIG_GET("msfem.solver.verbose", false), "bcgs",
DSC_CONFIG_GET("msfem.solver.preconditioner_type", std::string("sor")));*/
  inverse.apply(rhs.vector(), solution.vector());
  if (!solution.dofs_valid())
    DUNE_THROW(InvalidStateException, "Degrees of freedom of coarse solution are not valid!");
  DSC_PROFILER.stopTiming("msfem.solveCoarse");
  DSC_LOG_DEBUG << "Time to solve coarse MsFEM problem: " << DSC_PROFILER.getTiming("msfem.solveCoarse") << "ms."
                << std::endl;
}

CoarseScaleOperator::MatrixType& CoarseScaleOperator::system_matrix() { return global_matrix_; }

const CoarseScaleOperator::MatrixType& CoarseScaleOperator::system_matrix() const { return global_matrix_; }

} // namespace MsFEM {
} // namespace Multiscale {
} // namespace Dune {
