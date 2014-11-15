#include <config.h>

#include "coarse_scale_operator.hh"

#include <dune/stuff/common/exceptions.hh>
#include <dune/stuff/common/configuration.hh>
#include <dune/stuff/fem/functions/integrals.hh>
#include <dune/stuff/common/profiler.hh>
#include <dune/stuff/grid/walker.hh>
#include <dune/gdt/operators/projections.hh>
#include <dune/gdt/operators/prolongations.hh>
#include <dune/gdt/spaces/constraints.hh>
#include <dune/gdt/functionals/l2.hh>
#include <dune/multiscale/msfem/localproblems/localproblemsolver.hh>
#include <dune/multiscale/msfem/localproblems/localsolutionmanager.hh>
#include <dune/multiscale/msfem/localproblems/localgridlist.hh>
#include <dune/multiscale/msfem/msfem_traits.hh>
#include <dune/multiscale/problems/base.hh>
#include <dune/multiscale/problems/selector.hh>
#include <dune/multiscale/tools/discretefunctionwriter.hh>
#include <dune/multiscale/tools/misc.hh>
#include <dune/multiscale/msfem/coarse_rhs_functional.hh>
#include <sstream>

namespace Dune {
namespace Multiscale {

Stuff::LA::SparsityPatternDefault CoarseScaleOperator::pattern(const CoarseScaleOperator::RangeSpaceType& range_space,
                                                               const CoarseScaleOperator::SourceSpaceType& source_space,
                                                               const CoarseScaleOperator::GridViewType& grid_view) {
  return range_space.compute_volume_pattern(grid_view, source_space);
}

CoarseScaleOperator::CoarseScaleOperator(const CoarseScaleOperator::SourceSpaceType& source_space_in,
                                         LocalGridList& localGridList)
  : OperatorBaseType(global_matrix_, source_space_in)
  , AssemblerBaseType(source_space_in, source_space_in.grid_view().grid().template leafGridView<Interior_Partition>())
  , global_matrix_(coarse_space().mapper().size(), coarse_space().mapper().size(),
                   EllipticOperatorType::pattern(coarse_space()))
  , local_assembler_(local_operator_, localGridList)
  , msfem_rhs_(coarse_space(), "MsFEM right hand side")
  , dirichlet_projection_(coarse_space()) {
  DSC::Profiler::ScopedTiming st("msfem.coarse.assemble");
  msfem_rhs_.vector() *= 0;
  auto interior = coarse_space().grid_view().grid().template leafGridView<Interior_Partition>();
  CoarseRhsFunctional force_functional(msfem_rhs_.vector(), coarse_space(), localGridList, interior);

  const auto& dirichlet = DMP::getDirichletData();
  const auto& boundary_info = Problem::getModelData()->boundaryInfo();
  const auto& neumann = Problem::getNeumannData();

  typedef CommonTraits::InteriorGridViewType InteriorView;
  GDT::Operators::DirichletProjectionLocalizable<InteriorView, Problem::DirichletDataBase,
                                                 CommonTraits::DiscreteFunctionType>
  dirichlet_projection_operator(interior, boundary_info, *dirichlet, dirichlet_projection_);
  GDT::Functionals::L2Face<Problem::NeumannDataBase, CommonTraits::GdtVectorType, CommonTraits::SpaceType, InteriorView>
  neumann_functional(*neumann, msfem_rhs_.vector(), coarse_space(), interior);

  this->add_codim0_assembler(local_assembler_, this->matrix());
  this->add(force_functional);
  this->add(dirichlet_projection_operator, new DSG::ApplyOn::BoundaryEntities<CommonTraits::InteriorGridViewType>());
  this->add(neumann_functional, new DSG::ApplyOn::NeumannIntersections<CommonTraits::InteriorGridViewType>(boundary_info));
  AssemblerBaseType::assemble(true);

  // apply the dirichlet zero constraints to restrict the system to H^1_0
  GDT::Spaces::Constraints::Dirichlet<typename CommonTraits::GridViewType::Intersection, CommonTraits::RangeFieldType>
  dirichlet_constraints(boundary_info, coarse_space().mapper().maxNumDofs(), coarse_space().mapper().maxNumDofs());
  this->add(dirichlet_constraints, global_matrix_ /*, new GDT::ApplyOn::BoundaryEntities< GridViewType >()*/);
  this->add(dirichlet_constraints,
            force_functional.vector() /*, new GDT::ApplyOn::BoundaryEntities< GridViewType >()*/);
  AssemblerBaseType::assemble(false);

  // substract the operators action on the dirichlet values, since we assemble in H^1 but solve in H^1_0
  CommonTraits::GdtVectorType tmp(coarse_space().mapper().size());
  global_matrix_.mv(dirichlet_projection_.vector(), tmp);
  force_functional.vector() -= tmp;
}

void CoarseScaleOperator::assemble() { DUNE_THROW(Dune::InvalidStateException, "nobody should be calling this"); }

void CoarseScaleOperator::apply_inverse(CoarseScaleOperator::CoarseDiscreteFunction& solution) {

  DSC::Profiler::ScopedTiming st("msfem.coarse.solve");

  BOOST_ASSERT_MSG(msfem_rhs_.dofs_valid(), "Coarse scale RHS DOFs need to be valid!");
  DSC_PROFILER.startTiming("msfem.coarse.linearSolver");
  typedef typename BackendChooser<CoarseDiscreteFunctionSpace>::InverseOperatorType Inverse;
  const  Inverse inverse(
      global_matrix_, msfem_rhs_.space().communicator());

  auto options = Inverse::options("bicgstab.amg.ilu0");
  constexpr bool overwrite = true;
  options.set("preconditioner.anisotropy_dim", CommonTraits::world_dim, overwrite);
  options.set("preconditioner.isotropy_dim", CommonTraits::world_dim, overwrite);
  options.set("verbose", "2", overwrite);
  options.set("preconditioner.verbose", "2", overwrite);
  options.set("smoother.verbose", "2", overwrite);
  options.set("post_check_solves_system", "0", overwrite);
  inverse.apply(msfem_rhs_.vector(), solution.vector(), options);

  if (!solution.dofs_valid())
    DUNE_THROW(InvalidStateException, "Degrees of freedom of coarse solution are not valid!");

  solution.vector() += dirichlet_projection_.vector();

  DSC_PROFILER.stopTiming("msfem.coarse.linearSolver");
  DSC_LOG_INFO << "Time to solve coarse MsFEM problem: " << DSC_PROFILER.getTiming("msfem.coarse.linearSolver") << "ms."
               << std::endl;
}

const CoarseScaleOperator::SourceSpaceType& CoarseScaleOperator::coarse_space() const { return test_space(); }

} // namespace Multiscale {
} // namespace Dune {
