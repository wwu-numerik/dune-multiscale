#include <config.h>

#include "coarse_scale_operator.hh"

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/configuration.hh>
#include <dune/xt/common/timings.hh>
#include <dune/xt/grid/walker.hh>
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
#include <dune/multiscale/common/df_io.hh>
#include <dune/multiscale/tools/misc.hh>
#include <dune/multiscale/msfem/coarse_rhs_functional.hh>
#include <dune/xt/common/parallel/partitioner.hh>
#include <dune/grid/utility/partitioning/seedlist.hh>
#include <sstream>

namespace Dune {
namespace Multiscale {

Stuff::LA::SparsityPatternDefault CoarseScaleOperator::pattern(const CoarseScaleOperator::RangeSpaceType& range_space,
                                                               const CoarseScaleOperator::SourceSpaceType& source_space,
                                                               const CoarseScaleOperator::GridViewType& grid_view)
{
  return range_space.compute_volume_pattern(grid_view, source_space);
}

CoarseScaleOperator::CoarseScaleOperator(const DMP::ProblemContainer& problem,
                                         const CoarseScaleOperator::SourceSpaceType& source_space_in,
                                         LocalGridList& localGridList)
  : OperatorBaseType(global_matrix_, source_space_in)
  , AssemblerBaseType(source_space_in, source_space_in.grid_view().grid().leafGridView())
  , global_matrix_(
        coarse_space().mapper().size(), coarse_space().mapper().size(), EllipticOperatorType::pattern(coarse_space()))
  , local_operator_(problem.getDiffusion())
  , local_assembler_(local_operator_, &localGridList)
  , msfem_rhs_(coarse_space(), "MsFEM right hand side")
  , dirichlet_projection_(coarse_space())
  , problem_(problem)
{
  MS_LOG_INFO << "Assembling coarse system" << std::endl;
  Dune::XT::Common::ScopedTiming st("msfem.coarse.assemble");


  const auto& boundary_info = problem_.getModelData().boundaryInfo();
  const auto& neumann = problem_.getNeumannData();
  const auto& dirichlet = problem_.getDirichletData();

  msfem_rhs_.vector() *= 0;
  const auto interior = coarse_space().grid_view().grid().leafGridView();

  //  CoarseRhsFunctional force_functional(problem_, msfem_rhs_.vector(), coarse_space(), localGridList, interior);
  GDT::Functionals::L2Volume<Problem::SourceType, CommonTraits::GdtVectorType, CommonTraits::SpaceType>
      force_functional(problem_.getSource(), msfem_rhs_.vector(), coarse_space());

  GDT::Operators::
      DirichletProjectionLocalizable<UsedViewType, Problem::DirichletDataBase, CommonTraits::DiscreteFunctionType>
          dirichlet_projection_operator(interior, boundary_info, dirichlet, dirichlet_projection_);
  GDT::Functionals::L2Face<Problem::NeumannDataBase, CommonTraits::GdtVectorType, CommonTraits::SpaceType, UsedViewType>
      neumann_functional(neumann, msfem_rhs_.vector(), coarse_space(), interior);

  this->add_codim0_assembler(local_assembler_, this->matrix());
  this->add(force_functional);

  this->add(neumann_functional, new Dune::XT::Grid::ApplyOn::NeumannIntersections<UsedViewType>(boundary_info));
  this->add(dirichlet_projection_operator, new Dune::XT::Grid::ApplyOn::BoundaryEntities<UsedViewType>());
  AssemblerBaseType::assemble(true);

  // substract the operators action on the dirichlet values, since we assemble in H^1 but solve in H^1_0
  CommonTraits::GdtVectorType tmp(coarse_space().mapper().size());
  global_matrix_.mv(dirichlet_projection_.vector(), tmp);
  force_functional.vector() -= tmp;
  // apply the dirichlet zero constraints to restrict the system to H^1_0
  GDT::Spaces::DirichletConstraints<typename UsedViewType::Intersection> dirichlet_constraints(
      boundary_info, coarse_space().mapper().size(), true);
  this->add(dirichlet_constraints, new Dune::XT::Grid::ApplyOn::BoundaryEntities<UsedViewType>());
  //  if (problem.config().get("threading.smp_constraints", false))
  //    AssemblerBaseType::assemble(partitioning);
  //  else
  AssemblerBaseType::assemble(false);
  dirichlet_constraints.apply(global_matrix_, force_functional.vector());
}

void CoarseScaleOperator::assemble()
{
  DUNE_THROW(Dune::InvalidStateException, "nobody should be calling this");
}

void CoarseScaleOperator::apply_inverse(CoarseScaleOperator::CoarseDiscreteFunction& solution)
{
  // to synchronize timing:
  MPIHelper::getCollectiveCommunication().barrier();
  MS_LOG_INFO << "Assembling coarse system took "
              << std::lround(DXTC_TIMINGS.walltime("msfem.coarse.assemble") / 100.) / 10. << "s" << std::endl;
  Dune::XT::Common::ScopedTiming st("msfem.coarse.solve");

  BOOST_ASSERT_MSG(msfem_rhs_.dofs_valid(), "Coarse scale RHS DOFs need to be valid!");
  DXTC_TIMINGS.start("msfem.coarse.linearSolver");
  typedef typename BackendChooser<CoarseDiscreteFunctionSpace>::InverseOperatorType Inverse;
  const Inverse inverse(global_matrix_, coarse_space().communicator());

  const auto type = problem_.config().get("msfem.coarse_solver", "bicgstab.ilut");
  auto options = Inverse::options(type);
  constexpr bool overwrite = true;
  options.set("preconditioner.anisotropy_dim", CommonTraits::world_dim, overwrite);
  options.set("preconditioner.isotropy_dim", CommonTraits::world_dim, overwrite);
  options.set("verbose", problem_.config().get("msfem.coarse_solver.verbose", 2), overwrite);
  options.set("max_iter", problem_.config().get("msfem.coarse_solver.max_iter", 300u), overwrite);
  options.set("preconditioner.verbose", "2", overwrite);
  options.set("smoother.verbose", "2", overwrite);
  options.set("post_check_solves_system", problem_.config().get("msfem.coarse_solver.check", false), overwrite);
  try {
    inverse.apply(msfem_rhs_.vector(), solution.vector(), options);
  } catch (Dune::Stuff::Exceptions::linear_solver_failed& f) {
    // prevents all ranks from outputting the same detailed error message
    MS_LOG_ERROR_0 << f.what();
    DUNE_THROW(InvalidStateException, "Coarse solve failed.");
  }

  if (!solution.dofs_valid())
    DUNE_THROW(InvalidStateException, "Degrees of freedom of coarse solution are not valid!");

  solution.vector() += dirichlet_projection_.vector();

  DXTC_TIMINGS.stop("msfem.coarse.linearSolver");
  MS_LOG_INFO << "Time to solve coarse MsFEM problem: " << DXTC_TIMINGS.walltime("msfem.coarse.linearSolver") << "ms."
              << std::endl;
}

const CoarseScaleOperator::SourceSpaceType& CoarseScaleOperator::coarse_space() const
{
  return test_space();
}

} // namespace Multiscale {
} // namespace Dune {
