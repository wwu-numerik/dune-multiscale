#include <config.h>

#include <dune/common/exceptions.hh>
#include <dune/common/timer.hh>
#include <dune/multiscale/msfem/coarse_rhs_functional.hh>
#include <dune/multiscale/common/traits.hh>
#include <dune/multiscale/problems/selector.hh>
#include <dune/stuff/common/logging.hh>
#include <dune/stuff/common/profiler.hh>
#include <dune/stuff/common/configuration.hh>

#include <dune/gdt/spaces/continuouslagrange.hh>
#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/operators/elliptic-cg.hh>
#include <dune/gdt/functionals/l2.hh>
#include <dune/gdt/spaces/constraints.hh>
#include <dune/gdt/assembler/system.hh>
#include <dune/gdt/operators/projections.hh>

#include <limits>
#include <sstream>
#include <string>

#include "fem_solver.hh"

namespace Dune {
namespace Multiscale {

Elliptic_FEM_Solver::Elliptic_FEM_Solver(const CommonTraits::GdtSpaceType &space)
  : space_(space) {}

void Elliptic_FEM_Solver::apply(const Problem::DiffusionBase &diffusion,
                                const Problem::SourceType& force,
                                CommonTraits::DiscreteFunctionType &solution) const {
  DSC_LOG_DEBUG_0 << "Solving linear problem with standard FEM" << std::endl;

  DSC_PROFILER.startTiming("fem.apply");

  typedef CommonTraits::GridViewType GridViewType;

  const auto& boundary_info = Problem::getModelData()->boundaryInfo();
  const auto& neumann = Problem::getNeumannData();
  const auto& dirichlet = Problem::getDirichletData();

  const auto& space = space_;

  typedef GDT::Operators::EllipticCG< Problem::DiffusionBase, CommonTraits::LinearOperatorType,
                                      CommonTraits::GdtSpaceType > EllipticOperatorType;
  CommonTraits::LinearOperatorType system_matrix(space.mapper().size(), space.mapper().size(),
                                                 EllipticOperatorType::pattern(space));
  CommonTraits::GdtVectorType rhs_vector(space.mapper().size());
  auto& solution_vector = solution.vector();
  // left hand side (elliptic operator)
  EllipticOperatorType elliptic_operator(diffusion, system_matrix, space);
  // right hand side
  GDT::Functionals::L2Volume< Problem::SourceType, CommonTraits::GdtVectorType, CommonTraits::GdtSpaceType > force_functional(force, rhs_vector, space);
  GDT::Functionals::L2Face< Problem::NeumannDataBase, CommonTraits::GdtVectorType, CommonTraits::GdtSpaceType >
      neumann_functional(*neumann, rhs_vector, space);
  // dirichlet boundary values
  CommonTraits::DiscreteFunctionType dirichlet_projection(space);
  GDT::Operators::DirichletProjectionLocalizable< GridViewType, Problem::DirichletDataBase, CommonTraits::DiscreteFunctionType >
      dirichlet_projection_operator(*(space.grid_view()),
                                    boundary_info,
                                    *dirichlet,
                                    dirichlet_projection);
  DSC_PROFILER.startTiming("fem.assemble");
  // now assemble everything in one grid walk
  GDT::SystemAssembler< CommonTraits::GdtSpaceType > system_assembler(space);
  system_assembler.add(elliptic_operator);
  system_assembler.add(force_functional);
  system_assembler.add(neumann_functional, new GDT::ApplyOn::NeumannIntersections< GridViewType >(boundary_info));
  system_assembler.add(dirichlet_projection_operator,
                       new GDT::ApplyOn::BoundaryEntities< GridViewType >());
  system_assembler.assemble();
  DSC_PROFILER.stopTiming("fem.assemble");

  DSC_PROFILER.startTiming("fem.constraints");
  // substract the operators action on the dirichlet values, since we assemble in H^1 but solve in H^1_0
  CommonTraits::GdtVectorType tmp(space.mapper().size());
  system_matrix.mv(dirichlet_projection.vector(), tmp);
  rhs_vector -= tmp;
  // apply the dirichlet zero constraints to restrict the system to H^1_0
  GDT::Constraints::Dirichlet < typename GridViewType::Intersection, CommonTraits::RangeFieldType >
      dirichlet_constraints(boundary_info, space.mapper().maxNumDofs(), space.mapper().maxNumDofs());
  system_assembler.add(dirichlet_constraints, system_matrix/*, new GDT::ApplyOn::BoundaryEntities< GridViewType >()*/);
  system_assembler.add(dirichlet_constraints, rhs_vector/*, new GDT::ApplyOn::BoundaryEntities< GridViewType >()*/);
  system_assembler.assemble();
  DSC_PROFILER.stopTiming("fem.constraints");

  // solve the system
  DSC_PROFILER.startTiming("fem.solve");
  const Stuff::LA::Solver< CommonTraits::LinearOperatorType, typename CommonTraits::GdtSpaceType::CommunicatorType >
      linear_solver(system_matrix, space.communicator());
  auto linear_solver_options = linear_solver.options("bicgstab.amg.ilu0");
  linear_solver_options.set("max_iter",                 "5000", true);
  linear_solver_options.set("precision",                "1e-8", true);
  linear_solver_options.set("post_check_solves_system", "0",    true);
  linear_solver_options.set("preconditioner.anisotropy_dim", CommonTraits::world_dim, true);
  linear_solver_options.set("preconditioner.isotropy_dim", CommonTraits::world_dim, true);
  linear_solver_options.set("smoother.iterations", "1", true);
  linear_solver_options.set("smoother.relaxation_factor", "0.5", true);
  linear_solver_options.set("criterion.max_level", "100", true);
  linear_solver_options.set("criterion.coarse_target", "1000", true);
  linear_solver_options.set("criterion.min_coarse_rate", "1.2", true);
  linear_solver_options.set("criterion.prolong_damp", "1.6", true);
  linear_solver_options.set("criterion.anisotropy_dim", CommonTraits::world_dim, true);
  linear_solver_options.set("criterion.verbose", "1", true);
  linear_solver_options.set("smoother.verbose", "1", true);


  linear_solver.apply(rhs_vector, solution_vector, linear_solver_options);
  // add the dirichlet shift to obtain the solution in H^1
  DSC_PROFILER.stopTiming("fem.solve");

  solution_vector += dirichlet_projection.vector();

  DSC_PROFILER.stopTiming("fem.apply");
}

} // namespace Multiscale
}
