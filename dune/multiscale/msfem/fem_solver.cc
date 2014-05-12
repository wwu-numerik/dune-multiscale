#include <config.h>

#include <dune/common/exceptions.hh>
#include <dune/common/timer.hh>
#include <dune/fem/function/common/function.hh>
#include <dune/multiscale/common/righthandside_assembler.hh>
#include <dune/multiscale/common/traits.hh>
#include <dune/multiscale/fem/fem_traits.hh>
#include <dune/stuff/common/logging.hh>
#include <dune/stuff/common/profiler.hh>
#include <dune/stuff/functions/norm.hh>
#include <dune/stuff/common/parameter/configcontainer.hh>

#include <dune/gdt/spaces/continuouslagrange.hh>
#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/operators/elliptic.hh>
#include <dune/gdt/functionals/l2.hh>
#include <dune/gdt/spaces/constraints.hh>
#include <dune/gdt/assembler/system.hh>
#include <dune/gdt/operators/projections.hh>

#include <limits>
#include <sstream>
#include <string>

#include "dune/multiscale/common/dirichletconstraints.hh"
#include "fem_solver.hh"

namespace Dune {
namespace Multiscale {

Elliptic_FEM_Solver::Elliptic_FEM_Solver(const CommonTraits::GdtSpaceType &space)
  : space_(space) {}

void Elliptic_FEM_Solver::apply(const CommonTraits::DiffusionType& diffusion,
                                const CommonTraits::SourceType& force,
                                CommonTraits::GdtDiscreteFunctionType &solution) const {
  DSC_LOG_DEBUG << "Solving linear problem with standard FEM\n";
  DSC_PROFILER.startTiming("fem.apply");

  typedef CommonTraits::GridViewType GridViewType;
  typedef typename GridViewType::Intersection IntersectionType;

  const auto& boundary_info = Problem::getModelData()->boundaryInfo();
  const auto& neumann = Problem::getNeumannData();
  const auto& dirichlet = Problem::getDirichletData();

  const auto& space = space_;
  // elliptic operator (type only, for the sparsity pattern)
  typedef GDT::Operators::EllipticCG< CommonTraits::DiffusionType, CommonTraits::GdtMatrixType, CommonTraits::GdtSpaceType > EllipticOperatorType;
  // container
  CommonTraits::GdtMatrixType system_matrix(space.mapper().size(), space.mapper().size(), EllipticOperatorType::pattern(space));
  CommonTraits::GdtVectorType rhs_vector(space.mapper().size());
  CommonTraits::GdtVectorType dirichlet_shift_vector(space.mapper().size());
  auto& solution_vector = solution.vector();
  // left hand side (elliptic operator)
  EllipticOperatorType elliptic_operator(diffusion, system_matrix, space);
  // right hand side
  GDT::Functionals::L2Volume< CommonTraits::SourceType, CommonTraits::GdtVectorType, CommonTraits::GdtSpaceType > force_functional(force, rhs_vector, space);
  GDT::Functionals::L2Face< CommonTraits::NeumannDataType, CommonTraits::GdtVectorType, CommonTraits::GdtSpaceType >
      neumann_functional(*neumann, rhs_vector, space);
  // dirichlet boundary values
  CommonTraits::GdtDiscreteFunctionType dirichlet_projection(space, dirichlet_shift_vector);
  GDT::Operators::DirichletProjectionLocalizable< GridViewType, CommonTraits::DirichletDataType, CommonTraits::GdtDiscreteFunctionType >
      dirichlet_projection_operator(*(space.grid_view()),
                                    boundary_info,
                                    *dirichlet,
                                    dirichlet_projection);
  // now assemble everything in one grid walk
  GDT::SystemAssembler< CommonTraits::GdtSpaceType > system_assembler(space);
  system_assembler.add(elliptic_operator);
  system_assembler.add(force_functional);
  system_assembler.add(neumann_functional, new GDT::ApplyOn::NeumannIntersections< GridViewType >(boundary_info));
  system_assembler.add(dirichlet_projection_operator,
                       new GDT::ApplyOn::BoundaryEntities< GridViewType >());
  system_assembler.assemble();

  // substract the operators action on the dirichlet values, since we assemble in H^1 but solve in H^1_0
  auto tmp = rhs_vector.copy();
  system_matrix.mv(dirichlet_shift_vector, tmp);
  rhs_vector -= tmp;
  // apply the dirichlet zero constraints to restrict the system to H^1_0
  GDT::Constraints::Dirichlet < typename GridViewType::Intersection, CommonTraits::RangeFieldType >
      dirichlet_constraints(boundary_info, space.mapper().maxNumDofs(), space.mapper().maxNumDofs());
  system_assembler.add(dirichlet_constraints, system_matrix, new GDT::ApplyOn::BoundaryEntities< GridViewType >());
  system_assembler.add(dirichlet_constraints, rhs_vector, new GDT::ApplyOn::BoundaryEntities< GridViewType >());
  system_assembler.assemble();

  // solve the system
  const Stuff::LA::Solver< CommonTraits::GdtMatrixType > linear_solver(system_matrix);
  const auto linear_solver_type = linear_solver.options()[0];
  auto linear_solver_options = linear_solver.options(linear_solver_type);
  linear_solver_options.set("max_iter",                 "5000", true);
  linear_solver_options.set("precision",                "1e-8", true);
  linear_solver_options.set("post_check_solves_system", "0",    true);
  linear_solver.apply(rhs_vector, solution_vector, linear_solver_options);
  // add the dirichlet shift to obtain the solution in H^1
  solution_vector += dirichlet_shift_vector;

  DSC_PROFILER.stopTiming("fem.apply");
  DSC_LOG_DEBUG << "Standard FEM problem solved in " << DSC_PROFILER.getTiming("fem.apply") << "ms.\n";
}

} // namespace Multiscale
}
