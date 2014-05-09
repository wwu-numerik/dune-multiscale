#include <config.h>
// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "algorithm.hh"

#include <string>
#include <fstream>
#include <cmath>

#include <dune/common/timer.hh>

#include <dune/multiscale/hmm/cell_problem_numbering.hh>
#include <dune/multiscale/tools/misc/outputparameter.hh>
#include <dune/multiscale/common/righthandside_assembler.hh>
#include <dune/multiscale/common/output_traits.hh>
#include <dune/multiscale/common/error_calc.hh>
#include <dune/multiscale/fem/print_info.hh>
#include <dune/multiscale/problems/selector.hh>
#include <dune/multiscale/msfem/fem_solver.hh>
#include <dune/multiscale/tools/discretefunctionwriter.hh>

#define DUNE_STUFF_FUNCTIONS_DISABLE_CHECKS 1

#include <dune/stuff/common/ranges.hh>
#include <dune/stuff/common/profiler.hh>
#include <dune/stuff/common/logging.hh>
#include <dune/stuff/common/parameter/configcontainer.hh>
#include <dune/stuff/grid/boundaryinfo.hh>
#include <dune/stuff/la/container.hh>
#include <dune/stuff/la/solver.hh>
#include <dune/stuff/functions/combined.hh>
#include <dune/stuff/functions/constant.hh>
#include <dune/stuff/functions/global.hh>
#include <dune/stuff/grid/provider.hh>

#include <dune/gdt/spaces/continuouslagrange.hh>
#include <dune/gdt/operators/elliptic.hh>
#include <dune/gdt/functionals/l2.hh>
#include <dune/gdt/spaces/constraints.hh>

#include <dune/gdt/assembler/system.hh>
#include <dune/gdt/products/l2.hh>
#include <dune/gdt/products/h1.hh>
#include <dune/gdt/operators/projections.hh>
#include <dune/gdt/operators/prolongations.hh>

#include "fem_traits.hh"

namespace Dune {
namespace Multiscale {
namespace FEM {

//! the main FEM computation
void algorithm(const std::shared_ptr< const CommonTraits::GridType >& macro_grid_pointer, const std::string /*filename*/) {

  Dune::Timer timer;
  typedef CommonTraits::GridViewType GridViewType;
  typedef typename GridViewType::Intersection IntersectionType;

  const auto& diffusion = Problem::getDiffusion();
  const auto& force = Problem::getSource();
  const auto& boundary_info = Problem::getModelData()->boundaryInfo();
  const auto& neumann = Problem::getNeumannData();
  const auto& dirichlet = Problem::getDirichletData();

  Stuff::Grid::Providers::ConstDefault< CommonTraits::GridType > grid_provider(macro_grid_pointer);
  const auto grid_level = 0;
  const CommonTraits::GdtSpaceType space = CommonTraits::GdtSpaceProviderType::create(grid_provider, grid_level);
  const auto grid_view = space.grid_view();
  DSC_LOG_INFO << "assembling system (on a grid with " << grid_view->size(0) << " entities)... "
               << std::flush;

  // elliptic operator (type only, for the sparsity pattern)
  typedef GDT::Operators::EllipticCG< CommonTraits::DiffusionType, CommonTraits::GdtMatrixType, CommonTraits::GdtSpaceType > EllipticOperatorType;
  // container
  CommonTraits::GdtMatrixType system_matrix(space.mapper().size(), space.mapper().size(), EllipticOperatorType::pattern(space));
  CommonTraits::GdtVectorType rhs_vector(space.mapper().size());
  CommonTraits::GdtVectorType dirichlet_shift_vector(space.mapper().size());
  CommonTraits::GdtVectorType solution_vector(space.mapper().size());
  // left hand side (elliptic operator)
  EllipticOperatorType elliptic_operator(*diffusion, system_matrix, space);
  // right hand side
  GDT::Functionals::L2Volume< CommonTraits::SourceType, CommonTraits::GdtVectorType, CommonTraits::GdtSpaceType > force_functional(*force, rhs_vector, space);
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
  DSC_LOG_INFO << "done (took " << timer.elapsed() << "s)" << std::endl;
  timer.reset();
  // substract the operators action on the dirichlet values, since we assemble in H^1 but solve in H^1_0
  DSC_LOG_INFO << "applying dirichlet constraints... " << std::flush;
  auto tmp = rhs_vector.copy();
  system_matrix.mv(dirichlet_shift_vector, tmp);
  rhs_vector -= tmp;
  // apply the dirichlet zero constraints to restrict the system to H^1_0
  GDT::Constraints::Dirichlet < typename GridViewType::Intersection, CommonTraits::RangeFieldType >
      dirichlet_constraints(boundary_info, space.mapper().maxNumDofs(), space.mapper().maxNumDofs());
  system_assembler.add(dirichlet_constraints, system_matrix, new GDT::ApplyOn::BoundaryEntities< GridViewType >());
  system_assembler.add(dirichlet_constraints, rhs_vector, new GDT::ApplyOn::BoundaryEntities< GridViewType >());
  system_assembler.assemble();
  DSC_LOG_INFO << "done (took " << timer.elapsed() << "s)" << std::endl;
  timer.reset();
  // solve the system
  const Stuff::LA::Solver< CommonTraits::GdtMatrixType > linear_solver(system_matrix);
  const auto linear_solver_type = linear_solver.options()[0];
  auto linear_solver_options = linear_solver.options(linear_solver_type);
  linear_solver_options.set("max_iter",                 "5000", true);
  linear_solver_options.set("precision",                "1e-8", true);
  linear_solver_options.set("post_check_solves_system", "0",    true);
  DSC_LOG_INFO << "solving the linear system using '" << linear_solver_type << "'... " << std::flush;
  linear_solver.apply(rhs_vector, solution_vector, linear_solver_options);
  // add the dirichlet shift to obtain the solution in H^1
  solution_vector += dirichlet_shift_vector;
  DSC_LOG_INFO << "done (took " << timer.elapsed() << "s)" << std::endl;
  timer.reset();

  CommonTraits::GdtConstDiscreteFunctionType solution(space, solution_vector);
  ErrorCalculator(nullptr, &solution).print(DSC_LOG_INFO);
} // ... algorithm(...)

} // namespace FEM {
} // namespace Multiscale {
} // namespace Dune {
