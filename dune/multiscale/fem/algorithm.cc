#include <config.h>
// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "algorithm.hh"

#include <dune/common/timer.hh>

#include <dune/multiscale/hmm/cell_problem_numbering.hh>
#include <dune/multiscale/tools/misc/outputparameter.hh>
#include <dune/multiscale/common/righthandside_assembler.hh>
#include <dune/multiscale/common/newton_rhs.hh>
#include <dune/multiscale/common/output_traits.hh>
#include <dune/multiscale/common/error_calc.hh>
#include <dune/multiscale/fem/print_info.hh>
#include <dune/multiscale/problems/selector.hh>
#include <dune/multiscale/msfem/fem_solver.hh>
#include <dune/multiscale/tools/discretefunctionwriter.hh>

#include <dune/stuff/common/ranges.hh>
#include <dune/stuff/common/profiler.hh>
#include <dune/stuff/common/logging.hh>
#include <dune/stuff/common/parameter/configcontainer.hh>

#include <dune/stuff/grid/boundaryinfo.hh>
#include <dune/stuff/la/container/eigen.hh>
#include <dune/stuff/la/solver.hh>
#include <dune/stuff/functions/expression.hh>
#include <dune/stuff/functions/combined.hh>

#ifdef HAVE_DUNE_FEM // <- this is a temporary hack bc I am not up to date with dune-fem
# undef HAVE_DUNE_FEM
#endif

#include <dune/gdt/space/continuouslagrange/pdelab.hh>
#include <dune/gdt/localevaluation/elliptic.hh>
#include <dune/gdt/localoperator/codim0.hh>
#include <dune/gdt/localevaluation/product.hh>
#include <dune/gdt/localfunctional/codim0.hh>
#include <dune/gdt/localfunctional/codim1.hh>
#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/operator/projections.hh>
#include <dune/gdt/assembler/local/codim0.hh>
#include <dune/gdt/assembler/local/codim1.hh>
#include <dune/gdt/space/constraints.hh>
#include <dune/gdt/assembler/system.hh>
#include <dune/gdt/operator/products.hh>
#include <dune/gdt/operator/prolongations.hh>

#include <string>
#include <fstream>

#include "fem_traits.hh"

namespace Dune {
namespace Multiscale {
namespace FEM {

//! the main FEM computation
void algorithm(const std::shared_ptr< CommonTraits::GridType >& macro_grid_pointer, const std::string /*filename*/) {
  using namespace Dune;

  const auto problem_data = Problem::getModelData();
  print_info(*problem_data, DSC_LOG_INFO);
  typedef typename CommonTraits::GridType::LevelGridView GridViewType;
  const auto grid_view_ptr
      = std::make_shared< const GridViewType >(macro_grid_pointer->levelGridView(macro_grid_pointer->maxLevel()));

  // this is the detailed version of how to discretize with dune-gdt, for example the following problem
  // any function derived from Stuff::LocalizableFunctionInterface will serve as a data function
  // +==========================================================+
  // |+========================================================+|
  // ||  Testcase ER07: smooth data, nonhomogeneous dirichlet  ||
  // ||  (see page 858 in Epshteyn, Riviere, 2007)             ||
  // |+--------------------------------------------------------+|
  // ||  domain = [0, 1] x [0 , 1]                             ||
  // ||  diffusion = 1                                         ||
  // ||  force     = 64 pi^2 (cos(8 pi x) + cos(8 pi y))       ||
  // ||  dirichlet = cos(8 pi x) + cos(8 pi y)                 ||
  // ||  exact solution = cos(8 pi x) + cos(8 pi y)            ||
  // |+========================================================+|
  // +==========================================================+
  // * some types
  static const unsigned int polOrder = 1;
  typedef GridViewType::ctype DomainFieldType;
  static const unsigned int   dimDomain = GridViewType::dimension;
  typedef double            RangeFieldType;
  static const unsigned int dimRange = 1;
  typedef typename GridViewType::template Codim< 0 >::Entity EntityType;
  typedef Stuff::Function::Expression < EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange > FunctionType;
  typedef GDT::ContinuousLagrangeSpace::PdelabWrapper< GridViewType, polOrder, RangeFieldType, dimRange > SpaceType;
  typedef Stuff::LA::IstlRowMajorSparseMatrix< RangeFieldType > MatrixType; // <- change those to eigen, if available
  typedef Stuff::LA::IstlDenseVector< RangeFieldType > VectorType;          // <- change those to eigen, if available
  typedef GDT::DiscreteFunction< SpaceType, VectorType > DiscreteFunctionType;
  typedef Stuff::GridboundaryAllDirichlet< typename GridViewType::Intersection > BoundaryInfoType;
  const BoundaryInfoType boundary_info;
  // * create the space and the containers we need
  Dune::Timer timer;
  DSC_LOG_INFO << "assembling system... " << std::flush;
  const SpaceType space(grid_view_ptr);
  const std::unique_ptr< Stuff::LA::SparsityPatternDefault > sparsity_pattern(space.computePattern());
  MatrixType system_matrix(space.mapper().size(), space.mapper().size(), *sparsity_pattern);
  VectorType rhs_vector(space.mapper().size());
  VectorType dirichlet_shift_vector(space.mapper().size());
  // left hand side
  // * elliptic diffusion operator
  const FunctionType diffusion("x", "1.0", 0, "diffusion");
  typedef GDT::LocalOperator::Codim0Integral< GDT::LocalEvaluation::Elliptic< FunctionType > > EllipticOperatorType;
  const EllipticOperatorType diffusion_operator(diffusion);
  // right hand side
  // * L2 force functional
  const FunctionType force("x", "64.0 * pi * pi * (cos(8.0 * pi * x[0]) + cos(8.0 * pi * x[1]))", 3, "force");
  typedef GDT::LocalFunctional::Codim0Integral< GDT::LocalEvaluation::Product< FunctionType > > L2VolumeFunctionalType;
  const L2VolumeFunctionalType force_functional(force);
  // * L2 neumann functional
  const FunctionType neumann("x", "0", 1, "neumann");
  typedef GDT::LocalFunctional::Codim1Integral< GDT::LocalEvaluation::Product< FunctionType > > L2FaceFunctionalType;
  const L2FaceFunctionalType neumann_functional(neumann);
  // dirichlet boundary values
  const FunctionType dirichlet("x", "cos(8.0 * pi * x[0]) + cos(8.0 * pi * x[1])", 3, "dirichlet");
  DiscreteFunctionType dirichlet_projection(space, dirichlet_shift_vector, "discrete dirichlet");
  typedef GDT::ProjectionOperator::Dirichlet< GridViewType > DirichletProjectionOperatorType;
  const DirichletProjectionOperatorType dirichlet_projection_operator(*(space.gridView()), boundary_info);
  dirichlet_projection_operator.apply(dirichlet, dirichlet_projection);
  // local matrix assemBoundaryInfoTypebler
  typedef GDT::LocalAssembler::Codim0Matrix< EllipticOperatorType > LocalEllipticOperatorMatrixAssemblerType;
  const LocalEllipticOperatorMatrixAssemblerType diffusion_matrix_assembler(diffusion_operator);
  // local vector assemblers
  // * force vector
  typedef GDT::LocalAssembler::Codim0Vector< L2VolumeFunctionalType > LocalL2VolumeFunctionalVectorAssemblerType;
  const LocalL2VolumeFunctionalVectorAssemblerType force_vector_assembler(force_functional);
  // * neumann vector
  typedef GDT::LocalAssembler::Codim1Vector< L2FaceFunctionalType > LocalL2FaceFunctionalVectorAssemblerType;
  const LocalL2FaceFunctionalVectorAssemblerType neumann_vector_assembler(neumann_functional);
  // system assembler
  typedef GDT::SystemAssembler< SpaceType > SystemAssemblerType;
  SystemAssemblerType system_assembler(space);
  system_assembler.addLocalAssembler(diffusion_matrix_assembler, system_matrix);
  system_assembler.addLocalAssembler(force_vector_assembler, rhs_vector);
  system_assembler.addLocalAssembler(neumann_vector_assembler,
                                     typename SystemAssemblerType::AssembleOnNeumann(boundary_info),
                                     rhs_vector);
  system_assembler.assemble();
  DSC_LOG_INFO << "done (took " << timer.elapsed() << "s)" << std::endl;
  // substract the operators action on the dirichlet values, since we solve in H^1_0
  // (this is the version that acts on the container interface, if you would only use eigen you can use the backend()
  //  and get rid of the tmp)
  DSC_LOG_INFO << "applying dirichlet constraints... " << std::flush;
  auto tmp = rhs_vector.copy();
  system_matrix.mv(dirichlet_shift_vector, tmp);
  rhs_vector -= tmp;
  // apply the dirichlet zero constraints
  GDT::Constraints::Dirichlet < typename GridViewType::Intersection, RangeFieldType >
    dirichlet_constraints(boundary_info, space.mapper().maxNumDofs(), space.mapper().maxNumDofs());
  system_assembler.addLocalConstraints(dirichlet_constraints, system_matrix);
  system_assembler.addLocalConstraints(dirichlet_constraints, rhs_vector);
  system_assembler.applyConstraints();
  DSC_LOG_INFO << "done (took " << timer.elapsed() << "s)" << std::endl;
  // solve the system
  const Stuff::LA::Solver< MatrixType >linear_solver(system_matrix);
  const auto linear_solver_option = linear_solver.options()[0];
  DSC_LOG_INFO << "solving the linear system using '" << linear_solver_option << "'... " << std::flush;
  VectorType solution_vector(space.mapper().size());
  linear_solver.apply(rhs_vector, solution_vector, linear_solver_option);
  // add the dirichlet shift to obtain the solution in H^1
  solution_vector += dirichlet_shift_vector;
  DSC_LOG_INFO << "done (took " << timer.elapsed() << "s)" << std::endl;

  // visualize the solution
  DSC_LOG_INFO << "visualizing... " << std::flush;
  typedef GDT::ConstDiscreteFunction< SpaceType, VectorType > ConstDiscreteFunctionType;
  const ConstDiscreteFunctionType solution(space, solution_vector);
  diffusion.visualize(*grid_view_ptr, "diffusion");
  force.visualize(*grid_view_ptr, "force");
  dirichlet.visualize(*grid_view_ptr, "dirichlet");
  neumann.visualize(*grid_view_ptr, "neumann");
  solution.visualize("elliptic_fem_solution");
  DSC_LOG_INFO << "done (took " << timer.elapsed() << "s)" << std::endl;

  // measure the error against the exact solution
  DSC_LOG_INFO << "refining grid... " << std::flush;
  macro_grid_pointer->globalRefine(2);
  const auto finer_grid_view_ptr
      = std::make_shared< const GridViewType >(macro_grid_pointer->levelGridView(macro_grid_pointer->maxLevel()));
  DSC_LOG_INFO << "done (took " << timer.elapsed() << "s)" << std::endl;
  DSC_LOG_INFO << "prolonging solution... " << std::flush;
  const SpaceType finer_space(finer_grid_view_ptr);
  VectorType finer_solution_vector(finer_space.mapper().size());
  DiscreteFunctionType finer_solution(finer_space, finer_solution_vector);
  const GDT::ProlongationOperator::Generic< GridViewType > prolongation_operator(*finer_grid_view_ptr);
  prolongation_operator.apply(solution, finer_solution);
  DSC_LOG_INFO << "done (took " << timer.elapsed() << "s)" << std::endl;
  const FunctionType exact_solution("x",
                                    "cos(8.0 * pi * x[0]) + cos(8.0 * pi * x[1])",
                                    3,
                                    "exact solution",
                                    {{"-8.0 * pi * sin(8.0 * pi * x[0])",
                                      "-8.0 * pi * sin(8.0 * pi * x[1])"}});
  const Stuff::Function::Difference< FunctionType, DiscreteFunctionType > difference(exact_solution, finer_solution);
  DSC_LOG_INFO << "measuring relative L2 error: " << std::flush;
  const GDT::ProductOperator::L2Generic< GridViewType > l2_product(*finer_grid_view_ptr);
  DSC_LOG_INFO << std::sqrt(l2_product.apply2(difference, difference))
                  / std::sqrt(l2_product.apply2(exact_solution, exact_solution))
               << " (took " << timer.elapsed() << "s)" << std::endl;
  const GDT::ProductOperator::H1SemiGeneric< GridViewType > h1_semi_product(*finer_grid_view_ptr);
  DSC_LOG_INFO << "measuring relative H1 semi error: " << std::flush;
  DSC_LOG_INFO << std::sqrt(h1_semi_product.apply2(difference, difference))
                  / std::sqrt(h1_semi_product.apply2(exact_solution, exact_solution))
               << " (took " << timer.elapsed() << "s)" << std::endl;
}

} // namespace FEM {
} // namespace Multiscale {
} // namespace Dune {
