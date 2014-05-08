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

#include <dune/gdt/spaces/continuouslagrange/pdelab.hh>
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

template< class GridViewType >
class ProblemNineDiffusion
  : public Stuff::GlobalFunction< typename GridViewType::template Codim< 0 >::Entity
                                , typename GridViewType::ctype, GridViewType::dimension
                                , double, GridViewType::dimension, GridViewType::dimension >
{
  typedef Stuff::GlobalFunction< typename GridViewType::template Codim< 0 >::Entity
                               , typename GridViewType::ctype, GridViewType::dimension
                               , double, GridViewType::dimension, GridViewType::dimension > BaseType;
public:
  typedef typename BaseType::RangeType  RangeType;
  typedef typename BaseType::DomainType DomainType;

  ProblemNineDiffusion() {}

  virtual size_t order() const DS_FINAL DS_OVERRIDE
  {
    return 2;
  }

  virtual void evaluate(const DomainType& xx, RangeType& ret) const DS_FINAL DS_OVERRIDE
  {
    assert(ret.N() == 2);
    assert(ret.M() == 2);
    ret *= 0.0;
    ret[0][0] =
          2.0 * (1.0 / (8.0 * M_PI * M_PI)) * (1.0 / (2.0 + cos(2.0 * M_PI * (xx[0] / 0.05))));
    ret[1][1] = (1.0 / (8.0 * M_PI * M_PI)) * (1.0 + (0.5 * cos(2.0 * M_PI * (xx[0] / 0.05))));
  } // ... evaluate(...)
}; // ... ProblemNineDiffusion(...)

template< class GridViewType >
class ProblemNineForce
  : public Stuff::GlobalFunction< typename GridViewType::template Codim< 0 >::Entity
                                , typename GridViewType::ctype, GridViewType::dimension, double, 1 >
{
  typedef Stuff::GlobalFunction< typename GridViewType::template Codim< 0 >::Entity
                               , typename GridViewType::ctype, GridViewType::dimension, double, 1 > BaseType;
public:
  typedef typename BaseType::RangeType  RangeType;
  typedef typename BaseType::DomainType DomainType;

  ProblemNineForce() {}

  virtual size_t order() const DS_FINAL DS_OVERRIDE
  {
    return 1;
  }

  virtual void evaluate(const DomainType& xx, RangeType& ret) const DS_FINAL DS_OVERRIDE
  {
    const double eps = 0.05;
    const double pi_square = pow(M_PI, 2.0);
    const double x0_eps = (xx[0] / eps);
    const double cos_2_pi_x0_eps = cos(2.0 * M_PI * x0_eps);
    const double sin_2_pi_x0_eps = sin(2.0 * M_PI * x0_eps);
    const double coefficient_0 = 2.0 * (1.0 / (8.0 * M_PI * M_PI)) * (1.0 / (2.0 + cos_2_pi_x0_eps));
    const double coefficient_1 = (1.0 / (8.0 * M_PI * M_PI)) * (1.0 + (0.5 * cos_2_pi_x0_eps));
    const double sin_2_pi_x0 = sin(2.0 * M_PI * xx[0]);
    const double cos_2_pi_x0 = cos(2.0 * M_PI * xx[0]);
    const double sin_2_pi_x1 = sin(2.0 * M_PI * xx[1]);

    const double d_x0_coefficient_0 =
        pow(2.0 + cos_2_pi_x0_eps, -2.0) * (1.0 / (2.0 * M_PI)) * (1.0 / eps) * sin_2_pi_x0_eps;

    const RangeType grad_u =
        (2.0 * M_PI * cos_2_pi_x0 * sin_2_pi_x1) + ((-1.0) * eps * M_PI * (sin_2_pi_x0 * sin_2_pi_x1 * sin_2_pi_x0_eps)) +
        (M_PI * (cos_2_pi_x0 * sin_2_pi_x1 * cos_2_pi_x0_eps));

    const RangeType d_x0_x0_u =
        -(4.0 * pi_square * sin_2_pi_x0 * sin_2_pi_x1) -
        (2.0 * pi_square * (eps + (1.0 / eps)) * cos_2_pi_x0 * sin_2_pi_x1 * sin_2_pi_x0_eps) -
        (4.0 * pi_square * sin_2_pi_x0 * sin_2_pi_x1 * cos_2_pi_x0_eps);

    const RangeType d_x1_x1_u =
        -(4.0 * pi_square * sin_2_pi_x0 * sin_2_pi_x1) -
        (2.0 * pi_square * eps * cos_2_pi_x0 * sin_2_pi_x1 * sin_2_pi_x0_eps);

    ret = -(d_x0_coefficient_0 * grad_u) - (coefficient_0 * d_x0_x0_u) - (coefficient_1 * d_x1_x1_u);
  } // ... evaluate(...)
}; // ... ProblemNineForce(...)


template< class GridViewType >
class ProblemNineExactSolution
  : public Stuff::GlobalFunction< typename GridViewType::template Codim< 0 >::Entity
                                , typename GridViewType::ctype, GridViewType::dimension, double, 1 >
{
  typedef Stuff::GlobalFunction< typename GridViewType::template Codim< 0 >::Entity
                               , typename GridViewType::ctype, GridViewType::dimension, double, 1 > BaseType;
public:
  typedef typename BaseType::RangeType  RangeType;
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::JacobianRangeType JacobianRangeType;

  ProblemNineExactSolution() {}

  virtual size_t order() const DS_FINAL DS_OVERRIDE
  {
    return 3;
  }

  virtual void evaluate(const DomainType& xx, RangeType& ret) const DS_FINAL DS_OVERRIDE
  {
    const double eps = 0.05;
    static const double M_TWOPI = M_PI * 2.0;
    const double x0_eps = (xx[0] / eps);
    const double sin_2_pi_x0_eps = sin(M_TWOPI * x0_eps);
    const double x0_2_pi = M_TWOPI * xx[0];
    const double x1_2_pi = M_TWOPI * xx[1];
    const double sin_2_pi_x0 = sin(x0_2_pi);
    const double cos_2_pi_x0 = cos(x0_2_pi);
    const double sin_2_pi_x1 = sin(x1_2_pi);

    ret = sin_2_pi_x0 * sin_2_pi_x1
        + (0.5 * eps * cos_2_pi_x0 * sin_2_pi_x1 * sin_2_pi_x0_eps);
  } // ... evaluate(...)

  virtual void jacobian(const DomainType& xx, JacobianRangeType& ret) const DS_FINAL DS_OVERRIDE
  {
    const double eps = 0.05;
    static const double M_TWOPI = M_PI * 2.0;
    const double x0_eps = (xx[0] / eps);
    const double cos_2_pi_x0_eps = cos(M_TWOPI * x0_eps);
    const double sin_2_pi_x0_eps = sin(M_TWOPI * x0_eps);
    const double x0_2_pi = M_TWOPI * xx[0];
    const double x1_2_pi = M_TWOPI * xx[1];
    const double sin_2_pi_x0 = sin(x0_2_pi);
    const double cos_2_pi_x0 = cos(x0_2_pi);
    const double sin_2_pi_x1 = sin(x1_2_pi);
    const double cos_2_pi_x1 = cos(x1_2_pi);

    ret[0][1] = (M_TWOPI * sin_2_pi_x0 * cos_2_pi_x1) + (eps * M_PI * cos_2_pi_x0 * cos_2_pi_x1 * sin_2_pi_x0_eps);

    ret[0][0] = (M_TWOPI * cos_2_pi_x0 * sin_2_pi_x1) - (eps * M_PI * sin_2_pi_x0 * sin_2_pi_x1 * sin_2_pi_x0_eps) +
                   (M_PI * cos_2_pi_x0 * sin_2_pi_x1 * cos_2_pi_x0_eps);
  }
}; // ... ProblemNineExactSolution(...)


//! the main FEM computation
void algorithm(const std::shared_ptr< CommonTraits::GridType >& macro_grid_pointer, const std::string /*filename*/) {
  // some types
  typedef typename CommonTraits::GridType::LevelGridView GridViewType;
  typedef typename GridViewType::template Codim< 0 >::Entity EntityType;
  typedef GridViewType::ctype DomainFieldType;
  static const unsigned int   dimDomain = GridViewType::dimension;
  typedef double            RangeFieldType;
  static const unsigned int dimRange = 1;
  const auto grid_view
      = std::make_shared< const GridViewType >(macro_grid_pointer->levelGridView(macro_grid_pointer->maxLevel()));
  typedef DSG::BoundaryInfos::AllDirichlet< typename GridViewType::Intersection > BoundaryInfoType;
  const BoundaryInfoType boundary_info;
  Dune::Timer timer;
  // analytical data (should be problem nine)
  typedef ProblemNineDiffusion< GridViewType > DiffusionType;
  const DiffusionType diffusion;
  typedef ProblemNineForce< GridViewType > ForceType;
  const ForceType force;
  typedef Stuff::Functions::Constant< EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange >
      ConstantFunctionType;
  const ConstantFunctionType neumann(1.0);
  const ConstantFunctionType dirichlet(0.0);
  typedef ProblemNineExactSolution< GridViewType > ExactSolutionType;
  const ExactSolutionType exact_solution;
  DSC_LOG_INFO << "assembling system (on a grid view with " << grid_view->size(0) << " entities)... "
               << std::flush;
  static const unsigned int polOrder = 1;
  typedef GDT::Spaces::ContinuousLagrange::PdelabBased< GridViewType, polOrder, RangeFieldType, dimRange > SpaceType;
  const SpaceType space(grid_view);
  typedef Stuff::LA::IstlRowMajorSparseMatrix< RangeFieldType > MatrixType;
  typedef Stuff::LA::IstlDenseVector< RangeFieldType >          VectorType;
  typedef GDT::DiscreteFunction< SpaceType, VectorType >      DiscreteFunctionType;
  typedef GDT::ConstDiscreteFunction< SpaceType, VectorType > ConstDiscreteFunctionType;
  // elliptic operator (type only, for the sparsity pattern)
  typedef GDT::Operators::EllipticCG< DiffusionType, MatrixType, SpaceType > EllipticOperatorType;
  // container
  MatrixType system_matrix(space.mapper().size(), space.mapper().size(), EllipticOperatorType::pattern(space));
  VectorType rhs_vector(space.mapper().size());
  VectorType dirichlet_shift_vector(space.mapper().size());
  VectorType solution_vector(space.mapper().size());
  // left hand side (elliptic operator)
  EllipticOperatorType elliptic_operator(diffusion, system_matrix, space);
  // right hand side
  GDT::Functionals::L2Volume< ForceType, VectorType, SpaceType > force_functional(force, rhs_vector, space);
  GDT::Functionals::L2Face< ConstantFunctionType, VectorType, SpaceType >
      neumann_functional(neumann, rhs_vector, space);
  // dirichlet boundary values
  DiscreteFunctionType dirichlet_projection(space, dirichlet_shift_vector);
  GDT::Operators::DirichletProjectionLocalizable< GridViewType, ConstantFunctionType, DiscreteFunctionType >
      dirichlet_projection_operator(*(space.grid_view()),
                                    boundary_info,
                                    dirichlet,
                                    dirichlet_projection);
  // now assemble everything in one grid walk
  GDT::SystemAssembler< SpaceType > system_assembler(space);
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
  GDT::Constraints::Dirichlet < typename GridViewType::Intersection, RangeFieldType >
    dirichlet_constraints(boundary_info, space.mapper().maxNumDofs(), space.mapper().maxNumDofs());
  system_assembler.add(dirichlet_constraints, system_matrix, new GDT::ApplyOn::BoundaryEntities< GridViewType >());
  system_assembler.add(dirichlet_constraints, rhs_vector, new GDT::ApplyOn::BoundaryEntities< GridViewType >());
  system_assembler.assemble();
  DSC_LOG_INFO << "done (took " << timer.elapsed() << "s)" << std::endl;
  timer.reset();
  // solve the system
  const Stuff::LA::Solver< MatrixType > linear_solver(system_matrix);
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

  // now we measure the error (there are more intuitive ways to do this in dune-gdt, but this is the most efficient)
  // for the record: the easy way would be
  //   const Products::L2< GridViewType > l2_product(*grid_view);
  //   const double l2_norm = l2_product.induced_norm(exact_solution - solution);
  DSC_LOG_INFO << "computing errors... " << std::flush;
  const ConstDiscreteFunctionType solution(space, solution_vector);
  typedef Stuff::Functions::Difference< ExactSolutionType, ConstDiscreteFunctionType > DifferenceType;
  const DifferenceType difference(exact_solution, solution);
  // therefore we predefine all products
  const size_t over_integrate = 0; // <- would let the product use a higher quadrature oder than needed
  GDT::Products::L2Localizable< GridViewType, DifferenceType >
      l2_error_product(*grid_view, difference, over_integrate);
  GDT::Products::L2Localizable< GridViewType, ExactSolutionType >
      l2_reference_product(*grid_view, exact_solution, over_integrate);
  GDT::Products::H1SemiLocalizable< GridViewType, DifferenceType >
      h1_semi_error_product(*grid_view, difference, over_integrate);
  GDT::Products::H1SemiLocalizable< GridViewType, ExactSolutionType >
      h1_semi_reference_product(*grid_view, exact_solution, over_integrate);
  // so we can apply them all in one grid walk
  system_assembler.add(l2_error_product);
  system_assembler.add(l2_reference_product);
  system_assembler.add(h1_semi_error_product);
  system_assembler.add(h1_semi_reference_product);
  system_assembler.assemble();
  DSC_LOG_INFO << "done (took " << timer.elapsed() << "s)" << std::endl;
  // and access the result in apply2()
  DSC_LOG_INFO << "L2 error      (abs/rel): "
               << std::sqrt(l2_error_product.apply2()) << " / "
               << std::sqrt(l2_error_product.apply2()) / std::sqrt(l2_reference_product.apply2()) << std::endl;
  DSC_LOG_INFO << "semi H1 error (abs/rel): "
               << std::sqrt(h1_semi_error_product.apply2()) << " / "
               << std::sqrt(h1_semi_error_product.apply2()) / std::sqrt(h1_semi_reference_product.apply2())
               << std::endl;
} // ... algorithm(...)

} // namespace FEM {
} // namespace Multiscale {
} // namespace Dune {
