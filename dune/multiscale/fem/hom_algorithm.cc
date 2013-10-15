#include "hom_algorithm.hh"


#include <dune/multiscale/common/righthandside_assembler.hh>
#include <dune/multiscale/common/error_calc.hh>
#include <dune/multiscale/fem/constantdiffusionmatrix.hh>
#include <dune/multiscale/fem/print_info.hh>
#include <dune/multiscale/fem/fem_traits.hh>
#include <dune/multiscale/problems/selector.hh>
#include <dune/multiscale/common/elliptic_homogenizer.hh>

#include <dune/stuff/common/ranges.hh>
#include <dune/stuff/common/profiler.hh>
#include <dune/stuff/common/logging.hh>
#include <dune/stuff/common/parameter/configcontainer.hh>

#include <dune/fem/misc/l2norm.hh>

namespace Dune {
namespace Multiscale {
namespace FEM {

//! \TODO docme
void algorithm_hom_fem(typename CommonTraits::GridPointerType& macro_grid_pointer, const std::string filename) {
  using namespace Dune;

  const auto problem_data = Problem::getModelData();
  print_info(*problem_data, DSC_LOG_INFO);

  //! ---------------------------- grid parts ----------------------------------------------
  // grid part for the global function space, required for the finite element problem
  typename CommonTraits::GridPartType gridPart(*macro_grid_pointer);
  //! --------------------------------------------------------------------------------------

  //! ------------------------- discrete function spaces -----------------------------------
  // the global-problem function space:
  typename CommonTraits::DiscreteFunctionSpaceType discreteFunctionSpace(gridPart);
  //! --------------------------------------------------------------------------------------

  // defines the matrix A^{\epsilon} in our global problem  - div ( A^{\epsilon}(\nabla u^{\epsilon} ) = f
  const auto diffusion_op = Problem::getDiffusion();
  //! define the right hand side assembler tool
  // (for linear and non-linear elliptic and parabolic problems, for sources f and/or G )
  RightHandSideAssembler<typename CommonTraits::DiscreteFunctionType> rhsassembler;
  const auto f = Problem::getFirstSource(); // standard source f

  //! define the discrete (elliptic) operator that describes our problem
  // ( effect of the discretized differential operator on a certain discrete function )
  const typename FEMTraits::EllipticOperatorType discrete_elliptic_op(discreteFunctionSpace, *diffusion_op);

  // unit cube grid for the computations of cell problems
  const std::string unit_cell_location = "../dune/multiscale/grids/cell_grids/unit_cube.dgf";
  // descretized homogenizer:

  typedef Homogenizer HomogenizerType;

  // to create an empty diffusion matrix that can be filled with constant values
  typedef Dune::Multiscale::ConstantDiffusionMatrix<typename HomogenizerType::HomTensorType> HomDiffusionType;

  const HomogenizerType disc_homogenizer(unit_cell_location);
  const typename HomogenizerType::HomTensorType A_hom = disc_homogenizer.getHomTensor(*diffusion_op);
  const HomDiffusionType hom_diffusion_op(A_hom);

  //!TODO check: hatte nur 2 tmp parameter, Masse/CommonTraits::LowerOrderTermType hinzugefUGT
  typedef DiscreteEllipticOperator<typename CommonTraits::DiscreteFunctionType, HomDiffusionType>
  HomEllipticOperatorType;

  HomEllipticOperatorType hom_discrete_elliptic_op(discreteFunctionSpace, hom_diffusion_op);

  typename CommonTraits::LinearOperatorType hom_stiff_matrix("homogenized stiffness matrix", discreteFunctionSpace,
                                                             discreteFunctionSpace);

  typename CommonTraits::DiscreteFunctionType hom_rhs("homogenized rhs", discreteFunctionSpace);
  hom_rhs.clear();

  //! solution vector
  // - By solution, we denote the (discrete) homogenized solution determined with FEM on the coarse scale and FEM for
  // the cell problems
  typename CommonTraits::DiscreteFunctionType homogenized_solution(filename + " Homogenized Solution",
                                                                   discreteFunctionSpace);
  homogenized_solution.clear();
  hom_discrete_elliptic_op.assemble_matrix(hom_stiff_matrix);

  constexpr int hmm_polorder = 2 * CommonTraits::DiscreteFunctionSpaceType::polynomialOrder + 2;
  rhsassembler.assemble<hmm_polorder>(*f, hom_rhs);

  // set Dirichlet Boundary to zero
  boundaryTreatment(hom_rhs);

  const typename FEMTraits::InverseOperatorType hom_biCGStab(hom_stiff_matrix, 1e-8, 1e-8, 20000,
                                                             DSC_CONFIG_GET("global.cgsolver_verbose", false));
  hom_biCGStab(hom_rhs, homogenized_solution);
  write_discrete_function(homogenized_solution, "hom_fem");
}

} // namespace FEM {
} // namespace Multiscale {
} // namespace Dune {
