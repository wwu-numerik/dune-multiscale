// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "algorithm.hh"

#include <dune/multiscale/hmm/cell_problem_numbering.hh>
#include <dune/multiscale/tools/misc/outputparameter.hh>
#include <dune/multiscale/common/elliptic_homogenizer.hh>
#include <dune/multiscale/common/righthandside_assembler.hh>
#include <dune/multiscale/common/output_traits.hh>
#include <dune/multiscale/common/error_calc.hh>
#include <dune/multiscale/fem/constantdiffusionmatrix.hh>
#include <dune/multiscale/problems/selector.hh>

#include <dune/stuff/common/ranges.hh>
#include <dune/stuff/common/profiler.hh>
#include <dune/stuff/common/logging.hh>
#include <dune/stuff/common/parameter/configcontainer.hh>

#include <dune/fem/misc/l2norm.hh>

#include <string>
#include <fstream>

#include "fem_traits.hh"

namespace {
const std::string seperator_line =
    "---------------------------------------------------------------------------------\n";
}

namespace Dune {
namespace Multiscale {
namespace FEM {

//! set the dirichlet points to zero
template <class DiscreteFunctionType>
void boundaryTreatment(DiscreteFunctionType& rhs) {
  using namespace Dune::Stuff;
  const auto& discreteFunctionSpace = rhs.space();
  static const unsigned int faceCodim = 1;
  for (const auto& entity : discreteFunctionSpace) {
    for (const auto& intersection : DSC::intersectionRange(discreteFunctionSpace.gridPart(), entity)) {
      if (!intersection.boundary())
        continue;
      if (intersection.boundary() && (intersection.boundaryId() != 1))
        continue;

      auto rhsLocal = rhs.localFunction(entity);
      const auto face = intersection.indexInInside();
      for (auto loc_point : DSC::lagrangePointSetRange<faceCodim>(rhs.space(), entity, face))
        rhsLocal[loc_point] = 0;
    }
  }
} // boundaryTreatment

//! set the dirichlet points to the Dirichlet BC
template <class DirichletBC, class DiscreteFunctionType>
void setDirichletValues(DirichletBC& dirichlet_func, DiscreteFunctionType& func) {
  using namespace Dune::Stuff;
  const auto& discreteFunctionSpace = func.space();
  static const unsigned int faceCodim = 1;
  for (const auto& entity : discreteFunctionSpace) {
    for (const auto& intersection : DSC::intersectionRange(discreteFunctionSpace.gridPart(), entity)) {
      if (Dune::Multiscale::Problem::isDirichletBoundary(intersection)) {
        auto funcLocal = func.localFunction(entity);
        const auto face = intersection.indexInInside();
        for (auto loc_point : DSC::lagrangePointSetRange<faceCodim>(func.space(), entity, face)) {
          const auto& global_point =
              entity.geometry().global(discreteFunctionSpace.lagrangePointSet(entity).point(loc_point));
          CommonTraits::RangeType dirichlet_value(0.0);
          dirichlet_func.evaluate(global_point, dirichlet_value);
          funcLocal[loc_point] = dirichlet_value;
        }
      }
    }
  }
} // setDirichletValues

//! write discrete function to a file + VTK Output
void write_discrete_function(typename CommonTraits::DiscreteFunctionType& discrete_solution) {
  // write the final (discrete) solution to a file
  std::string solution_file = (boost::format("fem_solution_refLevel_%d") % DSC_CONFIG_GET("fem.grid_level", 4)).str();
  DiscreteFunctionWriter(solution_file).append(discrete_solution);

  // writing paraview data output
  // general output parameters
  Dune::Multiscale::OutputParameters outputparam;

  // create and initialize output class
  typename OutputTraits::IOTupleType fem_solution_series(&discrete_solution);
  outputparam.set_prefix((boost::format("fem_solution")).str());
  typename OutputTraits::DataOutputType femsol_dataoutput(discrete_solution.space().gridPart().grid(),
                                                          fem_solution_series, outputparam);
  // write data
  if (DSC_CONFIG_GET("problem.linear", true))
    femsol_dataoutput.writeData(1.0 /*dummy*/, "fem-solution");
  else
    femsol_dataoutput.writeData(1.0 /*dummy*/, "fem-newton-solution");

  //! -------------------------- writing data output Exact Solution ------------------------
  if (Problem::getModelData()->hasExactSolution()) {
    auto u_ptr = Dune::Multiscale::Problem::getExactSolution();
    const auto& u = *u_ptr;
    const OutputTraits::DiscreteExactSolutionType discrete_exact_solution("discrete exact solution ", u,
                                                                          discrete_solution.space().gridPart());
    // create and initialize output class
    OutputTraits::ExSolIOTupleType exact_solution_series(&discrete_exact_solution);
    outputparam.set_prefix("exact_solution");
    OutputTraits::ExSolDataOutputType exactsol_dataoutput(discrete_solution.space().gridPart().grid(),
                                                          exact_solution_series, outputparam);
    // write data
    exactsol_dataoutput.writeData(1.0 /*dummy*/, "exact-solution");
    // -------------------------------------------------------
  }
  //! --------------------------------------------------------------------------------------
}

//! \TODO docme
void solve(typename CommonTraits::DiscreteFunctionType& solution,
           const typename CommonTraits::DiscreteFunctionType& dirichlet_extension,
           const typename CommonTraits::DiscreteFunctionSpaceType& finerDiscreteFunctionSpace,
           const typename FEMTraits::EllipticOperatorType& discrete_elliptic_op,
           const typename CommonTraits::LowerOrderTermType& lower_order_term, // lower order term F(x, u(x), grad u(x) )
           const std::string& filename,
           const RightHandSideAssembler<typename CommonTraits::DiscreteFunctionType>& rhsassembler) {
  static const int fem_polorder = 2 * CommonTraits::DiscreteFunctionSpaceType::polynomialOrder + 2;

  //! *************************** Assembling the problem ****************************

  //! (stiffness) matrix
  typename CommonTraits::LinearOperatorType system_matrix("FEM Newton stiffness matrix", finerDiscreteFunctionSpace,
                                                 finerDiscreteFunctionSpace);

  //! right hand side vector
  // right hand side for the finite element method with Newton solver:
  // ( also right hand side for the finer discrete function space )
  typename CommonTraits::DiscreteFunctionType system_rhs("fem newton rhs", finerDiscreteFunctionSpace);
  system_rhs.clear();

  const auto f = Dune::Multiscale::Problem::getFirstSource();
  // assemble right hand side
  const auto diffusion_op = Problem::getDiffusion();
  // Neumann boundary condition
  const auto neumann_bc = Problem::getNeumannBC();

  if (DSC_CONFIG_GET("problem.linear", true)) {
    DSC_LOG_INFO << "Solving linear problem." << std::endl;
    DSC_LOG_INFO << "Solving linear problem with standard FEM and resolution level "
                 << DSC_CONFIG_GET("fem.grid_level", 4) << "." << std::endl;
    DSC_LOG_INFO << "------------------------------------------------------------------------------" << std::endl;

    // to assemble the computational time
    Dune::Timer assembleTimer;

    // assemble the stiffness matrix
    discrete_elliptic_op.assemble_matrix(system_matrix);

    DSC_LOG_INFO << "Time to assemble standard FEM stiffness matrix: " << assembleTimer.elapsed() << "s" << std::endl;

    // assemble right hand side
    rhsassembler.assemble<fem_polorder>(*f, *diffusion_op, dirichlet_extension, *neumann_bc, system_rhs);

    // set Dirichlet Boundary to zero
    boundaryTreatment(system_rhs);

    const typename FEMTraits::InverseOperatorType fem_biCGStab(
        system_matrix, 1e-8, 1e-8, 20000, DSC_CONFIG_GET("global.cgsolver_verbose", false),
        DSC_CONFIG_GET("fem.algebraic_solver", "bcgs"), DSC_CONFIG_GET("preconditioner_type", std::string("sor")));
    fem_biCGStab(system_rhs, solution);

    DSC_LOG_INFO << "---------------------------------------------------------------------------------" << std::endl;
    DSC_LOG_INFO << "Standard FEM problem solved in " << assembleTimer.elapsed() << "s." << std::endl << std::endl
                 << std::endl;
  } else {
    DSC_LOG_INFO << "Solving non-linear problem." << std::endl;
    DSC_LOG_INFO << "Solving nonlinear problem with FEM + Newton-Method. Resolution level of grid = "
                 << DSC_CONFIG_GET("fem.grid_level", 4) << "." << std::endl;
    DSC_LOG_INFO << "---------------------------------------------------------------------------------" << std::endl;

    Dune::Timer assembleTimer;
    //! residual vector
    // current residual
    typename CommonTraits::DiscreteFunctionType residual(filename + "FEM Newton Residual", finerDiscreteFunctionSpace);
    residual.clear();

    typename CommonTraits::RangeType relative_newton_error_finescale = 10000.0;
    typename CommonTraits::RangeType rhs_L2_norm = 10000.0;

    int iteration_step = 1;
    // the Newton step for the FEM reference problem (solved with Newton Method):
    // L2-Norm of residual < tolerance ?
    double tolerance = 1e-06;
    while (relative_newton_error_finescale > tolerance) {
      // (here: solution = solution from the last iteration step)
      DSC_LOG_INFO << "Newton iteration " << iteration_step << ":" << std::endl;
      Dune::Timer stepAssembleTimer;
      // assemble the stiffness matrix
      discrete_elliptic_op.assemble_jacobian_matrix(solution, dirichlet_extension, system_matrix);

      DSC_LOG_INFO << "Time to assemble FEM Newton stiffness matrix for current iteration: "
                   << stepAssembleTimer.elapsed() << "s" << std::endl;

      rhsassembler.assemble_for_Newton_method<fem_polorder>(*f, *diffusion_op, lower_order_term, solution,
                                                            dirichlet_extension, *neumann_bc, system_rhs);

      const Dune::Fem::L2Norm<typename CommonTraits::DiscreteFunctionType::GridPartType> l2norm(system_rhs.gridPart());
      rhs_L2_norm = l2norm.norm(system_rhs);
      if (rhs_L2_norm < 1e-10) {
        // residual solution almost identical to zero: break
        DSC_LOG_DEBUG << "Residual solution almost identical to zero. Therefore: break loop." << std::endl;
        DSC_LOG_DEBUG << "(L^2-Norm of current right hand side = " << rhs_L2_norm << " < 1e-10)" << std::endl;
        break;
      }
      // set Dirichlet Boundary to zero
      boundaryTreatment(system_rhs);

      const typename FEMTraits::InverseOperatorType fem_newton_biCGStab(system_matrix, 1e-8, 1e-8, 20000, true);
      fem_newton_biCGStab(system_rhs, residual);

      if (residual.dofsValid()) {
        solution += residual;
        relative_newton_error_finescale = l2norm.norm(residual);
        relative_newton_error_finescale /= l2norm.norm(solution);

        DSC_LOG_INFO << "Relative L2-Newton Error = " << relative_newton_error_finescale << std::endl;
        // residual solution almost identical to zero: break
        DSC_LOG_INFO << "Relative L2-Newton Error = " << relative_newton_error_finescale << std::endl;
        if (relative_newton_error_finescale <= tolerance) {
          DSC_LOG_INFO << "Since tolerance = " << tolerance << ": break loop." << std::endl;
        }
        residual.clear();
      } else {
        DSC_LOG_INFO << "WARNING! Invalid dofs in 'residual'." << std::endl;
        break;
      }
      iteration_step += 1;
    }

    if (DSC_CONFIG_GET("problem.linear", true)) {
      DSC_LOG_INFO << "Finite Element Problem solved in " << assembleTimer.elapsed() << "s." << std::endl << std::endl
                   << std::endl;
    } else {
      DSC_LOG_INFO << "Nonlinear Finite Element Problem solved with Newton Method in " << assembleTimer.elapsed()
                   << "s." << std::endl << std::endl << std::endl;
    }
    DSC_LOG_INFO << seperator_line;

  } // end 'problem.linear <-> else'

  //! ********************** End of assembling the reference problem ***************************
}

//! outputs Problem info to output stream
void print_info(const CommonTraits::ModelProblemDataType& info, std::ostream& out) {
  // epsilon is specified in the parameter file
  // 'epsilon' in for instance A^{epsilon}(x) = A(x,x/epsilon)
  const double epsilon_ = DSC_CONFIG_GET("problem.epsilon", 1.0f);
  const int refinement_level_ = DSC_CONFIG_GET("fem.grid_level", 4);
  out << "Log-File for Elliptic Model Problem " << Problem::name() << "." << std::endl << std::endl;
  if (DSC_CONFIG_GET("problem.linear", true))
    out << "Problem is declared as being LINEAR." << std::endl;
  else
    out << "Problem is declared as being NONLINEAR." << std::endl;

  if (info.hasExactSolution()) {
    out << "Exact solution is available." << std::endl << std::endl;
  } else {
    out << "Exact solution is not available." << std::endl << std::endl;
  }
  out << "Computations were made for:" << std::endl << std::endl;
  out << "Refinement Level for Grid = " << refinement_level_ << std::endl << std::endl;

  out << "Epsilon = " << epsilon_ << std::endl << std::endl;
}

//! the main FEM computation
void algorithm(typename CommonTraits::GridPointerType& macro_grid_pointer, const std::string filename) {
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
  // lower order term F(x, u(x), grad u(x) )
  const auto lower_order_term = Problem::getLowerOrderTerm();
  // Dirichlet boundary condition
  const auto dirichlet_bc = Problem::getDirichletBC();

  // discrete function that takes the values of the Dirichlet BC on the Dirichlet Boundary nodes
  // and that is zero elsewhere
  typename CommonTraits::DiscreteFunctionType dirichlet_extension("Dirichlet extension", discreteFunctionSpace);
  dirichlet_extension.clear();
  setDirichletValues(*dirichlet_bc, dirichlet_extension);

  //! define the right hand side assembler tool
  // (for linear and non-linear elliptic and parabolic problems, for sources f and/or G )
  RightHandSideAssembler<typename CommonTraits::DiscreteFunctionType> rhsassembler;

  //! define the discrete (elliptic) operator that describes our problem
  // ( effect of the discretized differential operator on a certain discrete function )
  const typename FEMTraits::EllipticOperatorType discrete_elliptic_op(discreteFunctionSpace, *diffusion_op,
                                                                      lower_order_term);

  //! solution vector
  // - By solution, we denote the "discrete solution" determined with FEM in the linear case or FEM-Newton (nonlinear
  // case)
  //    ( if the elliptic problem is linear, the 'solution' is determined without the Newton method )
  // - solution of the finite element method, where we use the Newton method to solve the non-linear system of equations
  typename CommonTraits::DiscreteFunctionType discrete_solution(filename + " FEM(-Newton) Solution",
                                                                discreteFunctionSpace);
  discrete_solution.clear();

  solve(discrete_solution, dirichlet_extension, discreteFunctionSpace, discrete_elliptic_op, *lower_order_term,
        filename, rhsassembler);
  discrete_solution += dirichlet_extension;

  // write FEM solution to a file and produce a VTK output
  write_discrete_function(discrete_solution);

  ErrorCalculator(nullptr, &discrete_solution).print(DSC_LOG_INFO_0);
}

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

  typedef Homogenizer<typename CommonTraits::DiffusionType> HomogenizerType;

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

  // write FEM solution to a file and produce a VTK output
  // ---------------------------------------------------------------------------------

  // write the final (discrete) solution to a file
  std::string solution_file =
      (boost::format("homogenized_solution_macro_refLevel_%d") % DSC_CONFIG_GET("fem.grid_level", 4)).str();
  DiscreteFunctionWriter(solution_file).append(homogenized_solution);

  // writing paraview data output
  // general output parameters
  Dune::Multiscale::OutputParameters outputparam;

  // create and initialize output class
  typename OutputTraits::IOTupleType hom_fem_solution_series(&homogenized_solution);
  outputparam.set_prefix((boost::format("homogenized_solution")).str());
  typename OutputTraits::DataOutputType homfemsol_dataoutput(homogenized_solution.space().gridPart().grid(),
                                                             hom_fem_solution_series, outputparam);
  homfemsol_dataoutput.writeData(1.0 /*dummy*/, "homogenized-solution");

  // ---------------------------------------------------------------------------------
}

} // namespace FEM {
} // namespace Multiscale {
} // namespace Dune {
