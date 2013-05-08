#include "cell_problem_solver.hh"

#include <dune/multiscale/hmm/discrete_cell_operator.hh>

namespace Dune {
namespace Multiscale {
namespace HMM {

const std::string CellProblemSolver::subdir_ = "cell_problems";

void CellProblemSolver::solve_jacobiancorrector_cellproblem(
  const typename CommonTraits::DiscreteFunctionType::DiscreteFunctionSpaceType::JacobianRangeType& gradient_PHI_H,
  const typename CommonTraits::DiscreteFunctionType::DiscreteFunctionSpaceType::JacobianRangeType& grad_old_coarse_function,
  const PeriodicDiscreteFunctionType& corrector_of_old_coarse_function,
  const DomainType& globalQuadPoint,
  PeriodicDiscreteFunctionType& jac_cor_cell_problem_solution) const
{
  // set solution equal to zero:
  jac_cor_cell_problem_solution.clear();

  //! the matrix in our linear system of equations
  // system matrix for the jacobian corrector cell problem to solve:
  // entries:
  // - int_Y DA^{\eps}( x_T + \delta y, \nabla_x u_H^{(n-1)})(x_T) + \nabla_y Q_h( u_H^{(n-1)})(y) ) \nablay_y
  // \phi_h_i(y) \nablay_y \phi_h_j(y) dy
  // ( \phi_h_i and \phi_h_j denote microscopic base functions.)
  CellFEMMatrix jac_cor_cell_system_matrix("Jacobian Corrector Cell Problem System Matrix",
                                           periodicDiscreteFunctionSpace_,
                                           periodicDiscreteFunctionSpace_);

  //! define the discrete (elliptic) cell problem operator
  // ( effect of the discretized differential operator on a certain discrete function )
  DiscreteCellProblemOperator cell_problem_op(periodicDiscreteFunctionSpace_, diffusion_);
  // we are looking for the derivative of the operator in \nabla_x u_H^{(n-1)})(x_T) + \nabla_y Q_h( u_H^{(n-1)})(y)
  // and in direction of the macroscopic base function

  //! right hand side vector of the algebraic jacobian corrector cell problem
  // entries of the right hand side vector:
  // - int_Y DA^{\eps}( x_T + \delta y, \nabla_x u_H^{(n-1)})(x_T) + \nabla_y Q_h( u_H^{(n-1)})(y) )( \nabla_x
  // \Phi_H(x_T) ) \nablay_y \phi_h(y) dy
  // \Phi_H denotes the current macroscopic base function and
  // \phi_h denotes the current microscopic base function.
  PeriodicDiscreteFunctionType jac_cor_cell_problem_rhs("rhs of jacobian corrector cell problem",
                                                        periodicDiscreteFunctionSpace_);
  jac_cor_cell_problem_rhs.clear();

  // assemble the stiffness matrix
  cell_problem_op.assemble_jacobian_matrix(globalQuadPoint,
                                           grad_old_coarse_function,
                                           corrector_of_old_coarse_function,
                                           jac_cor_cell_system_matrix);

  // assemble right hand side of algebraic jacobian corrector cell problem
  cell_problem_op.assemble_jacobian_corrector_cell_prob_RHS
    (globalQuadPoint,
    grad_old_coarse_function,
    corrector_of_old_coarse_function,
    gradient_PHI_H,
    jac_cor_cell_problem_rhs);

  const double norm_rhs = cell_problem_op.normRHS(jac_cor_cell_problem_rhs);

  if ( !( jac_cor_cell_problem_rhs.dofsValid() ) )
  {
    DUNE_THROW(Dune::InvalidStateException, "Jacobian Corrector Cell Problem RHS invalid.");
  }

  // is the right hand side of the jacobian corrector cell problem equal to zero or almost identical to zero?
  // if yes, the solution of the cell problem is also identical to zero. The solver is getting a problem with this
  // situation, which is why we do not solve cell problems for zero-right-hand-side, since we already know the result.
  if (norm_rhs < 1e-10)
  {
    jac_cor_cell_problem_solution.clear();
    // std :: cout << "Jacobian Corrector Cell Problem with solution zero." << std :: endl;
  } else {
    const bool CELLSOLVER_VERBOSE = DSC_CONFIG.get<bool>("problem.cellsolver_verbose", false);
    InverseCellFEMMatrix jac_cor_cell_fem_biCGStab(jac_cor_cell_system_matrix, 1e-8, 1e-8, 20000, CELLSOLVER_VERBOSE);
    jac_cor_cell_fem_biCGStab(jac_cor_cell_problem_rhs, jac_cor_cell_problem_solution);
  }
} // solve_jacobiancorrector_cellproblem

void CellProblemSolver::solvecellproblem(const typename CommonTraits::DiscreteFunctionType::DiscreteFunctionSpaceType::JacobianRangeType& gradient_PHI_H,
                      // the barycenter x_T of a macro grid element 'T'
                      const DomainType& globalQuadPoint,
                      PeriodicDiscreteFunctionType& cell_problem_solution) const {
  // set solution equal to zero:
  cell_problem_solution.clear();

  //! the matrix in our linear system of equations
  // in the non-linear case, it is the matrix for each iteration step
  CellFEMMatrix cell_system_matrix("Cell Problem System Matrix",
                                   periodicDiscreteFunctionSpace_,
                                   periodicDiscreteFunctionSpace_);

  //! define the discrete (elliptic) cell problem operator
  // ( effect of the discretized differential operator on a certain discrete function )
  const DiscreteCellProblemOperator cell_problem_op(periodicDiscreteFunctionSpace_, diffusion_);

  //! right hand side vector of the algebraic cell problem
  // (in the non-linear setting it changes for every iteration step)
  PeriodicDiscreteFunctionType cell_problem_rhs("rhs of cell problem", periodicDiscreteFunctionSpace_);
  cell_problem_rhs.clear();
  const bool CELLSOLVER_VERBOSE = DSC_CONFIG.get<bool>("problem.cellsolver_verbose", false);

  // NOTE:
  // is the right hand side of the cell problem equal to zero or almost identical to zero?
  // if yes, the solution of the cell problem is also identical to zero. The solver is getting a problem with this
  // situation, which is why we do not solve cell problems for zero-right-hand-side, since we already know the result.

  if (DSC_CONFIG_GET("problem.linear", true)) {
    // assemble the stiffness matrix
    cell_problem_op.assemble_matrix(globalQuadPoint, cell_system_matrix);
    // assemble right hand side of algebraic cell problem
    cell_problem_op.assembleCellRHS_linear(globalQuadPoint, gradient_PHI_H, cell_problem_rhs);
    const double norm_rhs = cell_problem_op.normRHS(cell_problem_rhs);
    if ( !( cell_problem_rhs.dofsValid() ) )
    {
      DUNE_THROW(Dune::InvalidStateException,"Cell Problem RHS invalid.");
    }

    if (norm_rhs < /*1e-06*/ 1e-10)
    {
      cell_problem_solution.clear();
      // std :: cout << "Cell problem with solution zero." << std :: endl;
    } else {
      const InverseCellFEMMatrix cell_fem_biCGStab(cell_system_matrix, 1e-8, 1e-8, 20000, CELLSOLVER_VERBOSE);
      cell_fem_biCGStab(cell_problem_rhs, cell_problem_solution);
    }
  } else {
    // nonlinear case:
    // starting value for the Newton method
    PeriodicDiscreteFunctionType zero_func(" constant zero function ", periodicDiscreteFunctionSpace_);
    zero_func.clear();

    // residual vector (current residual)
    PeriodicDiscreteFunctionType cell_problem_residual("Cell Problem Residual", periodicDiscreteFunctionSpace_);
    cell_problem_residual.clear();

    const L2Error< PeriodicDiscreteFunctionType > l2error;
    RangeType relative_newton_error = 10000.0;

    int iteration_step = 0;

    // the Newton step for for solving the current cell problem (solved with Newton Method):
    // L2-Norm of residual < tolerance ?
    const double tolerance = DSC_CONFIG_GET("problem.stochastic_pertubation", false)
                             ? 1e-01 * DSC_CONFIG_GET("problem.stochastic_variance",  0.01)
                             : 1e-06;

    while (relative_newton_error > tolerance)
    {
      // (here: cellproblem_solution = solution from the last iteration step)

      // assemble the stiffness matrix
      cell_problem_op.assemble_jacobian_matrix(globalQuadPoint,
                                               gradient_PHI_H,
                                               cell_problem_solution,
                                               cell_system_matrix);

      // assemble right hand side of algebraic cell problem (for the current iteration step)
      cell_problem_op.assembleCellRHS_nonlinear(globalQuadPoint,
                                                gradient_PHI_H,
                                                cell_problem_solution,
                                                cell_problem_rhs);

      const double norm_rhs = cell_problem_op.normRHS(cell_problem_rhs);
      if ( !( cell_problem_rhs.dofsValid() ) )
      {
        DUNE_THROW(Dune::InvalidStateException, "Cell Problem RHS invalid.");
      }
      if ((norm_rhs < DSC_CONFIG_GET("max_norm_rhs", 1e-10)))
      {
        break;
      }
      double biCG_tolerance = 1e-8;
      bool cell_solution_convenient = false;
      while (!cell_solution_convenient)
      {
        cell_problem_residual.clear();
        const InverseCellFEMMatrix cell_fem_newton_biCGStab(cell_system_matrix,
                                                      1e-8, biCG_tolerance, 20000, CELLSOLVER_VERBOSE);

        cell_fem_newton_biCGStab(cell_problem_rhs, cell_problem_residual);

        if ( cell_problem_residual.dofsValid() )
        { cell_solution_convenient = true; }

        if (biCG_tolerance > 1e-4)
        {
          DSC_LOG_ERROR << "WARNING! Iteration step " << iteration_step
                    << ". Invalid dofs in 'cell_problem_residual', but '" << relative_newton_error
                    <<
          " = relative_newton_error > 1e-01' and 'biCG_tolerance > 1e-4'. L^2-Norm of right hand side of cell problem: "
                    << norm_rhs << ". Therefore possibly inaccurate solution." << std::endl
                    << "Information:" << std::endl
                    << "x_T = globalQuadPoint = " << globalQuadPoint << "." << std::endl
                    << "nabla u_H^{(n-1)} = gradient_PHI_H = " << gradient_PHI_H[0] << "." << std::endl
                    << "Print right hand side? y/n: ";
          char answer;
          std::cin >> answer;
          if ( !(answer == 'n') )
          { cell_problem_op.printCellRHS(cell_problem_rhs); }
          DUNE_THROW(Dune::InvalidStateException, "");
        }
        biCG_tolerance *= 10.0;
      }

      cell_problem_solution += cell_problem_residual;

      relative_newton_error = l2error.template norm2< (2* polynomialOrder) + 2 >(cell_problem_residual, zero_func);
      const RangeType norm_cell_solution = l2error.template norm2< (2* polynomialOrder) + 2 >(cell_problem_solution,
                                                                                        zero_func);
      relative_newton_error = relative_newton_error / norm_cell_solution;
      cell_problem_residual.clear();
      if (iteration_step > 10)
      {
        DSC_LOG_ERROR << "Warning! Algorithm already reached Newton-iteration step " << iteration_step
                      << " for computing the nonlinear cellproblem." << std::endl
                      << "relative_newton_error = " << relative_newton_error << std::endl << std::endl;
        #ifdef FORCE
        residual_L2_norm = 0.0;
        #endif
      }
      iteration_step += 1;
    }
  }
  // end nonlinear case.

  if ( !( cell_problem_solution.dofsValid() ) )
  {
    DUNE_THROW(Dune::InvalidStateException, "Current solution of the cell problem invalid!");
  }
} // solvecellproblem


void CellProblemSolver::saveTheSolutions_baseSet(
  const typename CommonTraits::DiscreteFunctionType::DiscreteFunctionSpaceType& discreteFunctionSpace,
  const CellProblemNumberingManager& cp_num_manager   // just to check, if we use the correct numeration
    ) const
{
  typedef CommonTraits::DiscreteFunctionType DiscreteFunctionType;

  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType
  DiscreteFunctionSpaceType;

  typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;

  typedef typename DiscreteFunctionSpaceType::GridType GridType;

  typedef typename DiscreteFunctionSpaceType::JacobianRangeType
  JacobianRangeType;

  typedef typename DiscreteFunctionSpaceType::BaseFunctionSetType
  BaseFunctionSetType;

  typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;

  typedef typename DiscreteFunctionType::LocalFunctionType
  LocalFunctionType;

  typedef typename GridType::template Codim< 0 >::Entity EntityType;
  typedef typename EntityType::Geometry                  EntityGeometryType;

  typedef CachingQuadrature< GridPartType, 0 > EntityQuadratureType;

  enum { dimension = GridType::dimension };
  enum { maxnumOfBaseFct = 100 };

  std::string cell_solution_location = subdir_ + "/_cellSolutions_baseSet";
  DiscreteFunctionWriter dfw(cell_solution_location);

  DSC_PROFILER.startTiming("hmm.solver.saveTheSolutions_baseSet");

  // we want to determine minimum, average and maxiumum time for solving a cell problem in the current method
  Dune::Stuff::Common::MinMaxAvg<double> cell_time;

  int number_of_cell_problem = 0;

  IteratorType endit = discreteFunctionSpace.end();
  for (IteratorType it = discreteFunctionSpace.begin(); it != endit; ++it)
  {
    // gradients of the macroscopic base functions
    std::vector<JacobianRangeType> gradientPhi(maxnumOfBaseFct);

    // entity
    const EntityType& entity = *it;
    const BaseFunctionSetType baseSet = discreteFunctionSpace.baseFunctionSet(entity);
    const EntityGeometryType& geometry = entity.geometry();
    const EntityQuadratureType quadrature(entity, 0);
    const DomainType barycenter_of_entity = geometry.global( quadrature.point(0) );

    // number of base functions on entity
    const auto numBaseFunctions = baseSet.size();

    // calc Jacobian inverse before volume is evaluated
    const auto& inv = geometry.jacobianInverseTransposed( quadrature.point(0 /*=quadraturePoint*/) );
    baseSet.jacobianAll(quadrature[0], inv, gradientPhi);

    PeriodicDiscreteFunctionType correctorPhi_i("corrector Phi_i", periodicDiscreteFunctionSpace_);
    for (size_t i = 0; i < numBaseFunctions; ++i)
    {
      correctorPhi_i.clear();

      // take time
      DSC_PROFILER.startTiming("none.solvecellproblem");

      solvecellproblem(gradientPhi[i], barycenter_of_entity, correctorPhi_i);

      cell_time(DSC_PROFILER.stopTiming("none.solvecellproblem"));
      DSC_PROFILER.resetTiming("none.solvecellproblem");

      dfw.append(correctorPhi_i);

      // check if we use a correct numeration of the cell problems:
      typename EntityType::EntityPointer entity_pointer(*it);
      if ( !(cp_num_manager.get_number_of_cell_problem(entity_pointer, i) == number_of_cell_problem) )
      {
        DUNE_THROW(Dune::InvalidStateException, "Numeration of cell problems incorrect.");
      }
      number_of_cell_problem++;
    }
  } // end: for-loop: IteratorType it


  const auto total_time = DSC_PROFILER.stopTiming("hmm.solver.saveTheSolutions_baseSet") / 1000.f;
  DSC_LOG_INFO << std::endl;
  DSC_LOG_INFO << "In method: saveTheSolutions_baseSet." << std::endl << std::endl;
  DSC_LOG_INFO << "Cell problems solved for " << discreteFunctionSpace.grid().size(0) << " leaf entities."
                << std::endl;
  DSC_LOG_INFO << "Minimum time for solving a cell problem = " << cell_time.min() << "s." << std::endl;
  DSC_LOG_INFO << "Maximum time for solving a cell problem = " << cell_time.max()  << "s." << std::endl;
  DSC_LOG_INFO << "Average time for solving a cell problem = "
                << cell_time.average() << "s." << std::endl;
  DSC_LOG_INFO << "Total time for computing and saving the cell problems = "
                << total_time << "s," << std::endl << std::endl;
} // saveTheSolutions_baseSet

void CellProblemSolver::saveTheSolutions_discFunc(const CommonTraits::DiscreteFunctionType& macro_discrete_function) const
{
  typedef CommonTraits::DiscreteFunctionType DiscreteFunctionType;

  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType
  DiscreteFunctionSpaceType;

  typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;

  typedef typename DiscreteFunctionSpaceType::GridType GridType;

  typedef typename DiscreteFunctionSpaceType::JacobianRangeType
  JacobianRangeType;

  typedef typename DiscreteFunctionSpaceType::BaseFunctionSetType
  BaseFunctionSetType;

  typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;

  typedef typename DiscreteFunctionType::LocalFunctionType
  LocalFunctionType;

  typedef typename GridType::template Codim< 0 >::Entity EntityType;
  typedef typename EntityType::Geometry                  EntityGeometryType;

  typedef CachingQuadrature< GridPartType, 0 > EntityQuadratureType;

  enum { dimension = GridType::dimension };
  enum { maxnumOfBaseFct = 100 };
  std::string cell_solution_location = subdir_ + "/_cellSolutions_discFunc";
  DiscreteFunctionWriter dfw(cell_solution_location);

  DSC_PROFILER.startTiming("hmm.solver.saveTheSolutions_discFunc");

  // we want to determine minimum, average and maxiumum time for solving a cell problem in the current method
  Dune::Stuff::Common::MinMaxAvg<double> cell_time;

  const DiscreteFunctionSpaceType& discreteFunctionSpace = macro_discrete_function.space();

  int number_of_entity = 0;
  const IteratorType endit = discreteFunctionSpace.end();
  for (IteratorType it = discreteFunctionSpace.begin(); it != endit; ++it)
  {
    // entity
    const EntityType& entity = *it;
    const EntityGeometryType& geometry = entity.geometry();
    const EntityQuadratureType quadrature(entity, 0);
    const DomainType barycenter_of_entity = geometry.global( quadrature.point(0) );

    LocalFunctionType local_macro_disc = macro_discrete_function.localFunction(entity);
    JacobianRangeType grad_macro_discrete_function;
    local_macro_disc.jacobian(quadrature[0], grad_macro_discrete_function);

    PeriodicDiscreteFunctionType cell_solution_on_entity("corrector of macro discrete function",
                                                         periodicDiscreteFunctionSpace_);

    // take time
    DSC_PROFILER.startTiming("hmm.solver.saveTheSolutions_discFunc.solvecellproblem");

    solvecellproblem
      (grad_macro_discrete_function, barycenter_of_entity, cell_solution_on_entity);

    // min/max time
    cell_time(DSC_PROFILER.stopTiming("hmm.solver.saveTheSolutions_discFunc.solvecellproblem"));
    DSC_PROFILER.resetTiming("hmm.solver.saveTheSolutions_discFunc.solvecellproblem");

    dfw.append(cell_solution_on_entity);
    number_of_entity += 1;
  } // end: for-loop: IteratorType it

  const auto total_time = DSC_PROFILER.stopTiming("hmm.solver.saveTheSolutions_discFunc");
  DSC_LOG_INFO << std::endl;
  DSC_LOG_INFO << "In method: saveTheSolutions_discFunc." << std::endl << std::endl;
  DSC_LOG_INFO << "Cell problems solved for " << discreteFunctionSpace.grid().size(0) << " leaf entities."
                << std::endl;
  DSC_LOG_INFO << "Minimum time for solving a cell problem = " << cell_time.min() << "s." << std::endl;
  DSC_LOG_INFO << "Maximum time for solving a cell problem = " << cell_time.max() << "s." << std::endl;
  DSC_LOG_INFO << "Average time for solving a cell problem = "
                << cell_time.average() << "s." << std::endl;
  DSC_LOG_INFO << "Total time for computing and saving the cell problems = "
                << total_time << "s," << std::endl << std::endl;
} // saveTheSolutions_discFunc

void CellProblemSolver::saveTheJacCorSolutions_baseSet_discFunc(
  const CommonTraits::DiscreteFunctionType& macro_discrete_function,
  const CellProblemNumberingManager& cp_num_manager   // just to check, if we use the correct numeration
  ) const
{
  typedef CommonTraits::DiscreteFunctionType DiscreteFunctionType;

  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType
  DiscreteFunctionSpaceType;

  typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;

  typedef typename DiscreteFunctionSpaceType::GridType GridType;

  typedef typename DiscreteFunctionSpaceType::JacobianRangeType
  JacobianRangeType;

  typedef typename DiscreteFunctionSpaceType::BaseFunctionSetType
  BaseFunctionSetType;

  typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;

  typedef typename DiscreteFunctionType::LocalFunctionType
  LocalFunctionType;

  typedef typename GridType::template Codim< 0 >::Entity EntityType;
  typedef typename EntityType::EntityPointer
  EntityPointerType;
  typedef typename EntityType::Geometry EntityGeometryType;

  typedef CachingQuadrature< GridPartType, 0 > EntityQuadratureType;

  enum { dimension = GridType::dimension };
  enum { maxnumOfBaseFct = 100 };
  // where we save the solutions:
  const std::string cell_solution_location = subdir_ + "/_JacCorCellSolutions_baseSet_discFunc";
  DiscreteFunctionWriter dfw(cell_solution_location);
  // where we saved the solutions for the discrete function
  // NOTE: they already need to be assembled, i.e. we already applied the method saveSolutions_discFunc!
  const std::string cell_solution_discFunc_location = subdir_ + "/_cellSolutions_discFunc";

  // reader for the cell problem data file (discrete functions):
  DiscreteFunctionReader discrete_function_reader(cell_solution_discFunc_location);

  DSC_PROFILER.startTiming("hmm.solver.saveTheJacCorSolutions_baseSet_discFunc");

  // we want to determine minimum, average and maxiumum time for solving a cell problem in the current method
  Dune::Stuff::Common::MinMaxAvg<double> cell_time;

  int number_of_cell_problem = 0;

  const DiscreteFunctionSpaceType& discreteFunctionSpace = macro_discrete_function.space();

  int number_of_entity = 0;
  const IteratorType endit = discreteFunctionSpace.end();
  for (IteratorType it = discreteFunctionSpace.begin(); it != endit; ++it)
  {
    // gradients of the macroscopic base functions
    std::vector<JacobianRangeType> gradientPhi(maxnumOfBaseFct);

    // entity
    const EntityType& entity = *it;
    const BaseFunctionSetType baseSet = discreteFunctionSpace.baseFunctionSet(entity);
    const EntityGeometryType& geometry = entity.geometry();
    const EntityQuadratureType quadrature(entity, 0);
    const DomainType barycenter_of_entity = geometry.global( quadrature.point(0) );

    // number of base functions on entity
    const int numBaseFunctions = baseSet.size();

    // calc Jacobian inverse before volume is evaluated
    const auto& inv = geometry.jacobianInverseTransposed( quadrature.point(0 /*=quadraturePoint*/) );

    LocalFunctionType local_macro_disc = macro_discrete_function.localFunction(entity);
    JacobianRangeType grad_macro_discrete_function;
    local_macro_disc.jacobian(quadrature[0], grad_macro_discrete_function);

    PeriodicDiscreteFunctionType corrector_macro_discrete_function("corrector of macro discrete function",
                                                                   periodicDiscreteFunctionSpace_);
    discrete_function_reader.read(number_of_entity, corrector_macro_discrete_function);

    // the solution that we want to save to the data file
    PeriodicDiscreteFunctionType jac_corrector_Phi_i("jacobian corrector of Phi_i", periodicDiscreteFunctionSpace_);

    baseSet.jacobianAll(quadrature[0 /*=quadraturePoint*/], inv, gradientPhi);

    for (int i = 0; i < numBaseFunctions; ++i)
    {
      jac_corrector_Phi_i.clear();

      // take time
      DSC_PROFILER.startTiming("hmm.solver.saveTheJacCorSolutions_baseSet_discFunc.solve_jacobiancorrector_cellproblem");

      solve_jacobiancorrector_cellproblem
        (gradientPhi[i],
        grad_macro_discrete_function,
        corrector_macro_discrete_function,
        barycenter_of_entity,
        jac_corrector_Phi_i);

      // min/max time
      cell_time(DSC_PROFILER.stopTiming("hmm.solver.saveTheJacCorSolutions_baseSet_discFunc.solve_jacobiancorrector_cellproblem"));

      dfw.append(jac_corrector_Phi_i);

      // check if we use a correct numeration of the cell problems:
      const EntityPointerType entity_pointer(*it);
      if ( !(cp_num_manager.get_number_of_cell_problem(entity_pointer, i) == number_of_cell_problem) )
      {
        DUNE_THROW(Dune::InvalidStateException, "Numeration of cell problems incorrect.");
      }

      number_of_cell_problem++;
    }

    number_of_entity += 1;
  } // end: for-loop: IteratorType it

  const auto total_time = DSC_PROFILER.stopTiming("hmm.solver.saveTheJacCorSolutions_baseSet_discFunc");
  DSC_LOG_INFO << std::endl;
  DSC_LOG_INFO << "In method: saveTheJacCorSolutions_baseSet_discFunc." << std::endl << std::endl;
  DSC_LOG_INFO << "Cell problems solved for " << discreteFunctionSpace.grid().size(0) << " leaf entities."
                << std::endl;
  DSC_LOG_INFO << "Minimum time for solving a cell problem = " << cell_time.min() << "s." << std::endl;
  DSC_LOG_INFO << "Maximum time for solving a cell problem = " << cell_time.max() << "s." << std::endl;
  DSC_LOG_INFO << "Average time for solving a cell problem = "
                << cell_time.average()<< "s." << std::endl;
  DSC_LOG_INFO << "Total time for computing and saving the cell problems = "
                << total_time << "s," << std::endl << std::endl;
} // saveTheJacCorSolutions_baseSet_discFunc

} //namespace HMM {
} //namespace Multiscale {
} //namespace Dune {