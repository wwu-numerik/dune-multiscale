#ifndef DUNEMS_HMM_CELL_SOLVER_HH
#define DUNEMS_HMM_CELL_SOLVER_HH

#include <dune/stuff/configcontainer.hh>
#include <dune/multiscale/tools/disc_func_writer/discretefunctionwriter.hh>
#include <dune/multiscale/tools/solver/HMM/cell_problem_solving/discreteoperator.hh>
#include <dune/fem/operator/2order/lagrangematrixsetup.hh>

// ! ------------------------------------------------------------------------------------------------
// ! --------------------- the essential cell problem solver class ----------------------------------
namespace Dune {
template< class Traits >
class CellProblemSolver
{
public:
  // ! type of discrete functions
  typedef typename Traits::PeriodicDiscreteFunctionType PeriodicDiscreteFunctionType;

  // ! type of discrete function space
  typedef typename PeriodicDiscreteFunctionType::DiscreteFunctionSpaceType
  PeriodicDiscreteFunctionSpaceType;

  // ! type of grid partition
  typedef typename PeriodicDiscreteFunctionSpaceType::GridPartType PeriodicGridPartType;

  // ! type of grid
  typedef typename PeriodicDiscreteFunctionSpaceType::GridType PeriodicGridType;

  // ! type of range vectors
  typedef typename PeriodicDiscreteFunctionSpaceType::RangeType RangeType;

  // ! type of range vectors
  typedef typename PeriodicDiscreteFunctionSpaceType::DomainType DomainType;

  // ! polynomial order of base functions
  enum { polynomialOrder = PeriodicDiscreteFunctionSpaceType::polynomialOrder };

  // ! type of the (possibly non-linear) diffusion operator
  typedef typename Traits::DiffusionType DiffusionType;

  struct CellMatrixTraits
  {
    typedef PeriodicDiscreteFunctionSpaceType                          RowSpaceType;
    typedef PeriodicDiscreteFunctionSpaceType                          ColumnSpaceType;
    typedef LagrangeMatrixSetup< false >     StencilType;
    typedef ParallelScalarProduct< PeriodicDiscreteFunctionSpaceType > ParallelScalarProductType;

    template< class M >
    struct Adapter
    {
      typedef typename Traits::template LagrangeParallelMatrixAdapter< M > MatrixAdapterType;
    };
  };

  typedef SparseRowMatrixOperator< PeriodicDiscreteFunctionType, PeriodicDiscreteFunctionType,
                                   CellMatrixTraits > CellFEMMatrix;

  // OEMGMRESOp //OEMBICGSQOp // OEMBICGSTABOp
  typedef OEMBICGSTABOp< PeriodicDiscreteFunctionType, CellFEMMatrix > InverseCellFEMMatrix;

  // discrete elliptic operator describing the elliptic cell problems
  typedef DiscreteCellProblemOperator< PeriodicDiscreteFunctionType, DiffusionType > CellProblemOperatorType;

private:
  const PeriodicDiscreteFunctionSpaceType& periodicDiscreteFunctionSpace_; // Referenz &, wenn & verwendet, dann unten:
  DiffusionType& diffusion_;

  std::ofstream* data_file_;

public:
  // ! constructor - with diffusion operator A^{\epsilon}(x)
  CellProblemSolver(const PeriodicDiscreteFunctionSpaceType& periodicDiscreteFunctionSpace,
                    DiffusionType& diffusion_operator)
    : periodicDiscreteFunctionSpace_(periodicDiscreteFunctionSpace)
      , diffusion_(diffusion_operator)
      , data_file_(NULL)
  {}

  // ! constructor - with diffusion operator A^{\epsilon}(x)
  CellProblemSolver(const PeriodicDiscreteFunctionSpaceType& periodicDiscreteFunctionSpace,
                    DiffusionType& diffusion_operator,
                    std::ofstream& data_file)
    : periodicDiscreteFunctionSpace_(periodicDiscreteFunctionSpace)
      , diffusion_(diffusion_operator)
      , data_file_(&data_file)
  {}

  // ! ----------- method: solve cell problem ------------------------------------------

  template< class JacobianRangeImp >
  void solvecellproblem(JacobianRangeImp& gradient_PHI_H,
                        // the barycenter x_T of a macro grid element 'T'
                        const DomainType& globalQuadPoint,
                        PeriodicDiscreteFunctionType& cell_problem_solution) {
    // set solution equal to zero:
    cell_problem_solution.clear();

    // ! the matrix in our linear system of equations
    // in the non-linear case, it is the matrix for each iteration step
    CellFEMMatrix cell_system_matrix("Cell Problem System Matrix",
                                     periodicDiscreteFunctionSpace_,
                                     periodicDiscreteFunctionSpace_);

    // ! define the discrete (elliptic) cell problem operator
    // ( effect of the discretized differential operator on a certain discrete function )
    CellProblemOperatorType cell_problem_op(periodicDiscreteFunctionSpace_, diffusion_);

    // ! right hand side vector of the algebraic cell problem
    // (in the non-linear setting it changes for every iteration step)
    PeriodicDiscreteFunctionType cell_problem_rhs("rhs of cell problem", periodicDiscreteFunctionSpace_);
    cell_problem_rhs.clear();
    const bool CELLSOLVER_VERBOSE = Stuff::Config().get<bool>("problem.cellsolver_verbose", false);

    // NOTE:
    // is the right hand side of the cell problem equal to zero or almost identical to zero?
    // if yes, the solution of the cell problem is also identical to zero. The solver is getting a problem with this
    // situation, which is why we do not solve cell problems for zero-right-hand-side, since we already know the result.

    #ifdef LINEAR_PROBLEM

    // assemble the stiffness matrix
    cell_problem_op.assemble_matrix(globalQuadPoint, cell_system_matrix);

    // assemble right hand side of algebraic cell problem
    cell_problem_op.assembleCellRHS_linear(globalQuadPoint, gradient_PHI_H, cell_problem_rhs);

    const double norm_rhs = cell_problem_op.normRHS(cell_problem_rhs);

    if ( !( cell_problem_rhs.dofsValid() ) )
    {
      std::cout << "Cell Problem RHS invalid." << std::endl;
      abort();
    }

    if (norm_rhs < /*1e-06*/ 1e-10)
    {
      cell_problem_solution.clear();
      // std :: cout << "Cell problem with solution zero." << std :: endl;
    } else {
      InverseCellFEMMatrix cell_fem_biCGStab(cell_system_matrix, 1e-8, 1e-8, 20000, CELLSOLVER_VERBOSE);
      cell_fem_biCGStab(cell_problem_rhs, cell_problem_solution);
    }

    // end linear case.
    #else // ifdef LINEAR_PROBLEM
    // nonlinear case:

    // starting value for the Newton method
    PeriodicDiscreteFunctionType zero_func(" constant zero function ", periodicDiscreteFunctionSpace_);
    zero_func.clear();

    // residual vector (current residual)
    PeriodicDiscreteFunctionType cell_problem_residual("Cell Problem Residual", periodicDiscreteFunctionSpace_);
    cell_problem_residual.clear();

    L2Error< PeriodicDiscreteFunctionType > l2error;
    RangeType relative_newton_error = 10000.0;

    int iteration_step = 0;

    // the Newton step for for solving the current cell problem (solved with Newton Method):
    // L2-Norm of residual < tolerance ?
    #ifdef STOCHASTIC_PERTURBATION
    double tolerance = 1e-01 * VARIANCE;
    #else
    double tolerance = 1e-06;
    #endif // ifdef STOCHASTIC_PERTURBATION
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
        std::cout << "Cell Problem RHS invalid." << std::endl;
        abort();
      }

      if ( (norm_rhs < /*1e-06*/ 1e-10) /*|| ( (norm_rhs > 1e-06) && (norm_rhs < 1.5e-06) )*/ )
      {
        // residual solution almost identical to zero: break
        break;
      }

      double biCG_tolerance = 1e-8;
      bool cell_solution_convenient = false;
      while (cell_solution_convenient == false)
      {
        cell_problem_residual.clear();
        InverseCellFEMMatrix cell_fem_newton_biCGStab(cell_system_matrix,
                                                      1e-8, biCG_tolerance, 20000, CELLSOLVER_VERBOSE);

        cell_fem_newton_biCGStab(cell_problem_rhs, cell_problem_residual);

        if ( cell_problem_residual.dofsValid() )
        { cell_solution_convenient = true; }

        if (biCG_tolerance > 1e-4)
        {
          std::cout << "WARNING! Iteration step " << iteration_step
                    << ". Invalid dofs in 'cell_problem_residual', but '" << relative_newton_error
                    <<
          " = relative_newton_error > 1e-01' and 'biCG_tolerance > 1e-4'. L^2-Norm of right hand side of cell problem: "
                    << norm_rhs << ". Therefore possibly inaccurate solution." << std::endl;
          std::cout << "Information:" << std::endl;
          std::cout << "x_T = globalQuadPoint = " << globalQuadPoint << "." << std::endl;
          std::cout << "nabla u_H^{(n-1)} = gradient_PHI_H = " << gradient_PHI_H[0] << "." << std::endl;
          std::cout << "Print right hand side? y/n: ";
          char answer;
          std::cin >> answer;
          if ( !(answer == 'n') )
          { cell_problem_op.printCellRHS(cell_problem_rhs); }
          std::abort();
        }

        biCG_tolerance *= 10.0;
      }

      cell_problem_solution += cell_problem_residual;

      relative_newton_error = l2error.template norm2< (2* polynomialOrder) + 2 >(cell_problem_residual, zero_func);
      RangeType norm_cell_solution = l2error.template norm2< (2* polynomialOrder) + 2 >(cell_problem_solution,
                                                                                        zero_func);
      relative_newton_error = relative_newton_error / norm_cell_solution;

      // std :: cout << "L2-Norm of cell problem residual = " << residual_L2_norm << std :: endl;

      cell_problem_residual.clear();

      if (iteration_step > 10)
      {
        std::cout << "Warning! Algorithm already reached Newton-iteration step " << iteration_step
                  << " for computing the nonlinear cellproblem." << std::endl;
        std::cout << "relative_newton_error = " << relative_newton_error << std::endl << std::endl;
        #ifdef FORCE
        residual_L2_norm = 0.0;
        #endif
      }

      iteration_step += 1;
    }

    #endif // ifdef LINEAR_PROBLEM
    // end nonlinear case.

    if ( !( cell_problem_solution.dofsValid() ) )
    {
      std::cout << "Current solution of the cell problem invalid!" << std::endl;
      std::abort();
    }
  } // solvecellproblem

  // ! ----------- end method: solve cell problem ------------------------------------------

  // ! -------- method: solve the jacobian corrector cell problem ----------------------------
  // Problem to determine the Jacobian operator of the correction operator
  // ( the correction operator Q is defined via cell problems, here we determine the corresponding derivative DQ )
  // (see paper what it means and where it comes from. Note that it is only required for the nonlinear case)
  template< class JacobianRangeImp >
  void solve_jacobiancorrector_cellproblem(
    // gradient of macroscopic base function
    JacobianRangeImp& gradient_PHI_H,
    // gradient of the macroscopic function from the last iteration step
    JacobianRangeImp& grad_old_coarse_function,
    // gradient_y of the corrector of the macroscopic function from the last iteration step
    PeriodicDiscreteFunctionType& corrector_of_old_coarse_function,
    // the barycenter x_T of a macro grid element 'T'
    const DomainType& globalQuadPoint,
    PeriodicDiscreteFunctionType& jac_cor_cell_problem_solution) {
    // set solution equal to zero:
    jac_cor_cell_problem_solution.clear();

    // ! the matrix in our linear system of equations
    // system matrix for the jacobian corrector cell problem to solve:
    // entries:
    // - int_Y DA^{\eps}( x_T + \delta y, \nabla_x u_H^{(n-1)})(x_T) + \nabla_y Q_h( u_H^{(n-1)})(y) ) \nablay_y
    // \phi_h_i(y) \nablay_y \phi_h_j(y) dy
    // ( \phi_h_i and \phi_h_j denote microscopic base functions.)
    CellFEMMatrix jac_cor_cell_system_matrix("Jacobian Corrector Cell Problem System Matrix",
                                             periodicDiscreteFunctionSpace_,
                                             periodicDiscreteFunctionSpace_);

    // ! define the discrete (elliptic) cell problem operator
    // ( effect of the discretized differential operator on a certain discrete function )
    CellProblemOperatorType cell_problem_op(periodicDiscreteFunctionSpace_, diffusion_);
    // we are looking for the derivative of the operator in \nabla_x u_H^{(n-1)})(x_T) + \nabla_y Q_h( u_H^{(n-1)})(y)
    // and in direction of the macroscopic base function

    // ! right hand side vector of the algebraic jacobian corrector cell problem
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
      std::cout << "Jacobian Corrector Cell Problem RHS invalid." << std::endl;
      std::abort();
    }

    // is the right hand side of the jacobian corrector cell problem equal to zero or almost identical to zero?
    // if yes, the solution of the cell problem is also identical to zero. The solver is getting a problem with this
    // situation, which is why we do not solve cell problems for zero-right-hand-side, since we already know the result.
    if (norm_rhs < 1e-10)
    {
      jac_cor_cell_problem_solution.clear();
      // std :: cout << "Jacobian Corrector Cell Problem with solution zero." << std :: endl;
    } else {
      const bool CELLSOLVER_VERBOSE = Stuff::Config().get<bool>("problem.cellsolver_verbose", false);
      InverseCellFEMMatrix jac_cor_cell_fem_biCGStab(jac_cor_cell_system_matrix, 1e-8, 1e-8, 20000, CELLSOLVER_VERBOSE);
      jac_cor_cell_fem_biCGStab(jac_cor_cell_problem_rhs, jac_cor_cell_problem_solution);
    }
  } // solve_jacobiancorrector_cellproblem

  // ! ---------- method: solve the jacobian corrector cell problem -------------------------

  // two methods for solving and saving the solutions of the cell problems.
  // 1. save solutions for the whole set of macroscopic base function
  // (in general for the case of a linear diffusion operator)
  // 2. save the solutions for a fixed discrete function
  // (in general for the case of a nonlinear diffusion operator)

  // ! ---- method: solve and save the cell problems for the set of macroscopic base functions -----

  // here we need a 'cell problem numbering manager' to determine the number of the cell problem
  // (a combination of number of entity and number of local base function)
  // Structure:
  // Struktur der Indizierung fuer das Abspeichern der Loesungen der Zellprobleme:
  // wir loesen Zellprobleme fuer jede Entity des Makro-Grids und jede Basisfunktion, die einen nichtleeren support auf
  // dieser Entity besitzt, also schematisch:
  // Sei n=0,...,N einer Durchnumerierung der Entitys und i=0,...I_n eine zu einer festen Entity gehoerende Nummerierung
  // der Basisfunktionen mit nicht-leeren support.
  // Die Durchnummerierung der Loesungen der Zellprobleme k=0,...,K ist dann gegeben durch: k(n,i_n) = ( sum_(l=0)^(n-1)
  // ( I_l + 1) ) + i_n
  // NOTE: es verhaelt sich NICHT wie die vorhandene Methode mapToGlobal(entity,i) ! (die gibt die globale Nummer der
  // Basisfunktion zurueck, es gibt aber  deutlich mehr Zellprobleme zum Loesen!
  // (das wird aber alles im Hintergrund vom 'cell problem numbering manager')

  // compute and save solutions of the cell problems for the base function set of the 'discreteFunctionSpace'
  // requires cell problem numbering manager
  template< class DiscreteFunctionImp, class CellProblemNumberingManagerImp >
  void saveTheSolutions_baseSet(
    const typename DiscreteFunctionImp::DiscreteFunctionSpaceType& discreteFunctionSpace,
    const CellProblemNumberingManagerImp& cp_num_manager,   // just to check, if we use the correct numeration
    const std::string& filename) {
    typedef DiscreteFunctionImp DiscreteFunctionType;

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

    std::string cell_solution_location = "data/HMM/" + filename + "_cellSolutions_baseSet";
    DiscreteFunctionWriter dfw( (cell_solution_location).c_str() );

    const bool writer_is_open = dfw.is_open();

    long double starting_time = clock();

    // we want to determine minimum, average and maxiumum time for solving a cell problem in the current method
    double minimum_time_c_p = 1000000;
    double DUNE_UNUSED(average_time_c_p) = 0;
    double maximum_time_c_p = 0;

    int number_of_cell_problem = 0;

    if (writer_is_open)
    {
      IteratorType endit = discreteFunctionSpace.end();
      for (IteratorType it = discreteFunctionSpace.begin(); it != endit; ++it)
      {
        // gradients of the macroscopic base functions
        JacobianRangeType gradientPhi[maxnumOfBaseFct];

        // entity
        const EntityType& entity = *it;

        const BaseFunctionSetType baseSet
          = discreteFunctionSpace.baseFunctionSet(entity);

        const EntityGeometryType& geometry = entity.geometry();

        EntityQuadratureType quadrature(entity, 0);

        DomainType barycenter_of_entity = geometry.global( quadrature.point(0) );

        // number of base functions on entity
        const int numBaseFunctions = baseSet.numBaseFunctions();

        // calc Jacobian inverse before volume is evaluated
        const FieldMatrix< double, dimension, dimension >& inv
          = geometry.jacobianInverseTransposed( quadrature.point(0 /*=quadraturePoint*/) );

        PeriodicDiscreteFunctionType correctorPhi_i("corrector Phi_i", periodicDiscreteFunctionSpace_);

        for (int i = 0; i < numBaseFunctions; ++i)
        {
          baseSet.jacobian(i, quadrature[0 /*=quadraturePoint*/], gradientPhi[i]);
          // multiply with transpose of jacobian inverse
          gradientPhi[i][0] = FMatrixHelp::mult(inv, gradientPhi[i][0]);
        }

        for (int i = 0; i < numBaseFunctions; ++i)
        {
          correctorPhi_i.clear();

          // take time
          long double time_now = clock();

          solvecellproblem< JacobianRangeType >
            (gradientPhi[i], barycenter_of_entity, correctorPhi_i);

          // min/max time
          if ( (clock() - time_now) / CLOCKS_PER_SEC > maximum_time_c_p )
          { maximum_time_c_p = (clock() - time_now) / CLOCKS_PER_SEC; }
          if ( (clock() - time_now) / CLOCKS_PER_SEC < minimum_time_c_p )
          { minimum_time_c_p = (clock() - time_now) / CLOCKS_PER_SEC; }

          dfw.append(correctorPhi_i);

          // check if we use a correct numeration of the cell problems:
          if ( !(cp_num_manager.get_number_of_cell_problem(it, i) == number_of_cell_problem) )
          {
            std::cout << "Numeration of cell problems incorrect." << std::endl;
            std::abort();
          }

          number_of_cell_problem++;
        }
      } // end: for-loop: IteratorType it
    } // end: 'if ( writer_is_open )'

    if (data_file_)
    {
      if ( data_file_->is_open() )
      {
        (*data_file_) << std::endl;
        (*data_file_) << "In method: saveTheSolutions_baseSet." << std::endl << std::endl;
        (*data_file_) << "Cell problems solved for " << discreteFunctionSpace.grid().size(0) << " leaf entities."
                      << std::endl;
        (*data_file_) << "Minimum time for solving a cell problem = " << minimum_time_c_p << "s." << std::endl;
        (*data_file_) << "Maximum time for solving a cell problem = " << maximum_time_c_p << "s." << std::endl;
        (*data_file_) << "Average time for solving a cell problem = "
                      << ( (clock() - starting_time) / CLOCKS_PER_SEC ) / number_of_cell_problem << "s." << std::endl;
        (*data_file_) << "Total time for computing and saving the cell problems = "
                      << ( (clock() - starting_time) / CLOCKS_PER_SEC ) << "s," << std::endl << std::endl;
      }
    }
  } // saveTheSolutions_baseSet

  // ! ---- method: solve and save the cell problems for a fixed macroscopic discrete function -----

  // this method does not require a cell problem numbering manager
  // (it uses the standard counter for entities, provided by the corrsponding iterator)

  // compute and save solutions of the cell problems for a fixed macroscopic discrete function
  // (in gerneral it is the macro solution from the last iteration step)
  template< class DiscreteFunctionImp >
  void saveTheSolutions_discFunc(
    const DiscreteFunctionImp& macro_discrete_function,
    const std::string& filename) {
    typedef DiscreteFunctionImp DiscreteFunctionType;

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
    std::string cell_solution_location = "data/HMM/" + filename + "_cellSolutions_discFunc";
    DiscreteFunctionWriter dfw( (cell_solution_location).c_str() );

    long double starting_time = clock();

    // we want to determine minimum, average and maxiumum time for solving a cell problem in the current method
    double minimum_time_c_p = 1000000;
    double DUNE_UNUSED(average_time_c_p) = 0;
    double maximum_time_c_p = 0;

    const DiscreteFunctionSpaceType& discreteFunctionSpace = macro_discrete_function.space();

    if (dfw.is_open())
    {
      int number_of_entity = 0;
      IteratorType endit = discreteFunctionSpace.end();
      for (IteratorType it = discreteFunctionSpace.begin(); it != endit; ++it)
      {
        // entity
        const EntityType& entity = *it;

        const EntityGeometryType& geometry = entity.geometry();

        EntityQuadratureType quadrature(entity, 0);

        DomainType barycenter_of_entity = geometry.global( quadrature.point(0) );

        LocalFunctionType local_macro_disc = macro_discrete_function.localFunction(entity);
        JacobianRangeType grad_macro_discrete_function;
        local_macro_disc.jacobian(quadrature[0], grad_macro_discrete_function);

        PeriodicDiscreteFunctionType cell_solution_on_entity("corrector of macro discrete function",
                                                             periodicDiscreteFunctionSpace_);

        // take time
        long double time_now = clock();

        solvecellproblem< JacobianRangeType >
          (grad_macro_discrete_function, barycenter_of_entity, cell_solution_on_entity);

        // min/max time
        if ( (clock() - time_now) / CLOCKS_PER_SEC > maximum_time_c_p )
        { maximum_time_c_p = (clock() - time_now) / CLOCKS_PER_SEC; }
        if ( (clock() - time_now) / CLOCKS_PER_SEC < minimum_time_c_p )
        { minimum_time_c_p = (clock() - time_now) / CLOCKS_PER_SEC; }

        dfw.append(cell_solution_on_entity);

        #if 0
        // !LOESCHEN:
        if (number_of_entity == 341)
        {
          // in case you want to save the solutions of the two cell problems:
          #if 1
          typedef Tuple< PeriodicDiscreteFunctionImp* > IOTupleType;
          typedef DataOutput< GridType, IOTupleType >   DataOutputType;

          // general output parameters
          CellProblemDataOutputParameters outputparam;

          // sequence stamp
          std::stringstream outstring;

          // ------- cell problem -------------

          // create and initialize output class
          IOTupleType cellproblem_tuple(&cell_solution_on_entity);

          outputparam.set_prefix("cellSolution_saved_");
          outputparam.set_path("data/HMM/");
          DataOutputType cellSolution_dataoutput(periodicDiscreteFunctionSpace_.grid(), cellproblem_tuple, outputparam);

          // write data
          outstring << "cellSolution_saved_";
          cellSolution_dataoutput.writeData( 1.0 /*dummy*/, outstring.str() );
          // clear the std::stringstream:
          outstring.str( std::string() );

          #endif // if 1
        }
        // !-------------------------------------
        #endif // if 0

        number_of_entity += 1;
      } // end: for-loop: IteratorType it
    } // end: 'if ( writer_is_open )'

    if (data_file_)
    {
      if ( data_file_->is_open() )
      {
        (*data_file_) << std::endl;
        (*data_file_) << "In method: saveTheSolutions_discFunc." << std::endl << std::endl;
        (*data_file_) << "Cell problems solved for " << discreteFunctionSpace.grid().size(0) << " leaf entities."
                      << std::endl;
        (*data_file_) << "Minimum time for solving a cell problem = " << minimum_time_c_p << "s." << std::endl;
        (*data_file_) << "Maximum time for solving a cell problem = " << maximum_time_c_p << "s." << std::endl;
        (*data_file_) << "Average time for solving a cell problem = "
                      << ( (clock()
              - starting_time) / CLOCKS_PER_SEC ) / discreteFunctionSpace.grid().size(0) << "s." << std::endl;
        (*data_file_) << "Total time for computing and saving the cell problems = "
                      << ( (clock() - starting_time) / CLOCKS_PER_SEC ) << "s," << std::endl << std::endl;
      }
    }
  } // saveTheSolutions_discFunc

  #if 1

  // compute and save solutions of the jacobian corrector cell problems for the base function set of the
  // 'discreteFunctionSpace' and for a fixed macroscopic discrete function
  // requires cell problem numbering manager
  template< class DiscreteFunctionImp, class CellProblemNumberingManagerImp >
  void saveTheJacCorSolutions_baseSet_discFunc(
    const DiscreteFunctionImp& macro_discrete_function,
    const CellProblemNumberingManagerImp& cp_num_manager,   // just to check, if we use the correct numeration
    const std::string& filename) {
    typedef DiscreteFunctionImp DiscreteFunctionType;

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
    std::string cell_solution_location = "data/HMM/" + filename + "_JacCorCellSolutions_baseSet_discFunc";
    DiscreteFunctionWriter dfw( (cell_solution_location).c_str() );
    // where we saved the solutions for the discrete function
    // NOTE: they already need to be assembled, i.e. we already applied the method saveSolutions_discFunc!
    std::string cell_solution_discFunc_location = "data/HMM/" + filename + "_cellSolutions_discFunc";

    bool reader_is_open = false;

    // reader for the cell problem data file (discrete functions):
    DiscreteFunctionReader discrete_function_reader( (cell_solution_discFunc_location).c_str() );
    reader_is_open = discrete_function_reader.open();

    long double starting_time = clock();

    // we want to determine minimum, average and maxiumum time for solving a cell problem in the current method
    double minimum_time_c_p = 1000000;
//    double average_time_c_p = 0;
    double maximum_time_c_p = 0;

    int number_of_cell_problem = 0;

    const DiscreteFunctionSpaceType& discreteFunctionSpace = macro_discrete_function.space();

    if (dfw.is_open())
    {
      int number_of_entity = 0;
      IteratorType endit = discreteFunctionSpace.end();
      for (IteratorType it = discreteFunctionSpace.begin(); it != endit; ++it)
      {
        // gradients of the macroscopic base functions
        JacobianRangeType gradientPhi[maxnumOfBaseFct];

        // entity
        const EntityType& entity = *it;

        const BaseFunctionSetType baseSet
          = discreteFunctionSpace.baseFunctionSet(entity);

        const EntityGeometryType& geometry = entity.geometry();

        EntityQuadratureType quadrature(entity, 0);

        DomainType barycenter_of_entity = geometry.global( quadrature.point(0) );

        // number of base functions on entity
        const int numBaseFunctions = baseSet.size();

        // calc Jacobian inverse before volume is evaluated
        const FieldMatrix< double, dimension, dimension >& inv
          = geometry.jacobianInverseTransposed( quadrature.point(0 /*=quadraturePoint*/) );

        LocalFunctionType local_macro_disc = macro_discrete_function.localFunction(entity);
        JacobianRangeType grad_macro_discrete_function;
        local_macro_disc.jacobian(quadrature[0], grad_macro_discrete_function);

        PeriodicDiscreteFunctionType corrector_macro_discrete_function("corrector of macro discrete function",
                                                                       periodicDiscreteFunctionSpace_);
        if (reader_is_open)
        { discrete_function_reader.read(number_of_entity, corrector_macro_discrete_function); }

        // the solution that we want to save to the data file
        PeriodicDiscreteFunctionType jac_corrector_Phi_i("jacobian corrector of Phi_i", periodicDiscreteFunctionSpace_);

        for (int i = 0; i < numBaseFunctions; ++i)
        {
          baseSet.jacobian(i, quadrature[0 /*=quadraturePoint*/], gradientPhi[i]);
          // multiply with transpose of jacobian inverse
          gradientPhi[i][0] = FMatrixHelp::mult(inv, gradientPhi[i][0]);
        }

        for (int i = 0; i < numBaseFunctions; ++i)
        {
          jac_corrector_Phi_i.clear();

          // take time
          long double time_now = clock();

          solve_jacobiancorrector_cellproblem< JacobianRangeType >
            (gradientPhi[i],
            grad_macro_discrete_function,
            corrector_macro_discrete_function,
            barycenter_of_entity,
            jac_corrector_Phi_i);

          // min/max time
          if ( (clock() - time_now) / CLOCKS_PER_SEC > maximum_time_c_p )
          { maximum_time_c_p = (clock() - time_now) / CLOCKS_PER_SEC; }
          if ( (clock() - time_now) / CLOCKS_PER_SEC < minimum_time_c_p )
          { minimum_time_c_p = (clock() - time_now) / CLOCKS_PER_SEC; }

          dfw.append(jac_corrector_Phi_i);

          // check if we use a correct numeration of the cell problems:
          EntityPointerType entity_pointer(*it);
          if ( !(cp_num_manager.get_number_of_cell_problem(entity_pointer, i) == number_of_cell_problem) )
          {
            std::cout << "Numeration of cell problems incorrect." << std::endl;
            std::abort();
          }

          number_of_cell_problem++;
        }

        number_of_entity += 1;
      } // end: for-loop: IteratorType it
    } // end: 'if ( writer_is_open )'

    if (data_file_)
    {
      if ( data_file_->is_open() )
      {
        (*data_file_) << std::endl;
        (*data_file_) << "In method: saveTheJacCorSolutions_baseSet_discFunc." << std::endl << std::endl;
        (*data_file_) << "Cell problems solved for " << discreteFunctionSpace.grid().size(0) << " leaf entities."
                      << std::endl;
        (*data_file_) << "Minimum time for solving a cell problem = " << minimum_time_c_p << "s." << std::endl;
        (*data_file_) << "Maximum time for solving a cell problem = " << maximum_time_c_p << "s." << std::endl;
        (*data_file_) << "Average time for solving a cell problem = "
                      << ( (clock() - starting_time) / CLOCKS_PER_SEC ) / number_of_cell_problem << "s." << std::endl;
        (*data_file_) << "Total time for computing and saving the cell problems = "
                      << ( (clock() - starting_time) / CLOCKS_PER_SEC ) << "s," << std::endl << std::endl;
      }
    }
  } // saveTheJacCorSolutions_baseSet_discFunc

  #endif // if 1
}; // end class
} //namespace Dune {
#endif // DUNEMS_HMM_CELL_SOLVER_HH
