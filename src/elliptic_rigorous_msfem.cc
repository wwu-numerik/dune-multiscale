#include "common.hh"

#if HAVE_GRAPE
 #include <dune/grid/io/visual/grapedatadisplay.hh>
#endif

#include <dune/fem/gridpart/common/gridpart.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/space/lagrangespace.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/operator/lagrangeinterpolation.hh>

// to display data with ParaView:
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/fem/io/file/dataoutput.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/io/file/datawriter.hh>
#include <dune/fem/space/common/adaptmanager.hh>
#include <dune/fem/misc/l2error.hh>
#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/misc/h1norm.hh>
//#include <dune/fem/space/reducedbasisspace/reducedbasisspace.hh>

#include <dune/stuff/common/filesystem.hh>
//! local (dune-multiscale) includes
#include <dune/multiscale/problems/elliptic_problems/selector.hh>
#include <dune/multiscale/tools/solver/FEM/fem_solver.hh>

#include <dune/multiscale/tools/solver/MsFEM/msfem_localproblems/subgrid-list.hh>
#include <dune/multiscale/tools/solver/MsFEM/msfem_solver.hh>
#include <dune/multiscale/tools/meanvalue.hh>
#include <dune/multiscale/tools/improved_l2error.hh>
#include <dune/multiscale/tools/misc/h1error.hh>
#include <dune/multiscale/tools/disc_func_writer/discretefunctionwriter.hh>

//!-----------------------------------------------------------------------------

#include <dune/multiscale/tools/misc/outputparameter.hh>
#include <dune/multiscale/problems/elliptic_problems/selector.hh>

#include <dune/multiscale/msfem/rigorous_msfem_traits.hh>



#if 0

void solution_output(const Dune::MsfemTraits::DiscreteFunctionType& msfem_solution,
                     const Dune::MsfemTraits::DiscreteFunctionType& coarse_part_msfem_solution,
                     const Dune::MsfemTraits::DiscreteFunctionType& fine_part_msfem_solution,
                     Dune::myDataOutputParameters& outputparam,
                     const int loop_number,
                     int& total_refinement_level_,
                     int& coarse_grid_level_)
{
  using namespace Dune;
  //! ----------------- writing data output MsFEM Solution -----------------
  // --------- VTK data output for MsFEM solution --------------------------
  // create and initialize output class
  MsfemTraits::IOTupleType msfem_solution_series(&msfem_solution);
  const auto& gridPart = msfem_solution.space().gridPart();
  std::string outstring;
  outputparam.set_prefix("/msfem_solution");
  outstring = "msfem_solution";

  MsfemTraits::DataOutputType msfem_dataoutput(gridPart.grid(), msfem_solution_series, outputparam);
  msfem_dataoutput.writeData( 1.0 /*dummy*/, outstring );

  // create and initialize output class
  MsfemTraits::IOTupleType coarse_msfem_solution_series(&coarse_part_msfem_solution);

  outputparam.set_prefix("/coarse_part_msfem_solution");
  outstring = "coarse_part_msfem_solution";

  MsfemTraits::DataOutputType coarse_msfem_dataoutput(gridPart.grid(), coarse_msfem_solution_series, outputparam);
  coarse_msfem_dataoutput.writeData( 1.0 /*dummy*/, outstring );

  // create and initialize output class
  MsfemTraits::IOTupleType fine_msfem_solution_series(&fine_part_msfem_solution);

  outputparam.set_prefix("/fine_part_msfem_solution");
  // write data
  outstring = "fine_msfem_solution";

  MsfemTraits::DataOutputType fine_msfem_dataoutput(gridPart.grid(), fine_msfem_solution_series, outputparam);
  fine_msfem_dataoutput.writeData( 1.0 /*dummy*/, outstring);

  // ----------------------------------------------------------------------
  // ---------------------- write discrete msfem solution to file ---------
  const std::string location = (boost::format("msfem_solution_discFunc_refLevel_%d_%d")
                                %  total_refinement_level_ % coarse_grid_level_).str();
  DiscreteFunctionWriter(location).append(msfem_solution);
  //! --------------------------------------------------------------------
}

void data_output(const Dune::MsfemTraits::GridPartType& gridPart,
                 const Dune::MsfemTraits::DiscreteFunctionSpaceType& discreteFunctionSpace_coarse,
                 Dune::myDataOutputParameters& outputparam,
                 const int loop_number)
{
  using namespace Dune;
  // sequence stamp
  //! --------------------------------------------------------------------------------------

  //! -------------------------- writing data output Exact Solution ------------------------
  if (Problem::ModelProblemData::has_exact_solution)
  {
    const MsfemTraits::ExactSolutionType u;
    const MsfemTraits::DiscreteExactSolutionType discrete_exact_solution("discrete exact solution ", u, gridPart);
    // create and initialize output class
    MsfemTraits::ExSolIOTupleType exact_solution_series(&discrete_exact_solution);
    outputparam.set_prefix("/exact_solution");
    MsfemTraits::ExSolDataOutputType exactsol_dataoutput(gridPart.grid(), exact_solution_series, outputparam);
    // write data
    exactsol_dataoutput.writeData( 1.0 /*dummy*/, "exact-solution" );
    // -------------------------------------------------------
  }
  //! --------------------------------------------------------------------------------------

  //! --------------- writing data output for the coarse grid visualization ------------------
  MsfemTraits::DiscreteFunctionType coarse_grid_visualization("Visualization of the coarse grid",
                                                 discreteFunctionSpace_coarse);
  coarse_grid_visualization.clear();
  // -------------------------- data output -------------------------
  // create and initialize output class
  MsfemTraits::IOTupleType coarse_grid_series(&coarse_grid_visualization);

  const auto coarse_grid_fname = (boost::format("/coarse_grid_visualization_%d_") % loop_number).str();
  outputparam.set_prefix(coarse_grid_fname);
  MsfemTraits::DataOutputType coarse_grid_dataoutput(discreteFunctionSpace_coarse.gridPart().grid(), coarse_grid_series, outputparam);
  // write data
  coarse_grid_dataoutput.writeData( 1.0 /*dummy*/, coarse_grid_fname );
  // -------------------------------------------------------
  //! --------------------------------------------------------------------------------------
}

#endif


#if 1
//! algorithm
bool algorithm(const std::string& macroGridName,
               int& total_refinement_level_,
               int& coarse_grid_level_,
               int& number_of_layers_ ) {
  using namespace Dune;

  DSC_LOG_INFO << "loading dgf: " << macroGridName << std::endl;
  // we might use further grid parameters (depending on the grid type, e.g. Alberta), here we switch to default values
  // for the parameters:
  // create a grid pointer for the DGF file belongig to the macro grid:
  RigorousMsfemTraits::GridPointerType macro_grid_pointer(macroGridName);
  // refine the grid 'starting_refinement_level' times:
  macro_grid_pointer->globalRefine(coarse_grid_level_);
  //! ---- tools ----
  L2Error< RigorousMsfemTraits::DiscreteFunctionType > l2error;

  //! ---------------------------- grid parts ----------------------------------------------
  // grid part for the global function space, required for MsFEM-macro-problem
  RigorousMsfemTraits::GridPartType gridPart(*macro_grid_pointer);
  RigorousMsfemTraits::GridType& grid = gridPart.grid();
  //! --------------------------------------------------------------------------------------

  // coarse grid
  RigorousMsfemTraits::GridPointerType macro_grid_pointer_coarse(macroGridName);
  macro_grid_pointer_coarse->globalRefine(coarse_grid_level_);
  RigorousMsfemTraits::GridPartType gridPart_coarse(*macro_grid_pointer_coarse);
  RigorousMsfemTraits::GridType& grid_coarse = gridPart_coarse.grid();

  grid.globalRefine(total_refinement_level_ - coarse_grid_level_);

  //! ------------------------- discrete function spaces -----------------------------------
  // the global-problem function space:
  RigorousMsfemTraits::DiscreteFunctionSpaceType discreteFunctionSpace(gridPart);
  RigorousMsfemTraits::DiscreteFunctionSpaceType discreteFunctionSpace_coarse(gridPart_coarse);

  //! --------------------------- coefficient functions ------------------------------------

  // defines the matrix A^{\epsilon} in our global problem  - div ( A^{\epsilon}(\nabla u^{\epsilon} ) = f
  const RigorousMsfemTraits::DiffusionType diffusion_op;
  // define (first) source term:
  const RigorousMsfemTraits::FirstSourceType f; // standard source f


/// NEW RB TEST AREA
#if 1
  
  typedef BlockVector< FieldVector< double, 1> > VectorType; 
  typedef Matrix< FieldMatrix< double,1,1 > > MatrixType;
  
  RigorousMsfemTraits::RBSpace rb_space( discreteFunctionSpace_coarse );

  //! (stiffness) matrix
  int columns = 4;
  int rows = 4;
  int non_zero = columns * rows;
  MatrixType system_matrix( rows, columns );
  
  for (int row = 0; row != system_matrix.N(); ++row)
   for (int col = 0; col != system_matrix.M(); ++col)
     system_matrix[row][col] = 0.0;
   
  system_matrix[0][0] = 1.0;
  system_matrix[1][1] = 1.0;
  system_matrix[2][2] = 1.0;
  system_matrix[3][3] = 1.0;
  
  
  for (int row = 0; row != system_matrix.N(); ++row)
  { for (int col = 0; col != system_matrix.M(); ++col)
	 { std::cout << "[" << system_matrix[row][col] << "] "; }
    std::cout << std::endl;
  }
  
  std :: cout << "--------------------" << std ::endl;
  
  VectorType rhs( columns );
  for (int col = 0; col != rhs.N(); ++col)
    rhs[col] = 1.0;
  
  rhs[1] = 10.2723194;
  rhs[3] = -76.0000001;
  
  VectorType sol( columns );
  for (int col = 0; col != sol.N(); ++col)
    sol[col] = 0.0;
  
  typedef Dune::MatrixAdapter< MatrixType, VectorType, VectorType > MatrixOperatorType;
  MatrixOperatorType matrix_op( system_matrix );

  typedef Dune::SeqSSOR< MatrixType, VectorType, VectorType > PreconditionerType;
  PreconditionerType preconditioner( system_matrix, 100, 1.0 );
  
  typedef Dune::BiCGSTABSolver< VectorType > SolverType;
  Dune::InverseOperatorResult result_data;
  
  SolverType solver( matrix_op, preconditioner, 1e-10, 10000, true );
  //matrix_op.apply( rhs, sol);
  //std :: cout << "Done." << std ::endl;
  solver.apply(sol, rhs, result_data);
  
  
  for (int col = 0; col != sol.N(); ++col)
    { std::cout << "[" << sol[col] << "] "; }
  
  

//  SearchStrategyType search(source.gridPart().grid().leafView());
  const auto endit = discreteFunctionSpace.end();
  for(auto it = discreteFunctionSpace.begin(); it != endit ; ++it)
  {

    const auto& entity = *it;
    const auto& lagrangepoint_set = discreteFunctionSpace.lagrangePointSet(entity);
    
    const auto& geometry = entity.geometry();
    
    const auto number_of_points = lagrangepoint_set.nop();
    
    std::vector< RigorousMsfemTraits::DomainType > lagrange_points( number_of_points );
    for(int loc_point = 0; loc_point < number_of_points ; ++loc_point ) {
      int global_dof_number = mapToGlobal(entity, loc_point );
      lagrange_points[ loc_point ] = geometry.global(lagrangepoint_set.point( loc_point ) );
    }
#if 0

    auto target_local_function = target.localFunction(target_entity);



    typename TargetDiscreteFunctionSpaceType::RangeType source_value;


    const auto evaluation_entities = search(global_quads);
    assert(evaluation_entities.size() == global_quads.size());

    int k = 0;
    for(size_t qP = 0; qP < number_of_points ; ++qP)
    {
      if(std::isinf(target_local_function[ k ]))
      {
        const auto& global_point = global_quads[qP];
        // evaluate source function
        const auto source_entity = evaluation_entities[qP];
        const auto& source_geometry = source_entity->geometry();
        const auto& source_local_point = source_geometry.local(global_point);
        const auto& source_local_function = source.localFunction(*source_entity);
        source_local_function.evaluate(source_local_point, source_value);
        for(int i = 0; i < target_dimRange; ++i, ++k)
          target_local_function[k] = source_value[i];
      }
      else
        k += target_dimRange;
    }
#endif
  }

  
  abort();
  

#endif
  
#if 0
  //! ---------------------------- general output parameters ------------------------------
  // general output parameters
  Dune::myDataOutputParameters outputparam;
  data_output(gridPart, discreteFunctionSpace_coarse, outputparam, loop_number );
  
  //! ---------------------- solve MsFEM problem ---------------------------
  //! solution vector
  // solution of the standard finite element method
  MsfemTraits::DiscreteFunctionType msfem_solution("MsFEM Solution", discreteFunctionSpace);
  msfem_solution.clear();

  MsfemTraits::DiscreteFunctionType coarse_part_msfem_solution("Coarse Part MsFEM Solution", discreteFunctionSpace);
  coarse_part_msfem_solution.clear();

  MsfemTraits::DiscreteFunctionType fine_part_msfem_solution("Fine Part MsFEM Solution", discreteFunctionSpace);
  fine_part_msfem_solution.clear();

  const int number_of_level_host_entities = grid_coarse.size(0 /*codim*/);

  // number of layers per coarse grid entity T:  U(T) is created by enrichting T with n(T)-layers.
  MsfemTraits::MacroMicroGridSpecifierType specifier(discreteFunctionSpace_coarse, discreteFunctionSpace);
  for (int i = 0; i < number_of_level_host_entities; i += 1)
  {
    specifier.setLayer(i, number_of_layers_);
  }

  //default to false, error_estimation might change it
  bool repeat_algorithm = false;

  //! create subgrids:
  const bool silence = false;
  {//this scopes subgridlist
    MsfemTraits::SubGridListType subgrid_list(specifier, silence);

    // just for Dirichlet zero-boundary condition
    Elliptic_MsFEM_Solver< MsfemTraits::DiscreteFunctionType > msfem_solver(discreteFunctionSpace);
    msfem_solver.solve_dirichlet_zero(diffusion_op, f, specifier, subgrid_list,
                                      coarse_part_msfem_solution, fine_part_msfem_solution, msfem_solution);

    DSC_LOG_INFO << "Solution output for MsFEM Solution." << std::endl;
    solution_output(msfem_solution, coarse_part_msfem_solution,
                    fine_part_msfem_solution, outputparam, loop_number,
                    total_refinement_level_, coarse_grid_level_);

  }

  //! ---------------------- solve FEM problem with same (fine) resolution ---------------------------
  //! solution vector
  // solution of the standard finite element method
  MsfemTraits::DiscreteFunctionType fem_solution("FEM Solution", discreteFunctionSpace);
  fem_solution.clear();

  if ( DSC_CONFIG_GET("rigorous_msfem.fem_comparison",false) )
  {
    // just for Dirichlet zero-boundary condition
    const Elliptic_FEM_Solver< MsfemTraits::DiscreteFunctionType > fem_solver(discreteFunctionSpace);
    fem_solver.solve_dirichlet_zero(diffusion_op, f, fem_solution);

  //! ----------------------------------------------------------------------
  DSC_LOG_INFO << "Data output for FEM Solution." << std::endl;
  //! -------------------------- writing data output FEM Solution ----------

  // ------------- VTK data output for FEM solution --------------
  // create and initialize output class
  MsfemTraits::IOTupleType fem_solution_series(&fem_solution);
  outputparam.set_prefix("/fem_solution");
  MsfemTraits::DataOutputType fem_dataoutput(gridPart.grid(), fem_solution_series, outputparam);

    // write data
    fem_dataoutput.writeData( 1.0 /*dummy*/, "fem_solution" );
    // -------------------------------------------------------------
  }

  DSC_LOG_INFO << std::endl << "The L2 errors:" << std::endl << std::endl;
  //! ----------------- compute L2- and H1- errors -------------------
  if (Problem::ModelProblemData::has_exact_solution)
  {

    H1Error< MsfemTraits::DiscreteFunctionType > h1error;

    const MsfemTraits::ExactSolutionType u;

    MsfemTraits::RangeType msfem_error = l2error.norm< MsfemTraits::ExactSolutionType >(u,
                                                              msfem_solution,
                                                              2 * MsfemTraits::DiscreteFunctionSpaceType::polynomialOrder + 2);
    DSC_LOG_INFO << "|| u_msfem - u_exact ||_L2 =  " << msfem_error << std::endl << std::endl;

    MsfemTraits::RangeType h1_msfem_error(0.0);
    h1_msfem_error = h1error.semi_norm< MsfemTraits::ExactSolutionType >(u, msfem_solution);
    h1_msfem_error += msfem_error;
    DSC_LOG_INFO << "|| u_msfem - u_exact ||_H1 =  " << h1_msfem_error << std::endl << std::endl;

    if ( DSC_CONFIG_GET("rigorous_msfem.fem_comparison",false) )
    {
      MsfemTraits::RangeType fem_error = l2error.norm< MsfemTraits::ExactSolutionType >(u,
                                                            fem_solution,
                                                            2 * MsfemTraits::DiscreteFunctionSpaceType::polynomialOrder + 2);

      DSC_LOG_INFO << "|| u_fem - u_exact ||_L2 =  " << fem_error << std::endl << std::endl;

      MsfemTraits::RangeType h1_fem_error(0.0);

      h1_fem_error = h1error.semi_norm< MsfemTraits::ExactSolutionType >(u, fem_solution);
      h1_fem_error += fem_error;
      DSC_LOG_INFO << "|| u_fem - u_exact ||_H1 =  " << h1_fem_error << std::endl << std::endl;
    }
  } else if ( DSC_CONFIG_GET("rigorous_msfem.fem_comparison",false) )
  {
    DSC_LOG_ERROR << "Exact solution not available. Errors between MsFEM and FEM approximations for the same fine grid resolution."
                  << std::endl << std::endl;
    MsfemTraits::RangeType approx_msfem_error = l2error.norm2< 2* MsfemTraits::DiscreteFunctionSpaceType::polynomialOrder + 2 >(fem_solution,
                                                                                                      msfem_solution);
    DSC_LOG_INFO << "|| u_msfem - u_fem ||_L2 =  " << approx_msfem_error << std::endl << std::endl;
    H1Norm< MsfemTraits::GridPartType > h1norm(gridPart);
    MsfemTraits::RangeType h1_approx_msfem_error = h1norm.distance(fem_solution, msfem_solution);

    DSC_LOG_INFO << "|| u_msfem - u_fem ||_H1 =  " << h1_approx_msfem_error << std::endl << std::endl;
  }
  
  //! -------------------------------------------------------
  if (DSC_CONFIG_GET("adaptive", false)) {
    DSC_LOG_INFO << std::endl << std::endl;
    DSC_LOG_INFO << "---------------------------------------------" << std::endl;
  }
  return repeat_algorithm;
  
#endif
  return true;
} // function algorithm

#endif

int main(int argc, char** argv) {
  try {
    init(argc, argv);
    namespace DSC = Dune::Stuff::Common;
    //!TODO include base in config
    DSC_PROFILER.startTiming("msfem.all");

    const std::string datadir = DSC_CONFIG_GET("global.datadir", "data/");

    // generate directories for data output
    DSC::testCreateDirectory(datadir);

    DSC_LOG_INFO << boost::format("Data will be saved under: %s\nLogs will be saved under: %s/%s/ms.log.log\n")
                            % datadir % datadir % DSC_CONFIG_GET("logging.dir", "log");

    // syntax: info_from_par_file / default  / validation of the value

    // coarse_grid_level denotes the (starting) grid refinement level for the global coarse scale problem, i.e. it describes 'H'
    int coarse_grid_level_ = DSC_CONFIG_GETV( "rigorous_msfem.coarse_grid_level", 4, DSC::ValidateLess< int >( -1 ) );

    // syntax: info_from_par_file / default
    int number_of_layers_ = DSC_CONFIG_GET("rigorous_msfem.oversampling_layers", 4);
    
    // data for the model problem; the information manager
    // (see 'problem_specification.hh' for details)
    const Problem::ModelProblemData info;


    // total_refinement_level denotes the (starting) grid refinement level for the global fine scale problem, i.e. it describes 'h'
    int total_refinement_level_
      = DSC_CONFIG_GETV( "rigorous_msfem.fine_grid_level", 4, DSC::ValidateLess< int >(coarse_grid_level_-1) );

    // name of the grid file that describes the macro-grid:
    const std::string macroGridName = info.getMacroGridFile();


    DSC_LOG_INFO << "Error File for Elliptic Model Problem " << Dune::Stuff::Common::getTypename(info)
              << " with epsilon = " << DSC_CONFIG_GET("problem.epsilon", 1.0f) << "." << std::endl << std::endl;
    if ( DSC_CONFIG_GET("rigorous_msfem.petrov_galerkin", true) )
        DSC_LOG_INFO << "Use New Rigorous MsFEM in Petrov-Galerkin formulation with an uniform computation, i.e.:" << std::endl;
    else
        DSC_LOG_INFO << "Use New Rigorous MsFEM in classical (symmetric) formulation with an uniform computation, i.e.:" << std::endl;      
    DSC_LOG_INFO << "Uniformly refined coarse and fine mesh and" << std::endl;
    DSC_LOG_INFO << "the same number of layers for each (oversampled) local grid computation." << std::endl << std::endl;
    DSC_LOG_INFO << "Computations were made for:" << std::endl << std::endl;
    DSC_LOG_INFO << "Refinement Level for (uniform) Fine Grid = " << total_refinement_level_ << std::endl;
    DSC_LOG_INFO << "Refinement Level for (uniform) Coarse Grid = " << coarse_grid_level_ << std::endl;
    //DSC_LOG_INFO << "Oversampling Strategy = " << DSC_CONFIG_GET( "rigorous_msfem.oversampling_strategy", 1 ) << std::endl;
    DSC_LOG_INFO << "Number of layers for oversampling = " << number_of_layers_ << std::endl;
    if ( DSC_CONFIG_GET("rigorous_msfem.fem_comparison",false) )
       { DSC_LOG_INFO << std::endl << "Comparison with standard FEM computation on the MsFEM Fine Grid, i.e. on Refinement Level " << total_refinement_level_ << std::endl; }
    DSC_LOG_INFO << std::endl << std::endl;

    algorithm(macroGridName, total_refinement_level_, coarse_grid_level_, number_of_layers_ );

    // the reference problem generaly has a 'refinement_difference_for_referenceproblem' higher resolution than the
    //normal
    // macro problem

    const auto cpu_time = DSC_PROFILER.stopTiming("msfem.all");
    DSC_LOG_INFO << "Total runtime of the program: " << cpu_time << "ms" << std::endl;
    DSC_PROFILER.outputTimings("profiler");
    return 0;
  } catch (Dune::Exception& e) {
    std::cerr << e.what() << std::endl;
  }
  return 1;
} // main
