// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#ifdef HAVE_CMAKE_CONFIG
 #include "cmake_config.h"
#elif defined (HAVE_CONFIG_H)
 #include <config.h>
#endif // ifdef HAVE_CMAKE_CONFIG

#include "algorithm.hh"

#include <dune/multiscale/msfem/msfem_traits.hh>

#include <dune/fem/gridpart/common/gridpart.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/space/lagrange.hh>
//#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/operator/lagrangeinterpolation.hh>

// to display data with ParaView:
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/fem/io/file/dataoutput.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/io/file/datawriter.hh>
#include <dune/fem/space/common/adaptmanager.hh>
#include <dune/multiscale/common/error_calc.hh>

#include <dune/multiscale/common/righthandside_assembler.hh>
#include <dune/multiscale/problems/selector.hh>
#include <dune/multiscale/msfem/fem_solver.hh>
#include <dune/multiscale/msfem/localproblems/subgrid-list.hh>
#include <dune/multiscale/lod/lod_solver.hh>
#include <dune/multiscale/tools/meanvalue.hh>
#include <dune/multiscale/tools/improved_l2error.hh>
#include <dune/multiscale/hmm/cell_problem_numbering.hh>
#include <dune/multiscale/tools/misc/outputparameter.hh>
#include <dune/multiscale/msfem/msfem_grid_specifier.hh>
#include <dune/multiscale/common/output_traits.hh>

namespace Dune {
namespace Multiscale {
namespace MsFEM {

//! set the dirichlet points to the Dirichlet BC
template< class DirichletBC, class DiscreteFunctionType >
void setDirichletValues(MsFEMTraits::MacroMicroGridSpecifierType& specifier, DirichletBC &dirichlet_func, DiscreteFunctionType& func) {
  using namespace Dune::Stuff;
  const auto& discreteFunctionSpace = func.space();
  static const unsigned int faceCodim = 1;
  for (const auto& entity : discreteFunctionSpace)
  {

     // we might look for a better Dirichlet Extension in the case of a Neumann boundary part ..

#if 1
     auto funcLocal = func.localFunction(entity);

     // ------------- first, use a coarse grid interpolation ----------------
     bool intersects_dirichlet_boundary = false;

     const int level_difference = specifier.getLevelDifference();
     CommonTraits::DiscreteFunctionSpaceType& coarseSpace = specifier.coarseSpace();

     typedef typename CommonTraits::GridType::Traits::LeafIndexSet HostGridLeafIndexSet;
     const HostGridLeafIndexSet& coarseGridLeafIndexSet = coarseSpace.gridPart().grid().leafIndexSet();

     const auto& entity_ptr = coarseSpace.grid().entityPointer( entity.seed() );
     const CommonTraits::EntityPointerType& coarse_father = Stuff::Grid::make_father( coarseGridLeafIndexSet,
                                                                                      entity_ptr,
                                                                                      level_difference );
     CommonTraits::EntityType& coarse_grid_entity = *coarse_father;

     const auto& lagrangePointSet = specifier.coarseSpace().lagrangePointSet( coarse_grid_entity );
     const auto number_of_points = lagrangePointSet.nop();
     
     std::vector< std::size_t > indices;
     specifier.coarseSpace().mapper().map( coarse_grid_entity, indices );

#if 0
     int neumann_element = 0;
     int num_boundary_nodes = 0;
     for(size_t loc_point = 0; loc_point < number_of_points ; ++loc_point ) {
       
        if ( specifier.is_coarse_boundary_node( indices[ loc_point ] ) )
        {
          num_boundary_nodes += 1;
          if ( !specifier.is_coarse_dirichlet_node( indices[ loc_point ] ) )
           neumann_element += 1;
        }
      }
#if 0
     for (const auto& coarse_intersection
         : Dune::Stuff::Common::intersectionRange(coarseSpace.gridPart(), coarse_grid_entity ))
     {
        std::cout << "coarse_intersection.boundary() = " << coarse_intersection.boundary() << std::endl;
        std::cout << "coarse_intersection.boundaryId() = " << coarse_intersection.boundaryId() << std::endl;
        if ( coarse_intersection.boundary() && (coarse_intersection.boundaryId() == 2) )
        { neumann_element = true; }
     }
#endif
     if ( (neumann_element >= 1) && (num_boundary_nodes >= 2) )
     { continue; }
#endif

     assert( number_of_points == 3 );
     std::vector< CommonTraits::RangeType > dirichlet_values_coarse_corners( number_of_points );
     std::vector< CommonTraits::DomainType > corners( number_of_points );
     std::vector< CommonTraits::DomainType > dirichlet_corners;
     
     for(size_t loc_point = 0; loc_point < number_of_points ; ++loc_point ) {
        corners[ loc_point ] = coarse_grid_entity.geometry().global(lagrangePointSet.point( loc_point ) );
        
        if ( specifier.is_coarse_dirichlet_node( indices[ loc_point ] ) )
        {  dirichlet_func.evaluate( corners[ loc_point ], dirichlet_values_coarse_corners[ loc_point ] ); 
           intersects_dirichlet_boundary = true;
           dirichlet_corners.push_back( corners[ loc_point ] );
        }
        else
        { dirichlet_values_coarse_corners[ loc_point ] = 0.0; }
      }
    
     // LinearLagrangeInterpolation2D should be eventually replaced by
     // LinearLagrangeFunction2D< DiscreteFunctionSpace > dirichlet_interpolation
     LinearLagrangeInterpolation2D< CommonTraits::DiscreteFunctionSpaceType > dirichlet_interpolation
          ( corners[0], dirichlet_values_coarse_corners[0],
            corners[1], dirichlet_values_coarse_corners[1],
            corners[2], dirichlet_values_coarse_corners[2] );
     // --------------------------------------------------------------------------------
  
     if ( intersects_dirichlet_boundary == true ) 
     {
       for (const auto& intersection
         : Dune::Stuff::Common::intersectionRange(discreteFunctionSpace.gridPart(), entity))
       {

         const auto face = intersection.indexInInside();

         for(auto loc_point
           : Dune::Stuff::Common::lagrangePointSetRange<faceCodim>(func.space(), entity, face))
         {
           const auto& global_point = entity.geometry().global( discreteFunctionSpace.lagrangePointSet(entity).point( loc_point ) );
           CommonTraits::RangeType dirichlet_value(0.0);
           
           if ( intersection.boundary() && (intersection.boundaryId() == 2) )
           {
             bool point_is_dirichlet_point_on_neumann_boundary_piece = false;
             for ( int c = 0; c < dirichlet_corners.size(); c++ )
             {
                if ( global_point == dirichlet_corners[c] )
                { point_is_dirichlet_point_on_neumann_boundary_piece = true; }
             }
             if ( point_is_dirichlet_point_on_neumann_boundary_piece )
             { dirichlet_func.evaluate( global_point, dirichlet_value); }
             else
             { dirichlet_value = 0.0; }
             funcLocal[loc_point] = dirichlet_value;
           }
           else
           {
             // use the coarse grid interpolation, but on the Dirichlet boundary, use a fine grid interpolation
             if ( !intersection.boundary() )
             { 
               dirichlet_interpolation.evaluate( global_point, dirichlet_value);
             }
             else
             { dirichlet_func.evaluate( global_point, dirichlet_value); }
             funcLocal[loc_point] = dirichlet_value;
           }
         }
       }
     }
#endif
#if 0
    for (const auto& intersection
         : Dune::Stuff::Common::intersectionRange(discreteFunctionSpace.gridPart(), entity))
    {
      if ( !intersection.boundary() )
        continue;
      if ( intersection.boundary() && (intersection.boundaryId() != 1) )
        continue;

      auto funcLocal = func.localFunction(entity);
      const auto face = intersection.indexInInside();
      for(auto loc_point
          : Dune::Stuff::Common::lagrangePointSetRange<faceCodim>(func.space(), entity, face))
      {
        const auto& global_point = entity.geometry().global( discreteFunctionSpace.lagrangePointSet(entity).point( loc_point ) );
        CommonTraits::RangeType dirichlet_value(0.0);
        dirichlet_func.evaluate( global_point, dirichlet_value);
        funcLocal[loc_point] = dirichlet_value;
      }
    }
#endif
  }
} // setDirichletValues
  
//! \TODO docme
void solution_output(const CommonTraits::DiscreteFunctionType& lod_solution,
                     const CommonTraits::DiscreteFunctionType& coarse_part_lod_solution,
                     const CommonTraits::DiscreteFunctionType& fine_part_lod_solution,
                     const CommonTraits::DiscreteFunctionType& dirichlet_extension,
                     Dune::Multiscale::OutputParameters& outputparam,
                     int& total_refinement_level_,
                     int& coarse_grid_level_)
{
  using namespace Dune;
  //! ----------------- writing data output LOD Solution -----------------
  // --------- VTK data output for LOD solution --------------------------
  // create and initialize output class
  OutputTraits::IOTupleType lod_solution_series(&lod_solution);
  const auto& gridPart = lod_solution.space().gridPart();
  std::string outstring;
  outputparam.set_prefix("lod_solution");
  outstring = "lod_solution";

  OutputTraits::DataOutputType lod_dataoutput(gridPart.grid(), lod_solution_series, outputparam);
  lod_dataoutput.writeData( 1.0 /*dummy*/, outstring );

  // ---------------------------------------------------------------------------
  // create and initialize output class
  OutputTraits::IOTupleType coarse_lod_solution_series(&coarse_part_lod_solution);

  outputparam.set_prefix("coarse_part_lod_solution");
  outstring = "coarse_part_lod_solution";

  OutputTraits::DataOutputType coarse_lod_dataoutput(gridPart.grid(), coarse_lod_solution_series, outputparam);
  coarse_lod_dataoutput.writeData( 1.0 /*dummy*/, outstring );

  // ---------------------------------------------------------------------------
  // create and initialize output class
  OutputTraits::IOTupleType fine_lod_solution_series(&fine_part_lod_solution);

  outputparam.set_prefix("fine_part_lod_solution");
  // write data
  outstring = "fine_lod_solution";

  OutputTraits::DataOutputType fine_lod_dataoutput(gridPart.grid(), fine_lod_solution_series, outputparam);
  fine_lod_dataoutput.writeData( 1.0 /*dummy*/, outstring);
  
  // ---------------------------------------------------------------------------
  // create and initialize output class
  OutputTraits::IOTupleType dirichlet_extension_series(&dirichlet_extension);

  outputparam.set_prefix("dirichlet_extension");
  // write data
  outstring = "dirichlet_extension";

  OutputTraits::DataOutputType dirichlet_extension_dataoutput(gridPart.grid(), dirichlet_extension_series, outputparam);
  dirichlet_extension_dataoutput.writeData( 1.0 /*dummy*/, outstring);

  // ----------------------------------------------------------------------
  // ---------------------- write discrete lod solution to file ---------
  const std::string location = (boost::format("msfem_solution_discFunc_refLevel_%d_%d")
                                %  total_refinement_level_ % coarse_grid_level_).str();
  DiscreteFunctionWriter(location).append(lod_solution);
  //! --------------------------------------------------------------------
}

//! \TODO docme
void data_output(const CommonTraits::GridPartType& gridPart,
                 const CommonTraits::DiscreteFunctionSpaceType& discreteFunctionSpace_coarse,
                 Dune::Multiscale::OutputParameters& outputparam )
{
  using namespace Dune;
  // sequence stamp
  //! --------------------------------------------------------------------------------------

  //! -------------------------- writing data output Exact Solution ------------------------
  if (Problem::getModelData()->hasExactSolution())
  { 
    auto u_ptr = Dune::Multiscale::Problem::getExactSolution();
    const auto& u = *u_ptr;
    const OutputTraits::DiscreteExactSolutionType discrete_exact_solution("discrete exact solution ", u, gridPart);
    // create and initialize output class
    OutputTraits::ExSolIOTupleType exact_solution_series(&discrete_exact_solution);
    outputparam.set_prefix("exact_solution");
    OutputTraits::ExSolDataOutputType exactsol_dataoutput(gridPart.grid(), exact_solution_series, outputparam);
    // write data
    exactsol_dataoutput.writeData( 1.0 /*dummy*/, "exact-solution" );
    // -------------------------------------------------------
  }
  //! --------------------------------------------------------------------------------------

  //! --------------- writing data output for the coarse grid visualization ------------------
  CommonTraits::DiscreteFunctionType coarse_grid_visualization("Visualization of the coarse grid",
                                                 discreteFunctionSpace_coarse);
  coarse_grid_visualization.clear();
  // -------------------------- data output -------------------------
  // create and initialize output class
  OutputTraits::IOTupleType coarse_grid_series(&coarse_grid_visualization);

  const auto coarse_grid_fname = (boost::format("coarse_grid_visualization_")).str();
  outputparam.set_prefix(coarse_grid_fname);
  OutputTraits::DataOutputType coarse_grid_dataoutput(discreteFunctionSpace_coarse.gridPart().grid(), coarse_grid_series, outputparam);
  // write data
  coarse_grid_dataoutput.writeData( 1.0 /*dummy*/, coarse_grid_fname );
  // -------------------------------------------------------
  //! --------------------------------------------------------------------------------------
}


//! \TODO docme
void algorithm(const std::string& macroGridName,
               int& total_refinement_level_,
               int& coarse_grid_level_,
               int& number_of_layers_ ) {
  using namespace Dune;

  DSC_LOG_INFO << "loading dgf: " << macroGridName << std::endl;
  // we might use further grid parameters (depending on the grid type, e.g. Alberta), here we switch to default values
  // for the parameters:
  // create a grid pointer for the DGF file belongig to the macro grid:
  CommonTraits::GridPointerType macro_grid_pointer(macroGridName);
  // refine the grid 'starting_refinement_level' times:
  macro_grid_pointer->globalRefine(coarse_grid_level_);
  //! ---- tools ----
  L2Error< CommonTraits::DiscreteFunctionType > l2error;

  //! ---------------------------- grid parts ----------------------------------------------
  // grid part for the global function space, required for LOD-macro-problem
  CommonTraits::GridPartType gridPart(*macro_grid_pointer);
  CommonTraits::GridType& grid = gridPart.grid();
  //! --------------------------------------------------------------------------------------

  // coarse grid
  CommonTraits::GridPointerType macro_grid_pointer_coarse(macroGridName);
  macro_grid_pointer_coarse->globalRefine(coarse_grid_level_);
  CommonTraits::GridPartType gridPart_coarse(*macro_grid_pointer_coarse);
  CommonTraits::GridType& grid_coarse = gridPart_coarse.grid();

  grid.globalRefine(total_refinement_level_ - coarse_grid_level_);

  //! ------------------------- discrete function spaces -----------------------------------
  // the global-problem function space:
  CommonTraits::DiscreteFunctionSpaceType discreteFunctionSpace(gridPart);
  CommonTraits::DiscreteFunctionSpaceType discreteFunctionSpace_coarse(gridPart_coarse);

  //! --------------------------- coefficient functions ------------------------------------

  // defines the matrix A^{\epsilon} in our global problem
  //    - div ( A^{\epsilon}(\nabla u^{\epsilon} ) + F(x,u^{\epsilon},\nabla u^{\epsilon}) = f
  auto diffusion_op_ptr = Dune::Multiscale::Problem::getDiffusion();
  const auto& diffusion_op = *diffusion_op_ptr;
  // define (first) source term:
  auto f_ptr = Dune::Multiscale::Problem::getFirstSource();
  const auto& f = *f_ptr;
  const auto F_ptr = Dune::Multiscale::Problem::getLowerOrderTerm();
  const auto& F = *F_ptr; // lower term F(x,u ,\nabla u )
  // Dirichlet boundary condition
  const auto dirichlet_bc_ptr = Problem::getDirichletBC();
  const auto& dirichlet_bc = *dirichlet_bc_ptr; 
  // Neumann boundary condition
  const auto neumann_bc_ptr = Problem::getNeumannBC();
  const auto& neumann_bc = *neumann_bc_ptr;

  //! ---------------------------- general output parameters ------------------------------
  // general output parameters
  Dune::Multiscale::OutputParameters outputparam;
  data_output(gridPart, discreteFunctionSpace_coarse, outputparam );

  //! ---------------------- solve LOD problem ---------------------------
  //! solution vector
  // solution of the standard finite element method
  CommonTraits::DiscreteFunctionType msfem_solution("LOD Solution", discreteFunctionSpace);
  msfem_solution.clear();

  CommonTraits::DiscreteFunctionType coarse_part_msfem_solution("Coarse Part LOD Solution", discreteFunctionSpace);
  coarse_part_msfem_solution.clear();

  CommonTraits::DiscreteFunctionType fine_part_msfem_solution("Fine Part LOD Solution", discreteFunctionSpace);
  fine_part_msfem_solution.clear();

  const int number_of_level_host_entities = grid_coarse.size(0 /*codim*/);

  // number of layers per coarse grid entity T:  U(T) is created by enrichting T with n(T)-layers.
  MsFEMTraits::MacroMicroGridSpecifierType specifier(discreteFunctionSpace_coarse, discreteFunctionSpace);
  for (int i = 0; i < number_of_level_host_entities; i += 1)
  {
    specifier.setNoOfLayers(i, number_of_layers_);
  }
  //! \todo Important why? (Sven)
  specifier.setOversamplingStrategy( 3 ); //! Important!
  specifier.identify_coarse_dirichlet_nodes();
  //!specifier.identify_coarse_boundary_nodes();//!DELETE ME SOON!!!!!!!!!!!!

  // discrete function that takes the values of the Dirichlet BC on the Dirichlet Boundary nodes
  // and that is zero elsewhere
  CommonTraits::DiscreteFunctionType dirichlet_extension("Dirichlet extension", discreteFunctionSpace);
  dirichlet_extension.clear();
  setDirichletValues( specifier, dirichlet_bc, dirichlet_extension );

  //! create subgrids:
  {//this scopes subgridlist
    MsFEMTraits::SubGridListType subgrid_list(specifier, DSC_CONFIG_GET("logging.subgrid_silent", false));

    // just for Dirichlet zero-boundary condition
    Elliptic_Rigorous_MsFEM_Solver lod_solver(discreteFunctionSpace);
    lod_solver.solve(diffusion_op, f, dirichlet_extension, neumann_bc,
                     specifier, subgrid_list,
                     coarse_part_msfem_solution, fine_part_msfem_solution, msfem_solution);

    coarse_part_msfem_solution += dirichlet_extension;
    msfem_solution += dirichlet_extension;

    DSC_LOG_INFO << "Solution output for MsFEM Solution." << std::endl;
    solution_output(msfem_solution, coarse_part_msfem_solution,
                    fine_part_msfem_solution, dirichlet_extension, outputparam,
                    total_refinement_level_, coarse_grid_level_);
  }

  //! ---------------------- solve FEM problem with same (fine) resolution ---------------------------
  //! solution vector
  // solution of the standard finite element method
  CommonTraits::DiscreteFunctionType fem_solution("FEM Solution", discreteFunctionSpace);
  fem_solution.clear();

  if ( DSC_CONFIG_GET("rigorous_msfem.fem_comparison",false) )
  {
   
    // just for Dirichlet zero-boundary condition
    const Elliptic_FEM_Solver fem_solver(discreteFunctionSpace);
    fem_solver.solve(diffusion_op, F_ptr, f, dirichlet_extension, neumann_bc, fem_solution);
    fem_solution += dirichlet_extension;
    fem_solution.communicate();
    //! ----------------------------------------------------------------------
    DSC_LOG_INFO << "Data output for FEM Solution." << std::endl;
    //! -------------------------- writing data output FEM Solution ----------

    // ------------- VTK data output for FEM solution --------------
    // create and initialize output class
    OutputTraits::IOTupleType fem_solution_series(&fem_solution);
    outputparam.set_prefix("fem_solution");
    OutputTraits::DataOutputType fem_dataoutput(gridPart.grid(), fem_solution_series, outputparam);

    // write data
    fem_dataoutput.writeData( 1.0 /*dummy*/, "fem_solution" );
    // -------------------------------------------------------------

  }

  ErrorCalculator(&msfem_solution, &fem_solution).print(DSC_LOG_INFO_0);

} // function algorithm

} //namespace MsFEM {
} //namespace Multiscale {
} //namespace Dune {


