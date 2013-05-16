#ifdef HAVE_CMAKE_CONFIG
 #include "cmake_config.h"
#elif defined (HAVE_CONFIG_H)
 #include <config.h>
#endif // ifdef HAVE_CMAKE_CONFIG

#include "rigorous_msfem_solver.hh"

#include <dune/multiscale/common/righthandside_assembler.hh>
#include <dune/multiscale/msfem/localproblems/subgrid-list.hh>
#include <dune/multiscale/tools/misc/linear-lagrange-interpolation.hh>
#include <dune/multiscale/msfem/localproblems/localproblemsolver.hh>

#include <dune/common/fmatrix.hh>

#include <dune/fem/solver/oemsolver/oemsolver.hh>
#include <dune/fem/solver/inverseoperators.hh>
#include <dune/fem/operator/discreteoperatorimp.hh>
#include <dune/fem/operator/2order/lagrangematrixsetup.hh>
#include <dune/fem/operator/matrix/spmatrix.hh>
#include <dune/fem/space/common/adaptmanager.hh>


#include <dune/istl/matrix.hh>
#include <dune/stuff/fem/functions/checks.hh>

#include <dune/multiscale/msfem/msfem_grid_specifier.hh>

namespace Dune {
namespace Multiscale {
namespace MsFEM {

template< class MatrixImp >
void print_matrix( MatrixImp& system_matrix )
{
 std::cout << "---------------------------" << std::endl;
 std::cout << "Matrix:" << std::endl << std::endl;
 for (int row = 0; row != system_matrix.N(); ++row) {
   for (int col = 0; col != system_matrix.M(); ++col) {
     std::cout << system_matrix[row][col] << "  ";}
     std::cout << std::endl;
 }
 std::cout << "---------------------------" << std::endl;
 std::cout << std::endl << std::endl;
}

template< class VectorImp >
void print_vector( VectorImp& vector )
{
 std::cout << "---------------------------" << std::endl;
 std::cout << "Vector:" << std::endl << std::endl;
 for (int col = 0; col != vector.N(); ++col) {
     std::cout << vector[col] << "  ";
 }
 std::cout << std::endl << "---------------------------" << std::endl;
 std::cout << std::endl << std::endl;
}


Elliptic_Rigorous_MsFEM_Solver::Elliptic_Rigorous_MsFEM_Solver(const Elliptic_Rigorous_MsFEM_Solver::DiscreteFunctionSpace& discreteFunctionSpace)
: discreteFunctionSpace_(discreteFunctionSpace)
{}


// create a hostgrid function from a subgridfunction (projection for global continuity)
// Note: the maximum gride levels for both underlying grids must be the same
void Elliptic_Rigorous_MsFEM_Solver::subgrid_to_hostrid_projection(
        const Elliptic_Rigorous_MsFEM_Solver::SubgridDiscreteFunction& sub_func,
        Elliptic_Rigorous_MsFEM_Solver::DiscreteFunction& host_func) const
{
    host_func.clear();

    const SubgridDiscreteFunctionSpace& subDiscreteFunctionSpace = sub_func.space();
    const SubGridType& subGrid = subDiscreteFunctionSpace.grid();

    typedef typename SubgridDiscreteFunctionSpace::IteratorType SubgridIterator;
    typedef typename SubgridIterator::Entity                    SubgridEntity;
    typedef typename SubgridDiscreteFunction::LocalFunctionType SubgridLocalFunction;

    const SubgridIterator sub_endit = subDiscreteFunctionSpace.end();
    for (SubgridIterator sub_it = subDiscreteFunctionSpace.begin(); sub_it != sub_endit; ++sub_it)
    {
      const SubgridEntity& sub_entity = *sub_it;

      const HostEntityPointer host_entity_pointer = subGrid.getHostEntity< 0 >(*sub_it);
      const HostEntity& host_entity = *host_entity_pointer;

      const SubgridLocalFunction sub_loc_value = sub_func.localFunction(sub_entity);
      LocalFunction host_loc_value = host_func.localFunction(host_entity);

      const auto numBaseFunctions = sub_loc_value.baseFunctionSet().size();
      for (unsigned int i = 0; i < numBaseFunctions; ++i)
      {
        host_loc_value[i] = sub_loc_value[i];
      }
    }
} // subgrid_to_hostrid_projection

// vtk visualization of msfem basis functions
void Elliptic_Rigorous_MsFEM_Solver::vtk_output(
        Elliptic_Rigorous_MsFEM_Solver::MsFEMBasisFunctionType& msfem_basis_function_list,
        std::string basis_name ) const
{

  // general output parameters
  Dune::Multiscale::OutputParameters outputparam(DSC_CONFIG_GET("global.datadir", "data") + "/msfem_basis");

  typedef typename CommonTraits::IOTupleType IOTType;
  const auto& gridPart = msfem_basis_function_list[0]->space().gridPart();

  for ( size_t i = 0; i < msfem_basis_function_list.size(); i+=1 )
  {

    IOTType msfem_basis_series( &(*msfem_basis_function_list[i]) );

    const std::string ls_name_s = "/" + basis_name + (boost::format("_%d") % i).str();
    outputparam.set_prefix(ls_name_s);

    std::string outstring = basis_name;

    CommonTraits::DataOutputType msfem_basis_dataoutput(
      gridPart.grid(), msfem_basis_series, outputparam );
    msfem_basis_dataoutput.writeData( 1.0 /*dummy*/, outstring );

  }

  std::cout << "VTK Output for MsFEM basis functions successful." << std::endl << std::endl;

}



//! create standard coarse grid basis functions as discrete functions defined on the fine grid
// ------------------------------------------------------------------------------------
void Elliptic_Rigorous_MsFEM_Solver::add_coarse_basis_contribution(MacroMicroGridSpecifier& specifier,
                                    std::map<int,int>& global_id_to_internal_id,
                                    Elliptic_Rigorous_MsFEM_Solver::MsFEMBasisFunctionType& msfem_basis_function_list ) const
{

  DSC_LOG_INFO << "Create standard coarse grid basis functions as discrete functions on the fine grid... ";

  typedef typename HostEntity::Codim< 2 >::EntityPointer HostNodePointer;
  typedef typename GridPart::IntersectionIteratorType HostIntersectionIterator;
  typedef typename DiscreteFunctionSpace::BaseFunctionSetType CoarseBaseFunctionSet;

  const HostGridLeafIndexSet& coarseGridLeafIndexSet = specifier.coarseSpace().gridPart().grid().leafIndexSet();

  for (HostgridIterator it = discreteFunctionSpace_.begin(); it != discreteFunctionSpace_.end(); ++it)
  {
    typedef typename HostEntity::Codim< 0 >::EntityPointer
        HostEntityPointer;
    HostEntityPointer coarse_father = Stuff::Grid::make_father(coarseGridLeafIndexSet,
                                                               HostEntityPointer(*it),
                                                               specifier.getLevelDifference());

    const CoarseBaseFunctionSet coarseBaseSet = specifier.coarseSpace().baseFunctionSet( *coarse_father );
    const auto numBaseFunctions = coarseBaseSet.size();

    const auto& lagrangepoint_set = specifier.coarseSpace().lagrangePointSet(*coarse_father);
    const auto& coarse_geometry = (*coarse_father).geometry();

    const auto number_of_points = lagrangepoint_set.nop();

    std::vector< RangeType > phi( numBaseFunctions );
    //! TODO switch loops for more efficiency
    for(size_t loc_basis_number = 0; loc_basis_number < numBaseFunctions ; ++loc_basis_number ) {

      const int global_dof_number = specifier.coarseSpace().mapper().mapToGlobal(*coarse_father, loc_basis_number );
      if ( specifier.is_coarse_boundary_node( global_dof_number ) == true )
      { continue; }

      int global_interior_dof_number = global_id_to_internal_id[ global_dof_number ];

      // only implemented for 3 Lagrange Points, i.e. piecewise linear functions
      assert( number_of_points == 3 );
      std::vector< RangeType > phi_i( number_of_points );
      std::vector< DomainType > corners( number_of_points );

      for(size_t loc_point = 0; loc_point < number_of_points ; ++loc_point ) {

        coarseBaseSet.evaluateAll( lagrangepoint_set.point( loc_point ) , phi );
        phi_i[ loc_point ] = phi[ loc_basis_number ];
        corners[ loc_point ] = coarse_geometry.global(lagrangepoint_set.point( loc_point ) );
      }

      // LinearLagrangeInterpolation2D should be eventually replaced by
      // LinearLagrangeFunction2D< DiscreteFunctionSpace > coarse_basis_interpolation
      LinearLagrangeInterpolation2D< DiscreteFunctionSpace > coarse_basis_interpolation
          ( corners[0], phi_i[0], corners[1], phi_i[1], corners[2], phi_i[2] );

      LocalFunction loc_coarse_basis_function = (msfem_basis_function_list[global_interior_dof_number])->localFunction(*it);

      const int number_of_nodes_in_fine_entity = it->count< 2 >();
      if ( !( number_of_nodes_in_fine_entity == int( loc_coarse_basis_function.baseFunctionSet().size() ) ) )
      { DSC_LOG_ERROR << "Error! Inconsistency in 'rigorous_msfem_solver.hh'." << std::endl; }

      for (int i = 0; i < number_of_nodes_in_fine_entity; i += 1)
      {
        const HostNodePointer node = it->subEntity< 2 >(i);

        const DomainType coordinates_of_node = node->geometry().corner(0);
        if ( !( coordinates_of_node == it->geometry().corner(i) ) )
        { DSC_LOG_ERROR << "Error! Inconsistency in 'rigorous_msfem_solver.hh'." << std::endl; }

        RangeType coarse_value(0.0);
        coarse_basis_interpolation.evaluate(coordinates_of_node, coarse_value);

        // int global_index_node = gridPart.indexSet().index( *node );
        loc_coarse_basis_function[i] = coarse_value;
      }
    }
  }
  DSC_LOG_INFO << " done." << std::endl;
}

//! add corrector part to MsFEM basis functions
void Elliptic_Rigorous_MsFEM_Solver::add_corrector_contribution( MacroMicroGridSpecifier& specifier,
                                 std::map<int,int>& global_id_to_internal_id,
                                 MsFEMTraits::SubGridListType& subgrid_list,
                                 Elliptic_Rigorous_MsFEM_Solver::MsFEMBasisFunctionType& msfem_basis_function_list ) const
{

  DSC_LOG_INFO << "Add global corrector to create MsFEM basis functions from standard FEM basis functions... ";

  const HostGridLeafIndexSet& coarseGridLeafIndexSet = specifier.coarseSpace().gridPart().grid().leafIndexSet();

  typedef typename DiscreteFunctionSpace::IteratorType CoarseIterator;
  typedef typename CoarseIterator::Entity CoarseEntity;
  typedef typename CoarseEntity::Geometry CoarseGeometry;

  typedef typename DiscreteFunctionSpace::BaseFunctionSetType CoarseBaseFunctionSet;

  std::vector< JacobianRangeType > gradient_Phi(
     specifier.coarseSpace().mapper().maxNumDofs() );

  for (const CoarseEntity& coarse_grid_entity : specifier.coarseSpace())
  {

    const CoarseGeometry& coarse_grid_geometry = coarse_grid_entity.geometry();
    assert(coarse_grid_entity.partitionType() == InteriorEntity);

    const int global_index_entity = coarseGridLeafIndexSet.index(coarse_grid_entity);

    const CoarseBaseFunctionSet coarseBaseSet = specifier.coarseSpace().baseFunctionSet( coarse_grid_entity );
    const auto numBaseFunctions = coarseBaseSet.size();

    // the sub grid U(T) that belongs to the coarse_grid_entity T
    SubGridType& sub_grid_U_T = subgrid_list.getSubGrid(global_index_entity);
    SubGridPart subGridPart(sub_grid_U_T);

    const SubgridDiscreteFunctionSpace localDiscreteFunctionSpace(subGridPart);

    SubgridDiscreteFunction local_problem_solution_e0("Local problem Solution e_0", localDiscreteFunctionSpace);
    local_problem_solution_e0.clear();

    SubgridDiscreteFunction local_problem_solution_e1("Local problem Solution e_1", localDiscreteFunctionSpace);
    local_problem_solution_e1.clear();

    // --------- load local solutions -------
    // the file/place, where we saved the solutions of the cell problems
    const std::string local_solution_location = (boost::format("local_problems/_localProblemSolutions_%d")
                                                % global_index_entity).str();
    // reader for the cell problem data file:
    DiscreteFunctionReader discrete_function_reader(local_solution_location);
    // std::cout<< "... reading local problem solution " << global_index_entity << "/" << 0 << std::endl;
    discrete_function_reader.read(0, local_problem_solution_e0);
    // std::cout<< "... reading local problem solution " << global_index_entity << "/" << 1 << std::endl;
    discrete_function_reader.read(1, local_problem_solution_e1);

    // 1 point quadrature!! We only need the gradient of the base function,
    // which is constant on the whole entity.
    const CoarseQuadrature one_point_quadrature(coarse_grid_entity, 0);

    // the barycenter of the macro_grid_entity
    const typename CoarseQuadrature::CoordinateType& local_coarse_point
      = one_point_quadrature.point(0 /*=quadraturePoint*/);

    // transposed of the the inverse jacobian
    const auto& inverse_jac = coarse_grid_geometry.jacobianInverseTransposed(local_coarse_point);
    coarseBaseSet.jacobianAll(one_point_quadrature[0], inverse_jac, gradient_Phi);

    for (unsigned int i = 0; i < numBaseFunctions; ++i)
    {
      int global_dof_number = specifier.coarseSpace().mapper().mapToGlobal( coarse_grid_entity , i );
      if ( specifier.is_coarse_boundary_node( global_dof_number ) == true )
      { continue; }

      int global_interior_dof_number = global_id_to_internal_id[ global_dof_number ];

      DiscreteFunction correction_on_U_T("correction_on_U_T", discreteFunctionSpace_);

      correction_on_U_T.clear();
      subgrid_to_hostrid_projection(local_problem_solution_e0, correction_on_U_T);
      correction_on_U_T *= gradient_Phi[i][0][0];
      (*(msfem_basis_function_list[global_interior_dof_number])) += correction_on_U_T;

      correction_on_U_T.clear();
      subgrid_to_hostrid_projection(local_problem_solution_e1, correction_on_U_T);
      correction_on_U_T *= gradient_Phi[i][0][1];
      (*(msfem_basis_function_list[global_interior_dof_number])) += correction_on_U_T;

    }

  }

  DSC_LOG_INFO << " done." << std::endl;
}

void  Elliptic_Rigorous_MsFEM_Solver::solve_dirichlet_zero(const CommonTraits::DiffusionType& diffusion_op,
                          const CommonTraits::FirstSourceType& f,
                          // number of layers per coarse grid entity T:  U(T) is created by enrichting T with
                          // n(T)-layers.
                          MacroMicroGridSpecifier& specifier,
                          MsFEMTraits::SubGridListType& subgrid_list,
                          Elliptic_Rigorous_MsFEM_Solver::DiscreteFunction& coarse_scale_part,
                          Elliptic_Rigorous_MsFEM_Solver::DiscreteFunction& fine_scale_part,
                          Elliptic_Rigorous_MsFEM_Solver::DiscreteFunction& solution) const
{
  DSC::Profiler::ScopedTiming st("msfem.Elliptic_Rigorous_MsFEM_Solver.solve_dirichlet_zero");

  specifier.setOversamplingStrategy( 3 ); // for rigorous MsFEM!

  DiscreteFunctionSpace& coarse_space = specifier.coarseSpace();
  DiscreteFunctionSpace& fine_space = specifier.fineSpace();

  specifier.identify_coarse_boundary_nodes();

  const int number_of_internal_coarse_nodes = coarse_space.size() - specifier.get_number_of_coarse_boundary_nodes();

  // mapper: global_id_of_node -> new_id_of_node
  // ('new' means that we only count the internal nodes, boundary nodes do not receive an id)
  std::map<int,int> global_id_to_internal_id;
  for (int internal_id = 0, global_id = 0; global_id < coarse_space.size(); ++global_id)
  {
    if ( !specifier.is_coarse_boundary_node(global_id) )
    {
      global_id_to_internal_id[ global_id ] = internal_id;
      ++internal_id;
    }
  }

  MsFEMBasisFunctionType msfem_basis_function;
  for (int internal_id = 0; internal_id < number_of_internal_coarse_nodes; ++internal_id)
  {
    msfem_basis_function.emplace_back(new DiscreteFunction("MsFEM basis function", fine_space));
    msfem_basis_function[internal_id]->clear();
  }

  //! NOTE TODO for each MsFEM basis function save the support,
  //! i.e. a vector of entity points that describe the support of the basis
  //! function. This will save a lot of computational time when assembling the system matrix!

  MsFEMBasisFunctionType standard_basis_function;
  for (int internal_id = 0; internal_id < number_of_internal_coarse_nodes; internal_id += 1 )
   {
    standard_basis_function.emplace_back(new DiscreteFunction("Standard basis function", fine_space));
    standard_basis_function[internal_id]->clear();
   }

  add_coarse_basis_contribution( specifier, global_id_to_internal_id, msfem_basis_function );
  add_coarse_basis_contribution( specifier, global_id_to_internal_id, standard_basis_function );


  //! assemble all local problems (within constructor!)
  MsFEMLocalProblemSolver loc_prob_solver( specifier.fineSpace(), specifier, subgrid_list, diffusion_op,
                                           standard_basis_function, global_id_to_internal_id );
  loc_prob_solver.assemble_all(/*silence=*/false);

  // define the discrete (elliptic) operator that describes our problem
  add_corrector_contribution( specifier, global_id_to_internal_id, subgrid_list, msfem_basis_function );

  // just for VTK output for the basis function correctors
  /*
  MsFEMBasisFunctionType corrector_basis_function;
  for (int internal_id = 0; internal_id < number_of_internal_coarse_nodes; internal_id += 1 )
   {
    corrector_basis_function.emplace_back(new DiscreteFunction("Corrector basis function", fine_space));
    corrector_basis_function[internal_id]->clear();
   }

  add_corrector_contribution( specifier, global_id_to_internal_id, subgrid_list, corrector_basis_function );
  vtk_output( corrector_basis_function, "corrector_basis_function" );
  */

  if ( DSC_CONFIG_GET("rigorous_msfem.msfem_basis_vtk_output", 0) )
  {
     vtk_output( msfem_basis_function );
     vtk_output( standard_basis_function, "standard_basis_function" );
  }

  DSC_LOG_INFO << "Start assembling the stiffness matrix of the global problems.." << std::endl;

  DSC_LOG_INFO << "WARNING! Assembling the stiffness matrix of the global problems extremely expensive! Implementation is not yet efficient!" << std::endl;
  //! NOTE TODO for each MsFEM basis function save the support,
  //! i.e. a vector of entity points that describe the support of the basis
  //! function. This will save a lot of computational time when assembling the system matrix!

  //! (stiffness) matrix
  MatrixType system_matrix( number_of_internal_coarse_nodes, number_of_internal_coarse_nodes );

  if ( DSC_CONFIG_GET("rigorous_msfem.petrov_galerkin", true) )
  { assemble_matrix( diffusion_op, msfem_basis_function, standard_basis_function, system_matrix); }
  else
  { assemble_matrix( diffusion_op, msfem_basis_function, msfem_basis_function, system_matrix); }

  DSC_LOG_INFO << ".. assembling of the stiffness matrix done." << std::endl;

  //print_matrix( system_matrix );

  //! NOTE TODO: Assembling of right hand side is also quite expensive!
  VectorType rhs( number_of_internal_coarse_nodes );
  if ( DSC_CONFIG_GET("rigorous_msfem.petrov_galerkin", true) )
  { assemble_rhs( f, standard_basis_function, rhs ); }
  else
  { assemble_rhs( f, msfem_basis_function, rhs ); }

  //print_vector( rhs );

  double tol = DSC_CONFIG_GET("rigorous_msfem.macro_solver_tolerance", 1e-10 );
  int num_iterations = DSC_CONFIG_GET("rigorous_msfem.macro_solver_iterations", 10000 );

  MatrixOperatorType matrix_op( system_matrix );
  PreconditionerType preconditioner( system_matrix, 100, 0.9 );

  Dune::InverseOperatorResult result_data;
  VectorType solution_vector( number_of_internal_coarse_nodes );
  for (size_t col = 0; col != solution_vector.N(); ++col)
    solution_vector[col] = 0.0;

#ifdef SYMMETRIC_DIFFUSION_MATRIX
  if ( DSC_CONFIG_GET("rigorous_msfem.petrov_galerkin", true) )
  {
    typedef Dune::BiCGSTABSolver< VectorType > SolverType;
    SolverType solver( matrix_op, preconditioner, tol, num_iterations, true );
    solver.apply( solution_vector, rhs, result_data);
  }
  else
  {
    typedef Dune::CGSolver< VectorType > SolverType;
    SolverType solver( matrix_op, preconditioner, tol, num_iterations, true );
    solver.apply( solution_vector, rhs, result_data);
  }
#else
  typedef Dune::BiCGSTABSolver< VectorType > SolverType;

  SolverType solver( matrix_op, preconditioner, tol, num_iterations, true );
  solver.apply( solution_vector, rhs, result_data);
#endif

  coarse_scale_part.clear();
  for (int internal_id = 0; internal_id < number_of_internal_coarse_nodes; internal_id += 1 )
   {
     *(standard_basis_function[internal_id]) *= solution_vector[internal_id];
     coarse_scale_part += *(standard_basis_function[internal_id]);
   }

  solution.clear();
  for (int internal_id = 0; internal_id < number_of_internal_coarse_nodes; internal_id += 1 )
   {
     *(msfem_basis_function[internal_id]) *= solution_vector[internal_id];
     solution += *(msfem_basis_function[internal_id]);
   }

  fine_scale_part.assign(solution);
  fine_scale_part -= coarse_scale_part;
} // solve_dirichlet_zero

} //namespace MsFEM {
} //namespace Multiscale {
} //namespace Dune {


