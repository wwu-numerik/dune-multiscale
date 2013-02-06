#ifndef Elliptic_MSEM_Solver_HH
#define Elliptic_MSEM_Solver_HH

#include <dune/common/fmatrix.hh>

#include <dune/fem/solver/oemsolver/oemsolver.hh>
#include <dune/fem/solver/inverseoperators.hh>
#include <dune/fem/operator/discreteoperatorimp.hh>
#include <dune/fem/operator/2order/lagrangematrixsetup.hh>
#include <dune/fem/operator/matrix/spmatrix.hh>
#include <dune/fem/space/common/adaptmanager.hh>

#include <dune/multiscale/tools/assembler/righthandside_assembler.hh>
#include <dune/multiscale/tools/solver/MsFEM/msfem_localproblems/subgrid-list.hh>
#include <dune/multiscale/tools/assembler/matrix_assembler/elliptic_msfem_matrix_assembler.hh>
#include <dune/multiscale/tools/misc/linear-lagrange-interpolation.hh>
#include <dune/stuff/fem/functions/checks.hh>

#include <dune/multiscale/tools/solver/MsFEM/msfem_grid_specifier.hh>

// / done

namespace Dune {

template< class DiscreteFunctionType >
class Elliptic_Rigorous_MsFEM_Solver
{
public:
  typedef DiscreteFunctionType DiscreteFunction;

  typedef typename DiscreteFunction::FunctionSpaceType FunctionSpace;

  typedef typename DiscreteFunction::DiscreteFunctionSpaceType DiscreteFunctionSpace;

  typedef typename DiscreteFunction::LocalFunctionType LocalFunction;

  typedef typename DiscreteFunctionSpace::LagrangePointSetType LagrangePointSet;

  typedef typename DiscreteFunctionSpace::GridPartType GridPart;

  typedef typename DiscreteFunctionSpace::GridType HostGrid;

  typedef typename HostGrid::Traits::LeafIndexSet HostGridLeafIndexSet;

  typedef typename HostGrid::Traits::LeafIndexSet CoarseGridLeafIndexSet;

  typedef typename DiscreteFunctionSpace::DomainType        DomainType;
  typedef typename DiscreteFunctionSpace::RangeType         RangeType;
  typedef typename DiscreteFunctionSpace::JacobianRangeType JacobianRangeType;

  // typedef typename HostGrid ::template Codim< 0 > :: template Partition< All_Partition > :: LevelIterator
  // LevelEntityIteratorType;

  typedef typename DiscreteFunctionSpace::IteratorType HostgridIterator;

  typedef typename HostgridIterator::Entity HostEntity;

  typedef typename HostEntity::EntityPointer HostEntityPointer;

  // typedef typename HostGrid :: template Codim< 0 > :: template Partition< All_Partition > :: LevelIterator
  // HostGridLevelEntityIterator;

  enum { faceCodim = 1 };

  typedef typename GridPart::IntersectionIteratorType IntersectionIterator;

  typedef typename LagrangePointSet::template Codim< faceCodim >
    ::SubEntityIteratorType
  FaceDofIterator;

  // --------------------------- subgrid typedefs ------------------------------------

  typedef SubGrid< HostGrid::dimension, HostGrid > SubGridType;

  typedef LeafGridPart< SubGridType > SubGridPart;

  // typedef typename SubGridType ::template Codim< 0 > :: template Partition< All_Partition > :: LevelIterator
  // SubGridLevelEntityIteratorType;

  typedef LagrangeDiscreteFunctionSpace< FunctionSpace, SubGridPart, 1 >  // 1=POLORDER
  SubgridDiscreteFunctionSpace;

  typedef AdaptiveDiscreteFunction< SubgridDiscreteFunctionSpace > SubgridDiscreteFunction;

  typedef typename SubgridDiscreteFunctionSpace::IteratorType CoarseGridIterator;

  typedef typename CoarseGridIterator::Entity CoarseGridEntity;

  typedef typename CoarseGridEntity::EntityPointer CoarseGridEntityPointer;

  typedef typename SubgridDiscreteFunction::LocalFunctionType CoarseGridLocalFunction;

  typedef typename SubgridDiscreteFunctionSpace::LagrangePointSetType
  CoarseGridLagrangePointSet;

  typedef typename CoarseGridLagrangePointSet::template Codim< faceCodim >
    ::SubEntityIteratorType
  CoarseGridFaceDofIterator;

  //!-----------------------------------------------------------------------------------------

  //! --------------------- istl matrix and vector types -------------------------------------

  typedef BlockVector< FieldVector< double, 1> > VectorType;
  typedef Matrix< FieldMatrix< double,1,1 > > MatrixType;
  typedef MatrixAdapter< MatrixType, VectorType, VectorType > MatrixOperatorType;
  typedef SeqSSOR< MatrixType, VectorType, VectorType > PreconditionerType;
  typedef BiCGSTABSolver< VectorType > SolverType;
  typedef InverseOperatorResult InverseOperatorResultType;

  //! ----------------------------------------------------------------------------------------

private:
  const DiscreteFunctionSpace& discreteFunctionSpace_;

public:
  Elliptic_Rigorous_MsFEM_Solver(const DiscreteFunctionSpace& discreteFunctionSpace)
    : discreteFunctionSpace_(discreteFunctionSpace)
  {}

  template< class Stream >
  void oneLinePrint(Stream& stream, const DiscreteFunction& func) {
    typedef typename DiscreteFunction::ConstDofIteratorType
    DofIteratorType;
    DofIteratorType it = func.dbegin();
    stream << "\n" << func.name() << ": [ ";
    for ( ; it != func.dend(); ++it)
      stream << std::setw(5) << *it << "  ";

    stream << " ] " << std::endl;
  } // oneLinePrint

  template< class Stream >
  void oneLinePrint(Stream& stream, const SubgridDiscreteFunction& func) {
    typedef typename SubgridDiscreteFunction::ConstDofIteratorType
    DofIteratorType;
    DofIteratorType it = func.dbegin();
    stream << "\n" << func.name() << ": [ ";
    for ( ; it != func.dend(); ++it)
      stream << std::setw(5) << *it << "  ";

    stream << " ] " << std::endl;
  } // oneLinePrint

#if 0
  // create a hostgrid function from a subgridfunction (projection for global continuity)
  // Note: the maximum gride levels for both underlying grids must be the same
  void subgrid_to_hostrid_projection(const SubgridDiscreteFunction& sub_func, DiscreteFunction& host_func) const {
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

      const HostEntityPointer host_entity_pointer = subGrid.template getHostEntity< 0 >(*sub_it);
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

  
  //! copy coarse scale part of MsFEM solution into a function defined on the fine grid
  // ------------------------------------------------------------------------------------
  void identify_coarse_scale_part( MacroMicroGridSpecifier< DiscreteFunctionSpace >& specifier,
                                   const DiscreteFunction& coarse_msfem_solution,
                                   DiscreteFunction& coarse_scale_part ) const
  {
    
    DSC_LOG_INFO << "Indentifying coarse scale part of the MsFEM solution... ";
    
    coarse_scale_part.clear();
    typedef typename HostEntity::template Codim< 2 >::EntityPointer HostNodePointer;
    typedef typename GridPart::IntersectionIteratorType HostIntersectionIterator;
    const HostGridLeafIndexSet& coarseGridLeafIndexSet = coarse_msfem_solution.space().gridPart().grid().leafIndexSet();

    for (HostgridIterator it = discreteFunctionSpace_.begin(); it != discreteFunctionSpace_.end(); ++it)
    {
      typedef typename HostEntity::template Codim< 0 >::EntityPointer
      HostEntityPointer;
      HostEntityPointer coarse_father = Stuff::Grid::make_father(coarseGridLeafIndexSet,
                                                                 HostEntityPointer(*it),
                                                                 specifier.getLevelDifference());

      LinearLagrangeFunction2D< DiscreteFunctionSpace > interpolation_coarse(coarse_father);

      interpolation_coarse.set_corners(coarse_msfem_solution);

      LocalFunction host_loc_value = coarse_scale_part.localFunction(*it);

      const int number_of_nodes = (*it).template count< 2 >();
      if ( !( number_of_nodes == int( host_loc_value.baseFunctionSet().size() ) ) )
      { DSC_LOG_ERROR << "Error! Inconsistency in 'msfem_solver.hh'." << std::endl; }

      for (int i = 0; i < number_of_nodes; i += 1)
      {
        const HostNodePointer node = (*it).template subEntity< 2 >(i);

        const DomainType coordinates_of_node = node->geometry().corner(0);
        if ( !( coordinates_of_node == it->geometry().corner(i) ) )
        { DSC_LOG_ERROR << "Error! Inconsistency in 'msfem_solver.hh'." << std::endl; }

        RangeType coarse_value(0.0);
        interpolation_coarse.evaluate(coordinates_of_node, coarse_value);

        // int global_index_node = gridPart.indexSet().index( *node );
        host_loc_value[i] = coarse_value;
      }
    }
    DSC_LOG_INFO << " done." << std::endl;
  }
  // ------------------------------------------------------------------------------------

  


  //! identify fine scale part of MsFEM solution (including the projection!)
  // ------------------------------------------------------------------------------------
  template< class SubGridListType >
  void identify_fine_scale_part( MacroMicroGridSpecifier< DiscreteFunctionSpace >& specifier,
                                                          SubGridListType& subgrid_list,
                                                          const DiscreteFunction& coarse_msfem_solution,
                                                          DiscreteFunction& fine_scale_part ) const
  {

    fine_scale_part.clear();

    const HostGrid& grid = discreteFunctionSpace_.gridPart().grid();
    const GridPart& gridPart = discreteFunctionSpace_.gridPart();
    
    DiscreteFunctionSpace& coarse_space = specifier.coarseSpace();
    const HostGridLeafIndexSet& coarseGridLeafIndexSet = coarse_space.gridPart().grid().leafIndexSet();
    
    const int number_of_nodes = grid.size(2 /*codim*/);
    std::vector< std::vector< HostEntityPointer > > entities_sharing_same_node(number_of_nodes);

    // if the oversampling strategy is 1 or 2, we need identify the entities that share a certain node
    // (because a 'conforming projection' operator is required - for stragey 3, we directly get a conforming approximation)
    if ( ( specifier.getOversamplingStrategy() == 1 ) || ( specifier.getOversamplingStrategy() == 2 ) )
     {
       for (HostgridIterator it = discreteFunctionSpace_.begin(); it != discreteFunctionSpace_.end(); ++it)
       {
         const int number_of_nodes_in_entity = (*it).template count< 2 >();
         for (int i = 0; i < number_of_nodes_in_entity; i += 1)
         {
           const typename HostEntity::template Codim< 2 >::EntityPointer node = (*it).template subEntity< 2 >(i);
           const int global_index_node = gridPart.indexSet().index(*node);

           entities_sharing_same_node[global_index_node].push_back( HostEntityPointer(*it) );
         }
       }
     }

    DSC_LOG_INFO << "Indentifying fine scale part of the MsFEM solution... ";
    // iterator ueber coarse space
    for (HostgridIterator coarse_it = coarse_space.begin(); coarse_it != coarse_space.end(); ++coarse_it)
    {
      // the coarse entity 'T': *coarse_it

      // only required for oversampling strategy 1 and 2, where we need to identify the correction for each 
      DiscreteFunction correction_on_U_T("correction_on_U_T", discreteFunctionSpace_);
      correction_on_U_T.clear();

      const int index = coarseGridLeafIndexSet.index(*coarse_it);

      // the sub grid U(T) that belongs to the coarse_grid_entity T
      SubGridType& sub_grid_U_T = subgrid_list.getSubGrid(index);
      SubGridPart subGridPart(sub_grid_U_T);

      const SubgridDiscreteFunctionSpace localDiscreteFunctionSpace(subGridPart);

      SubgridDiscreteFunction local_problem_solution_e0("Local problem Solution e_0", localDiscreteFunctionSpace);
      local_problem_solution_e0.clear();

      SubgridDiscreteFunction local_problem_solution_e1("Local problem Solution e_1", localDiscreteFunctionSpace);
      local_problem_solution_e1.clear();

      // --------- load local solutions -------
      // the file/place, where we saved the solutions of the cell problems
      const std::string local_solution_location = (boost::format("local_problems/_localProblemSolutions_%d")
                                                  % index).str();
      // reader for the cell problem data file:
      DiscreteFunctionReader discrete_function_reader(local_solution_location);
      discrete_function_reader.read(0, local_problem_solution_e0);
      discrete_function_reader.read(1, local_problem_solution_e1);

      LocalFunction local_coarse_part = coarse_msfem_solution.localFunction(*coarse_it);

      // 1 point quadrature!! We only need the gradient of the coarse scale part on the element, which is a constant.
      const CachingQuadrature< GridPart, 0 > one_point_quadrature(*coarse_it, 0);

      JacobianRangeType grad_coarse_msfem_on_entity;
      local_coarse_part.jacobian(one_point_quadrature[0], grad_coarse_msfem_on_entity);

      //!
      // std :: cout << "grad_coarse_msfem_on_entity[ 0 ][ 1 ] = " << grad_coarse_msfem_on_entity[ 0 ][ 1 ] << std ::
      // endl;
      // std :: cout << "grad_coarse_msfem_on_entity[ 0 ][ 0 ] = " << grad_coarse_msfem_on_entity[ 0 ][ 0 ] << std ::
      // endl;
      // get the coarse gradient on T, multiply it with the local correctors and sum it up.
      local_problem_solution_e0 *= grad_coarse_msfem_on_entity[0][0];
      local_problem_solution_e1 *= grad_coarse_msfem_on_entity[0][1];
      local_problem_solution_e0 += local_problem_solution_e1;

      // oneLinePrint( DSC_LOG_DEBUG, local_problem_solution_e0 );
      
      // oversampling strategy 3: just sum up the local correctors:
      if ( (specifier.getOversamplingStrategy() == 3) )
       {
        subgrid_to_hostrid_projection(local_problem_solution_e0, correction_on_U_T);
       }

      // oversampling strategy 1 or 2: restrict the local correctors to the element T, sum them up and apply a conforming projection:
      if ( ( specifier.getOversamplingStrategy() == 1 ) || ( specifier.getOversamplingStrategy() == 2 ) )
       {
       
        if ( sub_grid_U_T.maxLevel() != discreteFunctionSpace_.gridPart().grid().maxLevel() )
        { DSC_LOG_ERROR << "Error: MaxLevel of SubGrid not identical to MaxLevel of FineGrid." << std::endl; }

        correction_on_U_T.clear();

        typedef typename SubgridDiscreteFunctionSpace::IteratorType SubgridIterator;
        typedef typename SubgridIterator::Entity                    SubgridEntity;
        typedef typename SubgridDiscreteFunction::LocalFunctionType SubgridLocalFunction;

        const SubgridIterator sub_endit = localDiscreteFunctionSpace.end();
        for (SubgridIterator sub_it = localDiscreteFunctionSpace.begin(); sub_it != sub_endit; ++sub_it)
        {
          const SubgridEntity& sub_entity = *sub_it;

          const HostEntityPointer fine_host_entity_pointer = sub_grid_U_T.template getHostEntity< 0 >(*sub_it);
          const HostEntity& fine_host_entity = *fine_host_entity_pointer;

          HostEntityPointer father = Stuff::Grid::make_father(coarseGridLeafIndexSet,
                                                              fine_host_entity_pointer,
                                                              specifier.getLevelDifference());
          if (!Stuff::Grid::entities_identical(*father,*coarse_it))
          {
            continue;
          }

          const SubgridLocalFunction sub_loc_value = local_problem_solution_e0.localFunction(sub_entity);
          LocalFunction host_loc_value = correction_on_U_T.localFunction(fine_host_entity);

          int number_of_nodes_entity = (*sub_it).template count< 2 >();
          for (int i = 0; i < number_of_nodes_entity; i += 1)
          {
            const typename HostEntity::template Codim< 2 >::EntityPointer node =
                fine_host_entity.template subEntity< 2 >(i);

            const int global_index_node = gridPart.indexSet().index(*node);

            // vector of coarse entities that share the above node
            std::vector< HostEntityPointer > coarse_entities;

            // count the number of different coarse-grid-entities that share the above node
            for (size_t j = 0; j < entities_sharing_same_node[global_index_node].size(); j += 1)
            {
              HostEntityPointer inner_it = Stuff::Grid::make_father(coarseGridLeafIndexSet,
                                                                    entities_sharing_same_node[global_index_node][j],
                                                                    specifier.getLevelDifference());
              bool new_entity_found = true;
              for (size_t k = 0; k < coarse_entities.size(); k += 1)
              {
                if (coarse_entities[k] == inner_it)
                { new_entity_found = false; }
              }
              if (new_entity_found)
              { coarse_entities.push_back(inner_it); }
            }

            host_loc_value[i] = ( sub_loc_value[i] / coarse_entities.size() );
          }
	}
       }
       
      fine_scale_part += correction_on_U_T;
    }
    DSC_LOG_INFO << " done." << std::endl;
  }
  // ------------------------------------------------------------------------------------
#endif

  //! create standard coarse grid basis functions as discrete functions defined on the fine grid
  // ------------------------------------------------------------------------------------
  template< class MsFEMBasisFunctionType >
  void add_coarse_basis_contribution( MacroMicroGridSpecifier< DiscreteFunctionSpace >& specifier,
				       std::map<int,int>& global_id_to_internal_id,
                                      MsFEMBasisFunctionType& msfem_basis_function_list ) const
  {
    
    DSC_LOG_INFO << "Create standard coarse grid basis functions as discrete functions on the fine grid... ";
    
    typedef typename HostEntity::template Codim< 2 >::EntityPointer HostNodePointer;
    typedef typename GridPart::IntersectionIteratorType HostIntersectionIterator;
    typedef typename DiscreteFunctionSpace::BaseFunctionSetType CoarseBaseFunctionSet;
  
    const HostGridLeafIndexSet& coarseGridLeafIndexSet = specifier.coarseSpace().gridPart().grid().leafIndexSet();

    for (HostgridIterator it = discreteFunctionSpace_.begin(); it != discreteFunctionSpace_.end(); ++it)
    {
      typedef typename HostEntity::template Codim< 0 >::EntityPointer
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
      //! TODO swich loops for more efficiency
      for(int loc_basis_number = 0; loc_basis_number < numBaseFunctions ; ++loc_basis_number ) {
      
        int global_dof_number = specifier.coarseSpace().mapToGlobal(*coarse_father, loc_basis_number );
	if ( specifier.is_coarse_boundary_node( global_dof_number ) == true )
	{ continue; }

        int global_interior_dof_number = global_id_to_internal_id[ global_dof_number ];
 
	// only implemented for 3 Lagrange Points, i.e. piecewise linear functions
	assert( number_of_points == 3 );
        std::vector< RangeType > phi_i( number_of_points );
        std::vector< DomainType > corners( number_of_points );
	
        for(int loc_point = 0; loc_point < number_of_points ; ++loc_point ) {

          coarseBaseSet.evaluateAll( lagrangepoint_set.point( loc_point ) , phi );
          phi_i[ loc_point ] = phi[ loc_basis_number ]; 
          corners[ loc_point ] = coarse_geometry.global(lagrangepoint_set.point( loc_point ) );
        }

        LinearLagrangeFunction2D< DiscreteFunctionSpace > coarse_basis_interpolation
          ( corners[0], phi_i[0], corners[1], phi_i[1], corners[2], phi_i[2] );
       
        LocalFunction loc_coarse_basis_function = (msfem_basis_function_list[global_interior_dof_number])->localFunction(*it);
      
        const int number_of_nodes_in_fine_entity = (*it).template count< 2 >();
        if ( !( number_of_nodes_in_fine_entity == int( loc_coarse_basis_function.baseFunctionSet().size() ) ) )
         { DSC_LOG_ERROR << "Error! Inconsistency in 'rigorous_msfem_solver.hh'." << std::endl; }
      
        for (int i = 0; i < number_of_nodes_in_fine_entity; i += 1)
         {
           const HostNodePointer node = (*it).template subEntity< 2 >(i);

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
  // ------------------------------------------------------------------------------------



  // - ∇ (A(x,∇u)) + b ∇u + c u = f - divG
  // then:
  // A --> diffusion operator ('DiffusionOperatorType')
  // b --> advective part ('AdvectionTermType')
  // c --> reaction part ('ReactionTermType')
  // f --> 'first' source term, scalar ('SourceTermType')
  // G --> 'second' source term, vector valued ('SecondSourceTermType')

  // homogenous Dirchilet boundary condition!:
  template< class DiffusionOperator, class SourceTerm, class SubGridListType >
  void solve_dirichlet_zero(const DiffusionOperator& diffusion_op,
                            const SourceTerm& f,
                            // number of layers per coarse grid entity T:  U(T) is created by enrichting T with
                            // n(T)-layers.
                            MacroMicroGridSpecifier< DiscreteFunctionSpace >& specifier,
                            SubGridListType& subgrid_list,
                            DiscreteFunction& coarse_scale_part,
                            DiscreteFunction& fine_scale_part,
                            DiscreteFunction& solution) const
  {
    DSC::Profiler::ScopedTiming st("msfem.Elliptic_Rigorous_MsFEM_Solver.solve_dirichlet_zero");
    
    specifier.setOversamplingStrategy( 3 ); // for rigorous MsFEM!
   
    DiscreteFunctionSpace& coarse_space = specifier.coarseSpace();
    DiscreteFunctionSpace& fine_space = specifier.fineSpace();
    
    specifier.identify_coarse_boundary_nodes();

    int number_of_internal_coarse_nodes = coarse_space.size() - specifier.get_number_of_coarse_boundary_nodes();

    // mapper: global_id_of_node -> new_id_of_node
    // ('new' means that we only count the internal nodes, boundary nodes do not receive an id)
    std::map<int,int> global_id_to_internal_id;
    int internal_id = 0;
    for (int global_id = 0; global_id < coarse_space.size(); global_id += 1 )
    {
      if ( specifier.is_coarse_boundary_node(global_id) == false )
      {
	global_id_to_internal_id[ global_id ] = internal_id;
	internal_id += 1;
      }
    }

    typedef std::vector< shared_ptr<DiscreteFunction> > MsFEMBasisFunctionType;
    MsFEMBasisFunctionType msfem_basis_function;
    for (int internal_id = 0; internal_id < number_of_internal_coarse_nodes; internal_id += 1 )
      msfem_basis_function.push_back(make_shared<DiscreteFunction>("MsFEM basis function", fine_space));

    add_coarse_basis_contribution( specifier, global_id_to_internal_id, msfem_basis_function );

    // VTK Output...

    // 0. Set up ID Mapper (identifiziere nicht-boundary LagrangePoints und verpasse diesen eine eigene Nummerierung)
    //    (in zwei Walks erledigen: 1. Walk: nur die Boundary Nodes sammeln, 2. Walk: den Rest identifizieren)
    //    eventuell einen Walk (oder sogar beide) in den Specifier stecken
    // 1. pointer array aus discrete functions auf dem feinen Gitter
    // 2. Befuelle diese mit den Macro Gitter Basis Funktionen
    // 3. addiere die globalen Korrektoren aufaddieren
    // 4. assembliere und loese das LGS
    
#if 0
    // discrete elliptic MsFEM operator (corresponds with MsFEM Matrix)
    typedef DiscreteEllipticMsFEMOperator< DiscreteFunction /*type of coarse space*/,
                                           MacroMicroGridSpecifier< DiscreteFunctionSpace >,
                                           DiscreteFunction /*type of fine space*/,
                                           DiffusionOperator > EllipticMsFEMOperatorType;



    // ------------------------------------------------------------

    DiscreteFunction coarse_msfem_solution("Coarse Part MsFEM Solution", coarse_space);
    coarse_msfem_solution.clear();

    //! define the right hand side assembler tool
    // (for linear and non-linear elliptic and parabolic problems, for sources f and/or G )
    typedef RightHandSideAssembler< DiscreteFunction > RhsAssembler;

    //! define the discrete (elliptic) operator that describes our problem
    // ( effect of the discretized differential operator on a certain discrete function )
    const EllipticMsFEMOperatorType elliptic_msfem_op(specifier,
                                                coarse_space,
                                                subgrid_list,
                                                diffusion_op);
    // discrete elliptic operator (corresponds with FEM Matrix)

    //! (stiffness) matrix
    MsFEMMatrix msfem_matrix("MsFEM stiffness matrix", coarse_space, coarse_space);

    //! right hand side vector
    // right hand side for the finite element method:
    DiscreteFunction msfem_rhs("MsFEM right hand side", coarse_space);
    msfem_rhs.clear();

    DSC_LOG_INFO  << std::endl << "Solving MsFEM problem." << std::endl;
    DSC_LOG_INFO << "Solving linear problem with MsFEM and maximum coarse grid level "
                << coarse_space.gridPart().grid().maxLevel() << "." << std::endl;
    DSC_LOG_INFO << "------------------------------------------------------------------------------" << std::endl;

    // to assemble the computational time
    Dune::Timer assembleTimer;

    // assemble the MsFEM stiffness matrix
    elliptic_msfem_op.assemble_matrix(msfem_matrix);   // einbinden!
    DSC_LOG_INFO << "Time to assemble MsFEM stiffness matrix: " << assembleTimer.elapsed() << "s" << std::endl;

    // assemble right hand side
    if ( DSC_CONFIG_GET("msfem.petrov_galerkin", 1 ) )
    { RhsAssembler::template assemble< 2* DiscreteFunctionSpace::polynomialOrder + 2 >(f, msfem_rhs); }
    else
    { RhsAssembler::template assemble_for_MsFEM_symmetric< 2* DiscreteFunctionSpace::polynomialOrder + 2 >(f, specifier, subgrid_list, msfem_rhs); }

    // oneLinePrint( DSC_LOG_DEBUG, fem_rhs );

    //! --- boundary treatment ---
    // set the dirichlet points to zero (in righ hand side of the fem problem)
    const HostgridIterator endit = coarse_space.end();
    for (HostgridIterator it = coarse_space.begin(); it != endit; ++it)
    {
      IntersectionIterator iit = coarse_space.gridPart().ibegin(*it);
      const IntersectionIterator endiit = coarse_space.gridPart().iend(*it);
      for ( ; iit != endiit; ++iit)
      {
        if ( !(*iit).boundary() )
          continue;

        LocalFunction rhsLocal = msfem_rhs.localFunction(*it);

        const LagrangePointSet& lagrangePointSet
          = coarse_space.lagrangePointSet(*it);

        const int face = (*iit).indexInInside();

        FaceDofIterator faceIterator
          = lagrangePointSet.template beginSubEntity< faceCodim >(face);
        const FaceDofIterator faceEndIterator
          = lagrangePointSet.template endSubEntity< faceCodim >(face);
        for ( ; faceIterator != faceEndIterator; ++faceIterator)
          rhsLocal[*faceIterator] = 0;
      }
    }
    //! --- end boundary treatment ---

    const InverseMsFEMMatrix msfem_biCGStab(msfem_matrix, 1e-8, 1e-8, 20000, true);
    msfem_biCGStab(msfem_rhs, coarse_msfem_solution);

    DSC_LOG_INFO << "---------------------------------------------------------------------------------" << std::endl;
    DSC_LOG_INFO << "MsFEM problem solved in " << assembleTimer.elapsed() << "s." << std::endl << std::endl
                << std::endl;

    // oneLinePrint( DSC_LOG_DEBUG, solution );
    // copy coarse grid function (defined on the subgrid) into a fine grid function
    solution.clear();

    //! copy coarse scale part of MsFEM solution into a function defined on the fine grid
    identify_coarse_scale_part( specifier, coarse_msfem_solution, coarse_scale_part );

    //! identify fine scale part of MsFEM solution (including the projection!)
    identify_fine_scale_part( specifier, subgrid_list, coarse_msfem_solution, fine_scale_part );

    // Auf Grobskalen MsFEM Anteil noch Feinksalen MsFEM Anteil aufaddieren.
    solution += coarse_scale_part;
    solution += fine_scale_part;
#endif
    

/// TEST AREA
#if 0
  

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
  

  MatrixOperatorType matrix_op( system_matrix );


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
    
    std::vector< MsfemTraits::DomainType > lagrange_points( number_of_points );
    for(int loc_point = 0; loc_point < number_of_points ; ++loc_point ) {
      int global_dof_number = discreteFunctionSpace.mapToGlobal(entity, loc_point );
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

#endif

    std::abort();
    
  } // solve_dirichlet_zero

  //! the following methods are not yet implemented, however note that the required tools are
  //! already available via 'righthandside_assembler.hh' and 'elliptic_fem_matrix_assembler.hh'!

  template< class DiffusionOperatorType, class ReactionTermType, class SourceTermType >
  void solve() {
    DSC_LOG_ERROR << "No implemented!" << std::endl;
  }

  template< class DiffusionOperatorType, class ReactionTermType, class SourceTermType, class SecondSourceTermType >
  void solve() {
    DSC_LOG_ERROR << "No implemented!" << std::endl;
  }
};
}

#endif // #ifndef Elliptic_MSEM_Solver_HH
