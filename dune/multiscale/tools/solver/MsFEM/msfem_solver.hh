#ifndef Elliptic_MSEM_Solver_HH
#define Elliptic_MSEM_Solver_HH

#include <dune/common/fmatrix.hh>

#include <dune/subgrid/subgrid.hh>


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

namespace Dune
{
  
  class MacroMicroGridSpecifier
  {

  public:
    
    MacroMicroGridSpecifier( int& number_of_level_host_entities, int& coarse_level_fine_level_difference )
    : coarse_level_fine_level_difference_( coarse_level_fine_level_difference ),
      number_of_level_host_entities_( number_of_level_host_entities )
     {
       for ( int i = 0; i < number_of_level_host_entities; i+=1 )
         {
	   // initialize with 0 layers:
	   number_of_layers.push_back(0);  
	 }
     }
     
    // get number of coarse grid entities
    int getNumOfCoarseEntities()
     {
       return number_of_level_host_entities_;
     }
     
    void setLayer( int i , int number_of_layers_for_entity )
     {
       if ( i < number_of_level_host_entities_ )
        { number_of_layers[ i ] = number_of_layers_for_entity; }
       else
        { std :: cout << "Error. Assertion (i < number_of_level_host_entities_) not filfilled." << std :: endl; abort(); }
     }

    int getLayer( int i )
     {
       if ( i < number_of_level_host_entities_ )
        { return number_of_layers[ i ]; }
       else
        { std :: cout << "Error. Assertion (i < number_of_level_host_entities_) not filfilled." << std :: endl; abort(); }
       return 0;
     }
     
  private:
    
    // level difference bettween coarse grid level and fine grid level
    int coarse_level_fine_level_difference_;
    
    // nomber of coarse grid entities
    int number_of_level_host_entities_;
    
    // layers for each coarse grid entity
    std :: vector < int > number_of_layers;
    
  };
  

  template< class DiscreteFunctionType >
  class Elliptic_MsFEM_Solver
  {

  public:

   typedef DiscreteFunctionType  DiscreteFunction;

   typedef typename DiscreteFunction :: FunctionSpaceType FunctionSpace;

   typedef typename DiscreteFunction :: DiscreteFunctionSpaceType DiscreteFunctionSpace;

   typedef typename DiscreteFunction :: LocalFunctionType LocalFunction;

   typedef typename DiscreteFunctionSpace :: LagrangePointSetType LagrangePointSet;

   typedef typename DiscreteFunctionSpace :: GridPartType GridPart;

   typedef typename DiscreteFunctionSpace :: GridType HostGrid;
   
   typedef typename DiscreteFunctionSpace :: DomainType DomainType;
   typedef typename DiscreteFunctionSpace :: RangeType RangeType;
   typedef typename DiscreteFunctionSpace :: JacobianRangeType JacobianRangeType;
   
   typedef typename HostGrid ::template Codim< 0 > :: template Partition< All_Partition > :: LevelIterator LevelEntityIteratorType;

   typedef typename DiscreteFunctionSpace :: IteratorType HostgridIterator;
   
   typedef typename HostgridIterator :: Entity HostEntity;
   
   typedef typename HostEntity :: EntityPointer HostEntityPointer; 

   typedef typename HostGrid :: template Codim< 0 > :: template Partition< All_Partition > :: LevelIterator HostGridLevelEntityIterator;

   enum { faceCodim = 1 };

   typedef typename GridPart :: IntersectionIteratorType IntersectionIterator;

   typedef typename LagrangePointSet :: template Codim< faceCodim > 
                                     :: SubEntityIteratorType
    FaceDofIterator;

   // --------------------------- subgrid typedefs ------------------------------------

   typedef SubGrid< HostGrid::dimension , HostGrid > SubGridType; 

   typedef LeafGridPart< SubGridType > SubGridPart;

   typedef typename SubGridType ::template Codim< 0 > :: template Partition< All_Partition > :: LevelIterator SubGridLevelEntityIteratorType;
   
   typedef LagrangeDiscreteFunctionSpace < FunctionSpace, SubGridPart, 1 > //1=POLORDER
      SubgridDiscreteFunctionSpace;

   typedef AdaptiveDiscreteFunction < SubgridDiscreteFunctionSpace > SubgridDiscreteFunction;
   
   typedef typename SubgridDiscreteFunctionSpace :: IteratorType CoarseGridIterator;
   
   typedef typename CoarseGridIterator :: Entity CoarseGridEntity;
   
   typedef typename CoarseGridEntity :: EntityPointer CoarseGridEntityPointer; 


   typedef typename SubgridDiscreteFunction :: LocalFunctionType CoarseGridLocalFunction;
   
   typedef typename SubgridDiscreteFunctionSpace :: LagrangePointSetType
    CoarseGridLagrangePointSet;
   
   typedef typename CoarseGridLagrangePointSet :: template Codim< faceCodim > 
                                               :: SubEntityIteratorType
    CoarseGridFaceDofIterator;
//!-----------------------------------------------------------------------------------------


   //! --------------------- the standard matrix traits -------------------------------------

   struct MatrixTraits
   {
     typedef SubgridDiscreteFunctionSpace RowSpaceType;
     typedef SubgridDiscreteFunctionSpace ColumnSpaceType;
     typedef LagrangeMatrixSetup< false > StencilType;
     typedef ParallelScalarProduct< SubgridDiscreteFunctionSpace > ParallelScalarProductType;

     template< class M >
     struct Adapter
     {
       typedef LagrangeParallelMatrixAdapter< M > MatrixAdapterType;
     };
   };

   //! --------------------------------------------------------------------------------------


   //! --------------------- type of fem stiffness matrix -----------------------------------

   typedef SparseRowMatrixOperator< SubgridDiscreteFunction, SubgridDiscreteFunction, MatrixTraits > MsFEMMatrix;

   //! --------------------------------------------------------------------------------------


   //! --------------- solver for the linear system of equations ----------------------------

   // use Bi CG Stab [OEMBICGSTABOp] or GMRES [OEMGMRESOp] for non-symmetric matrices and CG [CGInverseOp] for symmetric ones. GMRES seems to be more stable, but is extremely slow!
   typedef OEMBICGSQOp/*OEMBICGSTABOp*/< SubgridDiscreteFunction, MsFEMMatrix > InverseMsFEMMatrix;

   //! --------------------------------------------------------------------------------------


  private:
    const DiscreteFunctionSpace &discreteFunctionSpace_;

    std :: ofstream *data_file_;
    
    // path where to save the data output
    std :: string path_;

    
  public:
   Elliptic_MsFEM_Solver( const DiscreteFunctionSpace &discreteFunctionSpace, std :: string path = "" )
     : discreteFunctionSpace_( discreteFunctionSpace ),
       data_file_( NULL )
     { path_ = path; }

   Elliptic_MsFEM_Solver( const DiscreteFunctionSpace &discreteFunctionSpace, std :: ofstream& data_file, std :: string path = "" )
     : discreteFunctionSpace_( discreteFunctionSpace ),
       data_file_( &data_file )
     { path_ = path; }
     
   template < class Stream >
   void oneLinePrint( Stream& stream, const DiscreteFunction& func )
    {
      typedef typename DiscreteFunction::ConstDofIteratorType
         DofIteratorType;
      DofIteratorType it = func.dbegin();
      stream << "\n" << func.name() << ": [ ";
      for ( ; it != func.dend(); ++it )
         stream << std::setw(5) << *it << "  ";

      stream << " ] " << std::endl;
     }

   template < class Stream >
   void oneLinePrint( Stream& stream, const SubgridDiscreteFunction& func )
    {
      typedef typename SubgridDiscreteFunction::ConstDofIteratorType
         DofIteratorType;
      DofIteratorType it = func.dbegin();
      stream << "\n" << func.name() << ": [ ";
      for ( ; it != func.dend(); ++it )
         stream << std::setw(5) << *it << "  ";

      stream << " ] " << std::endl;
     }

   // create a hostgrid function from a subgridfunction (projection for global continuity)
   // Note: the maximum gride levels for both underlying grids must be the same
   void subgrid_to_hostrid_projection( const SubgridDiscreteFunction &sub_func, DiscreteFunction &host_func )
    {

       if ( sub_func.space().gridPart().grid().maxLevel() != host_func.space().gridPart().grid().maxLevel() )
         { std :: cout << "Error in method 'subgrid_to_hostrid_function': MaxLevel of SubGrid not identical to MaxLevel of FineGrid." << std :: endl; }

       host_func.clear();

       const SubgridDiscreteFunctionSpace &subDiscreteFunctionSpace = sub_func.space();
       const SubGridType &subGrid = subDiscreteFunctionSpace.grid();

       typedef typename SubgridDiscreteFunctionSpace :: IteratorType SubgridIterator;
       typedef typename SubgridIterator :: Entity SubgridEntity;
       typedef typename SubgridDiscreteFunction :: LocalFunctionType SubgridLocalFunction;

       SubgridIterator sub_endit = subDiscreteFunctionSpace.end();
       for( SubgridIterator sub_it = subDiscreteFunctionSpace.begin(); sub_it != sub_endit; ++sub_it )
          {

             const SubgridEntity &sub_entity = *sub_it;

             HostEntityPointer host_entity_pointer = subGrid.template getHostEntity<0>( *sub_it );
             const HostEntity& host_entity = *host_entity_pointer;

             SubgridLocalFunction sub_loc_value = sub_func.localFunction( sub_entity );
             LocalFunction host_loc_value = host_func.localFunction( host_entity );

             const unsigned int numBaseFunctions = sub_loc_value.baseFunctionSet().numBaseFunctions();
             for( unsigned int i = 0; i < numBaseFunctions; ++i )
               {
                 host_loc_value[ i ] = sub_loc_value[ i ];
               }

          }

    }



   // - ∇ (A(x,∇u)) + b ∇u + c u = f - divG
   // then:
   // A --> diffusion operator ('DiffusionOperatorType')
   // b --> advective part ('AdvectionTermType')
   // c --> reaction part ('ReactionTermType')
   // f --> 'first' source term, scalar ('SourceTermType')
   // G --> 'second' source term, vector valued ('SecondSourceTermType')


   // homogenous Dirchilet boundary condition!:
   template< class DiffusionOperator, class SourceTerm >
   void solve_dirichlet_zero( const DiffusionOperator &diffusion_op,
                              const SourceTerm &f,
                                    DiscreteFunctionSpace& coarse_space,
                                    DiscreteFunction& coarse_scale_part,
                                    DiscreteFunction& fine_scale_part,
                                    DiscreteFunction &solution )
   {
     HostGrid &grid = discreteFunctionSpace_.gridPart().grid();
     const GridPart &gridPart = discreteFunctionSpace_.gridPart();
     int number_of_level_host_entities = grid.size( coarse_space.gridPart().grid().maxLevel(), 0 /*codim*/ );
     int coarse_level_fine_level_difference = grid.maxLevel() - coarse_space.grid().maxLevel();
     MacroMicroGridSpecifier specifier( coarse_space.size() , discreteFunctionSpace_.gridPart().grid().maxLevel() - coarse_space.grid().maxLevel() );
     
     //! default - no layers
     // number of layers per coarse grid entity T:  U(T) is created by enrichting T with n(T)-layers.
     for ( int i = 0; i < number_of_level_host_entities; i+=1 )
       { specifier.setLayer( i , 0 ); }
    
     solve_dirichlet_zero( diffusion_op, f, coarse_space, specifier, coarse_scale_part, fine_scale_part, solution );
   }

   
   // homogenous Dirchilet boundary condition!:
   template< class DiffusionOperator, class SourceTerm >
   void solve_dirichlet_zero( const DiffusionOperator &diffusion_op,
                              const SourceTerm &f,
                              DiscreteFunctionSpace& coarse_space,
                              // number of layers per coarse grid entity T:  U(T) is created by enrichting T with n(T)-layers.
                              MacroMicroGridSpecifier& specifier,
                              DiscreteFunction& coarse_scale_part,
                              DiscreteFunction& fine_scale_part,
                              DiscreteFunction& solution )
   {
     
     // discrete elliptic MsFEM operator (corresponds with MsFEM Matrix)
     typedef DiscreteEllipticMsFEMOperator< SubgridDiscreteFunction, MacroMicroGridSpecifier,
                                            DiscreteFunction,
                                            DiffusionOperator > EllipticMsFEMOperatorType;

     int coarse_level = coarse_space.gridPart().grid().maxLevel(); 

     HostGrid &grid = discreteFunctionSpace_.gridPart().grid();
     const GridPart &gridPart = discreteFunctionSpace_.gridPart();

     // create subgrid:
     SubGridType subGrid( coarse_space.gridPart().grid() );
     subGrid.createBegin();

     //subGrid.insertLevel(coarse_level);
     //!!!! Diesen Iterator durch einen Iterator uber den SubSpace ersetzen??? (Klappt dann die Entity Identifikation?)
     LevelEntityIteratorType coarse_level_it = grid.template lbegin< 0 >( coarse_level );
     for( ; coarse_level_it != grid.template lend< 0 >( coarse_level ); ++coarse_level_it )
         subGrid.insertPartial( *coarse_level_it );

     subGrid.createEnd();

     subGrid.report();

     SubGridPart subGridPart( subGrid );

     SubgridDiscreteFunctionSpace coarseDiscreteFunctionSpace( subGridPart );

     // ----- check index sets --------------------

#if 1
      typedef typename HostGrid :: Traits :: LevelIndexSet HostGridLevelIndexSet;
      typedef typename SubGridType :: Traits :: LevelIndexSet SubGridLevelIndexSet;

      HostGridLevelEntityIterator level_iterator_end = grid.template lend< 0 >( coarse_level );
      HostGridLevelEntityIterator level_iterator_begin = grid.template lbegin< 0 >( coarse_level );
      
      const HostGridLevelIndexSet& hostGridLevelIndexSet = grid.levelIndexSet( coarse_level );

      const SubGridLevelIndexSet& subGridLevelIndexSet = subGrid.levelIndexSet( coarse_level );

      bool index_sets_correct = true;

      CoarseGridIterator coarse_it_ = coarseDiscreteFunctionSpace.begin();
      for ( HostGridLevelEntityIterator lit = level_iterator_begin;
            lit != level_iterator_end ; ++lit )
       {
          int index_host = hostGridLevelIndexSet.index( *lit );

          int index_sub = coarseDiscreteFunctionSpace.gridPart().indexSet().index( *coarse_it_ );

          if ( index_host != index_sub )
           { index_sets_correct = false; }

          int number_of_nodes = (*lit).template count<2>();

          for ( int i = 0; i < number_of_nodes; i += 1 )
           {
             if ( (*lit).geometry().corner(i) != (*coarse_it_).geometry().corner(i) )
              { index_sets_correct = false; }
           }
         ++coarse_it_;
       }

     if ( index_sets_correct == false )
      { std :: cout << "Index Sets not correct" << std :: endl; abort(); }
#endif

     // -------------------------------------------

     SubgridDiscreteFunction coarse_msfem_solution( "Coarse Part MsFEM Solution", coarseDiscreteFunctionSpace );
     coarse_msfem_solution.clear();

     //! create subgrids:
     bool silence = false;

     typedef SubGridList< DiscreteFunction, SubGridType, MacroMicroGridSpecifier > SubGridListType;
     SubGridListType subgrid_list( discreteFunctionSpace_ , specifier, coarse_level , silence );


     //! define the right hand side assembler tool
     // (for linear and non-linear elliptic and parabolic problems, for sources f and/or G )
     RightHandSideAssembler< SubgridDiscreteFunction > rhsassembler;

     //! define the discrete (elliptic) operator that describes our problem
     // ( effect of the discretized differential operator on a certain discrete function )
     EllipticMsFEMOperatorType elliptic_msfem_op( coarseDiscreteFunctionSpace,
                                                  discreteFunctionSpace_,
                                                  subgrid_list,
                                                  diffusion_op, *data_file_, path_ );
     // discrete elliptic operator (corresponds with FEM Matrix)

     //! (stiffness) matrix
     MsFEMMatrix msfem_matrix( "MsFEM stiffness matrix", coarseDiscreteFunctionSpace, coarseDiscreteFunctionSpace );

     //! right hand side vector
     // right hand side for the finite element method:
     SubgridDiscreteFunction msfem_rhs( "MsFEM right hand side", coarseDiscreteFunctionSpace );
     msfem_rhs.clear();

     std :: cout << "Solving MsFEM problem." << std :: endl;

     if ( data_file_ )
      {
        if (data_file_->is_open())
         {
           *data_file_ << "Solving linear problem with MsFEM and coarse grid level " << coarse_level << "." << std :: endl;
           *data_file_ << "------------------------------------------------------------------------------" << std :: endl;
         }
      }
      
     // to assemble the computational time
     Dune::Timer assembleTimer;
   
     // assemble the MsFEM stiffness matrix
     elliptic_msfem_op.assemble_matrix( msfem_matrix ); // einbinden!

     std::cout << "Time to assemble MsFEM stiffness matrix: " << assembleTimer.elapsed() << "s" << std::endl;

     if ( data_file_ )
      {
        if (data_file_->is_open())
         {
           *data_file_ << "Time to assemble MsFEM stiffness matrix: " << assembleTimer.elapsed() << "s" << std::endl;
         }
      }
      
     // assemble right hand side
     rhsassembler.template assemble< 2 * SubgridDiscreteFunctionSpace :: polynomialOrder + 2 >( f , msfem_rhs);

     //oneLinePrint( std::cout , fem_rhs );
    
     
     // --- boundary treatment ---
     // set the dirichlet points to zero (in righ hand side of the fem problem)
     CoarseGridIterator endit = coarseDiscreteFunctionSpace.end();
     for( CoarseGridIterator it = coarseDiscreteFunctionSpace.begin(); it != endit; ++it )
       {

         HostEntityPointer host_entity_pointer = subGrid.template getHostEntity<0>( *it );
         const HostEntity& host_entity = *host_entity_pointer;

         IntersectionIterator iit = gridPart.ibegin( host_entity );
         const IntersectionIterator endiit = gridPart.iend( host_entity );
         for( ; iit != endiit; ++iit )
           {

              if( !(*iit).boundary() )
                continue;

              CoarseGridLocalFunction rhsLocal = msfem_rhs.localFunction( *it );
	      
              const CoarseGridLagrangePointSet &lagrangePointSet
                = coarseDiscreteFunctionSpace.lagrangePointSet( *it );

              const int face = (*iit).indexInInside();

              CoarseGridFaceDofIterator faceIterator
                = lagrangePointSet.template beginSubEntity< faceCodim >( face );
              const CoarseGridFaceDofIterator faceEndIterator
                = lagrangePointSet.template endSubEntity< faceCodim >( face );
              for( ; faceIterator != faceEndIterator; ++faceIterator )
                rhsLocal[ *faceIterator ] = 0;

           }

       }
     // --- end boundary treatment ---

     InverseMsFEMMatrix msfem_biCGStab( msfem_matrix, 1e-8, 1e-8, 20000, VERBOSE );
     msfem_biCGStab( msfem_rhs, coarse_msfem_solution );

     if ( data_file_ )
      {
        if (data_file_->is_open())
         {
           *data_file_ << "---------------------------------------------------------------------------------" << std :: endl;
           *data_file_ << "MsFEM problem solved in " << assembleTimer.elapsed() << "s." << std :: endl << std :: endl << std :: endl;
         }
      }      
      
      
     // oneLinePrint( std::cout , solution );

     // copy coarse grid function (defined on the subgrid) into a fine grid function
     solution.clear();


#if 1
    int fine_level = grid.maxLevel();

    std :: cout << "Indentifying coarse scale part of the MsFEM solution... ";

    //! copy coarse scale part of MsFEM solution into a function defined on the fine grid
    // ------------------------------------------------------------------------------------
    typedef typename HostEntity :: template Codim< 2 > :: EntityPointer HostNodePointer;
    
    typedef typename GridPart :: IntersectionIteratorType HostIntersectionIterator;
    
    for( HostgridIterator it = discreteFunctionSpace_.begin(); it != discreteFunctionSpace_.end(); ++it )
      {

        typename HostEntity :: template Codim< 0 > :: EntityPointer coarse_father = it;
        for (int lev = 0; lev < ( fine_level - coarse_level) ; ++lev)
          coarse_father = coarse_father->father();
        if ( subGrid.template contains<0>( coarse_father ) == false )
          { std :: cout << "Error in msfem_solver.hh: Entity not in Subgrid!" << std :: endl; }

        CoarseGridEntityPointer coarse_it = subGrid.template getSubGridEntity<0>( *coarse_father );

        LinearLagrangeFunction2D< SubgridDiscreteFunctionSpace > interpolation_coarse( coarse_it );
        interpolation_coarse.set_corners( coarse_msfem_solution );

        LocalFunction host_loc_value = coarse_scale_part.localFunction( *it );

        int number_of_nodes = (*it).template count<2>();
        if (!( number_of_nodes == host_loc_value.baseFunctionSet().numBaseFunctions() ))
         { std :: cout << "Error! Inconsistency in 'msfem_solver.hh'." << std :: endl; }

        for ( int i = 0; i < number_of_nodes; i += 1 )
         {

            const HostNodePointer node = (*it).template subEntity<2>(i);

            DomainType coordinates_of_node = node->geometry().corner(0);
            if ( !( coordinates_of_node == it->geometry().corner(i)) )
             { std :: cout << "Error! Inconsistency in 'msfem_solver.hh'." << std :: endl; }

            RangeType coarse_value(0.0);
            interpolation_coarse.evaluate( coordinates_of_node, coarse_value );

            // int global_index_node = gridPart.indexSet().index( *node );
            host_loc_value[ i ] = coarse_value;

          }

      }
    std :: cout << " done." << std :: endl;
#endif
     // ------------------------------------------------------------------------------------

#if 1
     fine_scale_part.clear();
     

     int number_of_nodes = grid.size( grid.maxLevel(), 2 /*codim*/ );
     std :: vector< std :: vector < HostEntityPointer > > entities_sharing_same_node( number_of_nodes );
      
     for( HostgridIterator it = discreteFunctionSpace_.begin(); it != discreteFunctionSpace_.end(); ++it )
        {
	  int number_of_nodes_in_entity = (*it).template count<2>();
	  for ( int i = 0; i < number_of_nodes_in_entity; i += 1 )
	    {
	      const typename HostEntity :: template Codim< 2 > :: EntityPointer node = (*it).template subEntity<2>(i);
	      int global_index_node = gridPart.indexSet().index( *node );
	      
	      entities_sharing_same_node[ global_index_node ].push_back( it );
	    }
        }
        

    //! indentify fine scale part of MsFEM solution (including the projection!)
    // ------------------------------------------------------------------------------------

    std :: cout << "Indentifying fine scale part of the MsFEM solution... ";
    // iterator ueber coarse space
     for( CoarseGridIterator it = coarseDiscreteFunctionSpace.begin(); it != endit; ++it )
       {

         // the coarse entity 'T'	 
         HostEntityPointer coarse_host_entity_pointer = subGrid.template getHostEntity<0>( *it );
         const HostEntity& coarse_host_entity = *coarse_host_entity_pointer;

         DiscreteFunction correction_on_U_T( "correction_on_U_T", discreteFunctionSpace_ );
         correction_on_U_T.clear();

         const typename HostGrid :: Traits :: LevelIndexSet& hostGridLevelIndexSet
               = coarse_space.gridPart().grid().levelIndexSet( subGrid.maxLevel() );

         int index = hostGridLevelIndexSet.index( coarse_host_entity );

         // the sub grid U(T) that belongs to the coarse_grid_entity T
         SubGridType& sub_grid_U_T = subgrid_list.getSubGrid( index );
         SubGridPart subGridPart( sub_grid_U_T );

         SubgridDiscreteFunctionSpace localDiscreteFunctionSpace( subGridPart );

         SubgridDiscreteFunction local_problem_solution_e0( "Local problem Solution e_0", localDiscreteFunctionSpace );
         local_problem_solution_e0.clear();

         SubgridDiscreteFunction local_problem_solution_e1( "Local problem Solution e_1", localDiscreteFunctionSpace );
         local_problem_solution_e1.clear();

         // --------- load local solutions -------

         char location_lps[50];
         sprintf( location_lps, "/local_problems/_localProblemSolutions_%d", index );
         std::string location_lps_s( location_lps );

         std :: string local_solution_location;

         // the file/place, where we saved the solutions of the cell problems
         local_solution_location = path_ + location_lps_s;

         bool reader_is_open = false;
         // reader for the cell problem data file:
         DiscreteFunctionReader discrete_function_reader( (local_solution_location).c_str() );
         reader_is_open = discrete_function_reader.open();

         if (reader_is_open)
          { discrete_function_reader.read( 0, local_problem_solution_e0 ); }

         if (reader_is_open)
          { discrete_function_reader.read( 1, local_problem_solution_e1 ); }

         typedef typename SubgridDiscreteFunction :: LocalFunctionType SubgridLocalFunction;
	 SubgridLocalFunction local_coarse_part = coarse_msfem_solution.localFunction( *it );

         // 1 point quadrature!! We only need the gradient of the coarse scale part on the element, which is a constant.
         CachingQuadrature< SubGridPart, 0 > one_point_quadrature( *it, 0 );

         JacobianRangeType grad_coarse_msfem_on_entity;
         local_coarse_part.jacobian( one_point_quadrature[ 0 ], grad_coarse_msfem_on_entity );

         //!
         // std :: cout << "grad_coarse_msfem_on_entity[ 0 ][ 1 ] = " << grad_coarse_msfem_on_entity[ 0 ][ 1 ] << std :: endl;
         // std :: cout << "grad_coarse_msfem_on_entity[ 0 ][ 0 ] = " << grad_coarse_msfem_on_entity[ 0 ][ 0 ] << std :: endl;
         local_problem_solution_e0 *= grad_coarse_msfem_on_entity[ 0 ][ 0 ];
         local_problem_solution_e1 *= grad_coarse_msfem_on_entity[ 0 ][ 1 ];
         local_problem_solution_e0 += local_problem_solution_e1;

         // oneLinePrint( std::cout , local_problem_solution_e0 );

         subgrid_to_hostrid_projection( local_problem_solution_e0, correction_on_U_T );
         // hol die den Gradient und addiere.
#if 1
         if ( sub_grid_U_T.maxLevel() != discreteFunctionSpace_.gridPart().grid().maxLevel() )
           { std :: cout << "Error: MaxLevel of SubGrid not identical to MaxLevel of FineGrid." << std :: endl; }

         correction_on_U_T.clear();
 
	 typedef typename SubgridDiscreteFunctionSpace :: IteratorType SubgridIterator;
         typedef typename SubgridIterator :: Entity SubgridEntity;
         typedef typename SubgridDiscreteFunction :: LocalFunctionType SubgridLocalFunction;

         SubgridIterator sub_endit = localDiscreteFunctionSpace.end();
         for( SubgridIterator sub_it = localDiscreteFunctionSpace.begin(); sub_it != sub_endit; ++sub_it )
           {

             const SubgridEntity &sub_entity = *sub_it;

             HostEntityPointer fine_host_entity_pointer = sub_grid_U_T.template getHostEntity<0>( *sub_it );
             const HostEntity& fine_host_entity = *fine_host_entity_pointer;

             HostEntityPointer father = fine_host_entity_pointer;
             for (int lev = 0; lev < ( fine_level - coarse_level) ; ++lev)
                 father = father->father();

             bool entities_identical = true;
             int number_of_nodes = (*father).template count<2>();
             for ( int k = 0; k < number_of_nodes; k += 1 )
               {
		   if ( !(father->geometry().corner(k) == coarse_host_entity.geometry().corner(k)) )
                    { entities_identical = false; }
               }

             if ( entities_identical == false )
               {
                    // std :: cout << "father->geometry().corner(0) = " << father->geometry().corner(0) << std :: endl;
                    // std :: cout << "father->geometry().corner(1) = " << father->geometry().corner(1) << std :: endl;
                    // std :: cout << "father->geometry().corner(2) = " << father->geometry().corner(2) << std :: endl;
                    // std :: cout << "coarse_host_entity.geometry().corner(0) = " << coarse_host_entity.geometry().corner(0) << std :: endl;
                    // std :: cout << "coarse_host_entity.geometry().corner(1) = " << coarse_host_entity.geometry().corner(1) << std :: endl;
                    // std :: cout << "coarse_host_entity.geometry().corner(2) = " << coarse_host_entity.geometry().corner(2) << std :: endl << std :: endl;
                   continue;
               }

             SubgridLocalFunction sub_loc_value = local_problem_solution_e0.localFunction( sub_entity );
             LocalFunction host_loc_value = correction_on_U_T.localFunction( fine_host_entity );

	     int number_of_nodes_entity = (*sub_it).template count<2>();
             for ( int i = 0; i < number_of_nodes_entity; i += 1 )
               {
                 const typename HostEntity :: template Codim< 2 > :: EntityPointer node = fine_host_entity.template subEntity<2>(i);
	         int global_index_node = gridPart.indexSet().index( *node );

                 // vector of coarse entities that share the above node
                 std :: vector < HostEntityPointer > coarse_entities;

		 // count the number of different coarse-grid-entities that share the above node
		 for( int j = 0; j < entities_sharing_same_node[global_index_node].size(); j += 1 )
                   {
                      HostEntityPointer inner_it = entities_sharing_same_node[ global_index_node ][ j ];
                      for (int lev = 0; lev < ( fine_level - coarse_level) ; ++lev )
                         inner_it = inner_it->father();

                      bool new_entity_found = true;
                      for ( int k = 0; k < coarse_entities.size(); k += 1 )
                        {
                          if ( coarse_entities[k] == inner_it )
                           { new_entity_found = false; }
                        }
                      if ( new_entity_found == true )
                       { coarse_entities.push_back( inner_it ); }

                   }

                 host_loc_value[ i ] = ( sub_loc_value[ i ] / coarse_entities.size() );
               }

           }

#endif

        fine_scale_part+=correction_on_U_T;


      }
    std :: cout << " done." << std :: endl;
#endif
    // ------------------------------------------------------------------------------------

     // Auf Grobskalen MsFEM Anteil noch Feinksalen MsFEM Anteil aufaddieren.
     solution += coarse_scale_part;
     solution += fine_scale_part;

   }


   //! the following methods are not yet implemented, however note that the required tools are 
   //! already available via 'righthandside_assembler.hh' and 'elliptic_fem_matrix_assembler.hh'!

   template< class DiffusionOperatorType, class ReactionTermType, class SourceTermType >
   void solve( )
   {
     std :: cout << "No implemented!" << std :: endl;
   }


   template< class DiffusionOperatorType, class ReactionTermType, class SourceTermType, class SecondSourceTermType>
   void solve( )
   {
     std :: cout << "No implemented!" << std :: endl;
   }



  };



}

#endif // #ifndef Elliptic_MSEM_Solver_HH
