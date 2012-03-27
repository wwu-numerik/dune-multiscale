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
   
   typedef typename HostGrid ::template Codim< 0 > :: template Partition< All_Partition > :: LevelIterator LevelEntityIteratorType;

   typedef typename DiscreteFunctionSpace :: IteratorType HostgridIterator;
   
   typedef typename HostgridIterator :: Entity HostEntity;
   
   typedef typename HostEntity :: EntityPointer HostEntityPointer; 

   
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


   // - ∇ (A(x,∇u)) + b ∇u + c u = f - divG
   // then:
   // A --> diffusion operator ('DiffusionOperatorType')
   // b --> advective part ('AdvectionTermType')
   // c --> reaction part ('ReactionTermType')
   // f --> 'first' source term, scalar ('SourceTermType')
   // G --> 'second' source term, vector valued ('SecondSourceTermType')

//! wieder einbinden:
#if 0
   // homogenous Dirchilet boundary condition!:
   template< class DiffusionOperator, class SourceTerm >
   void solve_dirichlet_zero( const DiffusionOperator &diffusion_op,
                              const SourceTerm &f,
                              const int coarse_level,
                                    DiscreteFunction &solution )
   {
     HostGrid &grid = discreteFunctionSpace_.gridPart().grid();
     const GridPart &gridPart = discreteFunctionSpace_.gridPart();
     int number_of_level_host_entities = grid.size( coarse_level, 0 /*codim*/ );
     
     //! default - no layers
     // number of layers per coarse grid entity T:  U(T) is created by enrichting T with n(T)-layers.
     std :: vector < int > number_of_layers( number_of_level_host_entities );
     for ( int i = 0; i < number_of_level_host_entities; i+=1 )
        { number_of_layers[i] = 0; }
     
     solve_dirichlet_zero( diffusion_op, f, coarse_level, number_of_layers, solution );
   }
#endif
   
   // homogenous Dirchilet boundary condition!:
   template< class DiffusionOperator, class SourceTerm >
   void solve_dirichlet_zero( const DiffusionOperator &diffusion_op,
                              const SourceTerm &f,
                              const int coarse_level,
                              // number of layers per coarse grid entity T:  U(T) is created by enrichting T with n(T)-layers.
                              std :: vector < int >& number_of_layers,
DiscreteFunctionSpace& discreteFunctionSpace2, //! loeschen!!!!!!!!!!!
                              DiscreteFunction &solution )
   {

     // discrete elliptic MsFEM operator (corresponds with MsFEM Matrix)
     typedef DiscreteEllipticMsFEMOperator< SubgridDiscreteFunction,
                                            DiscreteFunction,
                                            DiffusionOperator > EllipticMsFEMOperatorType;

     HostGrid &grid = discreteFunctionSpace_.gridPart().grid();
     const GridPart &gridPart = discreteFunctionSpace_.gridPart();

     LevelEntityIteratorType coarse_level_it = grid.template lbegin< 0 >( coarse_level );


#if 1
     HostGrid &grid2 = discreteFunctionSpace2.gridPart().grid();
     const GridPart &gridPart2 = discreteFunctionSpace2.gridPart();
#endif


     // create subgrid:
     SubGridType subGrid( grid2 /*!*!/*/ );
     subGrid.createBegin();

     //subGrid.insertLevel(coarse_level);

     for( ; coarse_level_it != grid2.template lend< 0 >( coarse_level ); ++coarse_level_it ) /*!/!*/
         subGrid.insertPartial( *coarse_level_it );

     subGrid.createEnd();

     subGrid.report();

     SubGridPart subGridPart( subGrid );

     SubgridDiscreteFunctionSpace coarseDiscreteFunctionSpace( subGridPart );

     SubgridDiscreteFunction coarse_msfem_solution( "Coarse Part MsFEM Solution", coarseDiscreteFunctionSpace );
     coarse_msfem_solution.clear();

     //! create subgrids:
     bool silence = false;

     typedef SubGridList< DiscreteFunction, SubGridType > SubGridListType;
     SubGridListType subgrid_list( discreteFunctionSpace_ , number_of_layers, coarse_level , silence );


     //! define the right hand side assembler tool
     // (for linear and non-linear elliptic and parabolic problems, for sources f and/or G )
     RightHandSideAssembler< SubgridDiscreteFunction > rhsassembler;

     //! define the discrete (elliptic) operator that describes our problem
     // ( effect of the discretized differential operator on a certain discrete function )
     EllipticMsFEMOperatorType elliptic_msfem_op( coarseDiscreteFunctionSpace,
                                                  discreteFunctionSpace_, discreteFunctionSpace2/*!loeschen*/,
                                                  subgrid_list,
                                                  diffusion_op, *data_file_, path_ );
     // discrete elliptic operator (corresponds with FEM Matrix)

//!!! loeschen:
#if 1
     DiscreteEllipticOperator< SubgridDiscreteFunction, DiffusionOperator, DummyMass< SubgridDiscreteFunctionSpace > > discrete_elliptic_op( coarseDiscreteFunctionSpace, diffusion_op );
#endif

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
//!     elliptic_msfem_op.assemble_matrix( msfem_matrix ); // einbinden!
discrete_elliptic_op.assemble_matrix( msfem_matrix, false ); //!loeschen///!!!
//! loeschen!!!!!!!!!!!!!!!
#if 1
    // boundary treatment

    for( CoarseGridIterator coarse_grid_it = coarseDiscreteFunctionSpace.begin(); coarse_grid_it != coarseDiscreteFunctionSpace.end(); ++coarse_grid_it )
    {

      const CoarseGridEntity &coarse_grid_entity = *coarse_grid_it;
      HostEntityPointer fine_entity_pointer = coarseDiscreteFunctionSpace.grid().template getHostEntity<0>( coarse_grid_entity );

      const HostEntity& fine_entity = *fine_entity_pointer;

      const CoarseGridLagrangePointSet &lagrangePointSet = coarseDiscreteFunctionSpace.lagrangePointSet( coarse_grid_entity );

      const IntersectionIterator iend = discreteFunctionSpace2.gridPart().iend( fine_entity );
      for( IntersectionIterator iit = discreteFunctionSpace2.gridPart().ibegin( fine_entity ); iit != iend; ++iit ) /*! 2 entfernen loeschen! */
        {
	  
           if ( iit->neighbor() ) //if there is a neighbor entity
            {
              // check if the neighbor entity is in the subgrid
              const HostEntityPointer neighborFineEntityPointer = iit->outside();
              const HostEntity& neighborFineEntity = *neighborFineEntityPointer;
              if ( coarseDiscreteFunctionSpace.grid().template contains<0>( neighborFineEntity ) )
               {
                 continue;
               }

            }

           const int face = (*iit).indexInInside();
           const CoarseGridFaceDofIterator fdend = lagrangePointSet.template endSubEntity< 1 >( face );
           for( CoarseGridFaceDofIterator fdit = lagrangePointSet.template beginSubEntity< 1 >( face ); fdit != fdend; ++fdit )
              msfem_matrix.localMatrix( coarse_grid_entity, coarse_grid_entity ).unitRow( *fdit );	
        }
  
    }
#endif


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
     
#if 0
     //subGrid.globalRefine( 2 );
       
     LevelEntityIteratorType coarse_level_it_2 = grid.template lbegin< 0 >( grid.maxLevel() );
     SubGridType subGrid_2( grid ); 
     subGrid_2.createBegin();
     

     for( ; coarse_level_it_2 != grid.template lend< 0 >( grid.maxLevel() ); ++coarse_level_it_2 )
         subGrid_2.insertPartial( *coarse_level_it_2 );

     subGrid_2.createEnd();
     SubGridPart subGridPart_2( subGrid_2 );
     SubgridDiscreteFunctionSpace coarseDiscreteFunctionSpace_2( subGridPart_2 );
     SubgridDiscreteFunction coarse_msfem_solution_2( "Coarse Part MsFEM Solution", coarseDiscreteFunctionSpace_2 );
     coarse_msfem_solution_2.clear();
 
#if 0
//! type of restrict-prolong operator
typedef RestrictProlongDefault< SubgridDiscreteFunction >
  RestrictProlongOperatorType;
//! type of the adaption manager
typedef AdaptationManager< SubGridType, RestrictProlongOperatorType >
  AdaptationManagerType;
  
// one for the discreteFunctionSpace
RestrictProlongOperatorType rp( coarse_msfem_ );
AdaptationManagerType adaptationManager( grid, rp );

for( CoarseGridIterator it = coarseDiscreteFunctionSpace.begin(); it != endit; ++it )
       { subGrid.mark( 2 , *it ); }

adaptationManager.adapt();
#endif
#if 0
// type of restrict-prolong operator
typedef RestrictProlongDefault< DiscreteFunction >
  RestrictProlongOperatorType;
// type of the adaption manager
typedef AdaptationManager< HostGrid, RestrictProlongOperatorType >
  AdaptationManagerType;
  
// one for the discreteFunctionSpace
RestrictProlongOperatorType rp( solution );
AdaptationManagerType adaptationManager( grid, rp );

for( HostgridIterator it = discreteFunctionSpace_.begin(); it != discreteFunctionSpace_.end(); ++it )
       { grid.mark( -2 , *it ); }

adaptationManager.adapt();
#endif

//!eIt->geometry().corner(0)[0]

#if 0
for( CoarseGridIterator it = coarseDiscreteFunctionSpace.begin(); it != endit; ++it )
       { subGrid.mark( 2 , *it ); }
       
        subGrid.preAdapt();
        subGrid.adapt();
        subGrid.postAdapt();
#endif

#if 0
typedef typename SubgridDiscreteFunction::DofIteratorType SubgridDofIteratorType;
typedef typename DiscreteFunction::DofIteratorType DofIteratorType;

SubgridDofIteratorType sub_it = coarse_msfem_solution.dbegin();
DofIteratorType it = solution.dbegin();
#endif
#if 0
int number_dofs = 0;
for ( ; it != solution.dend(); ++it )
 number_dofs += 1;
std :: cout << " number_dofs = " << number_dofs << std :: endl;

int number_dofs_sub = 0;
for ( ; sub_it != coarse_msfem_solution.dend(); ++sub_it )
 number_dofs_sub += 1;
std :: cout << " number_dofs_sub = " << number_dofs_sub << std :: endl;
#endif
#if 0
for ( ; it != solution.dend(); ++it )
 {
    *it = *sub_it;
    ++sub_it;
 }
#endif  
#if 0
for( HostgridIterator it = discreteFunctionSpace_.begin(); it != discreteFunctionSpace_.end(); ++it )
       { grid.mark( 2 , *it ); }

adaptationManager.adapt();
#endif
#if 0
     for( CoarseGridIterator it = coarseDiscreteFunctionSpace.begin(); it != endit; ++it )
       {

         //HostEntityPointer host_entity_pointer = subGrid.template getHostEntity<0>( *it );
         //const HostEntity& host_entity = *host_entity_pointer;


         CoarseGridLocalFunction sub_loc_value = coarse_msfem_solution.localFunction( *it );
         CoarseGridLocalFunction host_loc_value = coarse_msfem_solution_2.localFunction( *it );

         const unsigned int numBaseFunctions = sub_loc_value.baseFunctionSet().numBaseFunctions();
         for( unsigned int i = 0; i < numBaseFunctions; ++i )
           {
             host_loc_value[ i ] = sub_loc_value[ i ];
           }

       }
#endif
#if 0
     for( CoarseGridIterator it = coarseDiscreteFunctionSpace_2.begin(); it != endit; ++it )
       {

         HostEntityPointer host_entity_pointer = subGrid_2.template getHostEntity<0>( *it );
         const HostEntity& host_entity = *host_entity_pointer;


         CoarseGridLocalFunction sub_loc_value = coarse_msfem_solution_2.localFunction( *it );
         LocalFunction host_loc_value = solution.localFunction( host_entity );

         const unsigned int numBaseFunctions = sub_loc_value.baseFunctionSet().numBaseFunctions();
         for( unsigned int i = 0; i < numBaseFunctions; ++i )
           {
             host_loc_value[ i ] = sub_loc_value[ i ];
           }

       }
#endif
#endif
#if 0  
     

#if 1
     for( FineGridIterator it = fineDiscreteFunctionSpace.begin(); it != endit; ++it )
       {
	 
      const Entity &entity = *it;
      const Geometry &geometry = entity.geometry();



      // for constant diffusion "2*discreteFunctionSpace_.order()" is sufficient, for the general case, it is better to use a higher order quadrature:
      Quadrature quadrature( entity, 2*discreteFunctionSpace_.order()+2 );
      const size_t numQuadraturePoints = quadrature.nop();
      for( size_t quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint )
      {
        // local (barycentric) coordinates (with respect to entity)
        const typename Quadrature::CoordinateType &local_point = quadrature.point( quadraturePoint );

        DomainType global_point = geometry.global( local_point );


         HostEntityPointer host_entity_pointer = subGrid.template getHostEntity<0>( *it );
         const HostEntity& host_entity = *host_entity_pointer;


         CoarseGridLocalFunction sub_loc_value = coarse_msfem_solution.localFunction( *it );
         LocalFunction host_loc_value = solution.localFunction( host_entity );

         const unsigned int numBaseFunctions = sub_loc_value.baseFunctionSet().numBaseFunctions();
         for( unsigned int i = 0; i < numBaseFunctions; ++i )
           {
             host_loc_value[ i ] = sub_loc_value[ i ];
           }
       }
#endif
//DomainType x(0.0);
//RangeType y(0.0);
//solution.evaluate(x,y);
#endif     
#if 0
std:: cout << "coarseDiscreteFunctionSpace.size() = " << coarseDiscreteFunctionSpace.size() << std :: endl;
std:: cout << "discreteFunctionSpace_.size() = " << discreteFunctionSpace_.size() << std :: endl;

     for( CoarseGridIterator it = coarseDiscreteFunctionSpace.begin(); it != endit; ++it )
       {

         HostEntityPointer host_entity_pointer = subGrid.template getHostEntity<0>( *it );
         const HostEntity& host_entity = *host_entity_pointer;


         CoarseGridLocalFunction sub_loc_value = coarse_msfem_solution.localFunction( *it );
         LocalFunction host_loc_value = solution.localFunction( host_entity );

         const unsigned int numBaseFunctions = sub_loc_value.baseFunctionSet().numBaseFunctions();
         for( unsigned int i = 0; i < numBaseFunctions; ++i )
           {
             RangeType coarse_value(0.0);
             DomainType point = it->geometry().corner(i);
             /*sub_loc_value*/coarse_msfem_solution.evaluate(point, coarse_value);
             std :: cout << "sub_loc_value[ i ] = " << sub_loc_value[ i ] << " und coarse_value = " << coarse_value << std :: endl;
             //host_loc_value[ i ] = sub_loc_value[ i ];
           }

       }
#endif




#if 1
typedef typename SubgridDiscreteFunction::DofIteratorType SubgridDofIteratorType;

SubgridDofIteratorType sub_it = coarse_msfem_solution.dbegin();
for ( ; sub_it != coarse_msfem_solution.dend(); ++sub_it )
 {
    std :: cout << "*sub_it = " << *sub_it << std :: endl;
 }
#endif
std :: cout << " " << std :: endl;
//abort();

#if 1

    typedef typename HostEntity :: template Codim< 2 > :: EntityPointer HostNodePointer;
    
    typedef typename GridPart :: IntersectionIteratorType HostIntersectionIterator;
    
    for( HostgridIterator it = discreteFunctionSpace_.begin(); it != discreteFunctionSpace_.end(); ++it )
      {

        int coarse_level = subGrid.maxLevel();
        int fine_level = grid.maxLevel();

        typename HostEntity :: template Codim< 0 > :: EntityPointer coarse_father = it;
        for (int lev = 0; lev < ( fine_level - coarse_level) ; ++lev)
          coarse_father = coarse_father->father();
        if ( subGrid.template contains<0>( coarse_father ) == false )
          { std :: cout << "Error in msfem_solver.hh: Entity not in Subgrid!" << std :: endl; }

        CoarseGridEntityPointer coarse_it = subGrid.template getSubGridEntity<0>( *coarse_father );

//        LinearLagrangeFunction2D< SubgridDiscreteFunctionSpace > interpolation_coarse( coarse_it );
//        interpolation_coarse.set_corners( coarse_msfem_solution );

#if 0
        CoarseGridLocalFunction sub_loc_value = coarse_msfem_solution.localFunction( *coarse_it );
        RangeType value(1.0);
        LinearLagrangeFunction2D< SubgridDiscreteFunctionSpace > interpolation_coarse
          ( coarse_it->geometry().corner(0), sub_loc_value[0],
            coarse_it->geometry().corner(1), sub_loc_value[1],
            coarse_it->geometry().corner(2), sub_loc_value[2] );
#endif

        LocalFunction host_loc_value = solution.localFunction( *it );

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
            //interpolation_coarse.evaluate( coordinates_of_node, coarse_value );

coarse_msfem_solution.evaluate( coordinates_of_node, coarse_value );
#if 0
if ( host_loc_value[ i ] != 0.0 )
{ std :: cout << "host_loc_value[ i ] = " << host_loc_value[ i ] << " und coarse_value = " << coarse_value << std :: endl; }
#endif

            // int global_index_node = gridPart.indexSet().index( *node );
            host_loc_value[ i ] = coarse_value;

          }

      }
#endif

     std :: cout << "Auf Grobskalen MsFEM Anteil noch Feinksalen MsFEM Anteil aufaddieren." << std :: endl << std :: endl;  
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
