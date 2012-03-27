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
   
   
   // homogenous Dirchilet boundary condition!:
   template< class DiffusionOperator, class SourceTerm >
   void solve_dirichlet_zero( const DiffusionOperator &diffusion_op,
                              const SourceTerm &f,
                              const int coarse_level,
			       // number of layers per coarse grid entity T:  U(T) is created by enrichting T with n(T)-layers.
			       std :: vector < int >& number_of_layers,
                              DiscreteFunction &solution )
   {

     // discrete elliptic MsFEM operator (corresponds with MsFEM Matrix)
     typedef DiscreteEllipticMsFEMOperator< SubgridDiscreteFunction,
                                            DiscreteFunction,
					     DiffusionOperator > EllipticMsFEMOperatorType;

     HostGrid &grid = discreteFunctionSpace_.gridPart().grid();
     const GridPart &gridPart = discreteFunctionSpace_.gridPart();

     LevelEntityIteratorType coarse_level_it = grid.template lbegin< 0 >( coarse_level );

     // create subgrid:
     SubGridType subGrid( grid );
     subGrid.createBegin();

     for( ; coarse_level_it != grid.template lend< 0 >( coarse_level ); ++coarse_level_it )
         subGrid.insert/*!Partial*/( *coarse_level_it );

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
                                                  discreteFunctionSpace_, subgrid_list,
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
     elliptic_msfem_op.assemble_matrix( msfem_matrix ); 
     
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

#if 1
for( CoarseGridIterator it = coarseDiscreteFunctionSpace.begin(); it != endit; ++it )
       { subGrid.mark( 2 , *it ); }
       
        subGrid.preAdapt();
        subGrid.adapt();
        subGrid.postAdapt();
#endif

typedef typename SubgridDiscreteFunction::DofIteratorType SubgridDofIteratorType;
typedef typename DiscreteFunction::DofIteratorType DofIteratorType;

SubgridDofIteratorType sub_it = coarse_msfem_solution.dbegin();
DofIteratorType it = solution.dbegin();
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
             host_loc_value[ i ] = sub_loc_value[ i ];
           }

       }
#endif



#if 1

    typedef typename HostEntity :: template Codim< 2 > :: EntityPointer HostNodePointer;
    
    typedef typename GridPart :: IntersectionIteratorType HostIntersectionIterator;
    
    for( HostgridIterator it = discreteFunctionSpace_.begin(); it != discreteFunctionSpace_.end(); ++it )
      {
	
        LocalFunction host_loc_value = solution.localFunction( *it );
	
	int number_of_nodes = (*it).template count<2>();
	if (!( number_of_nodes == host_loc_value.baseFunctionSet().numBaseFunctions() ))
	 { std :: cout << "Error!" << std :: endl; }
	  
        for ( int i = 0; i < number_of_nodes; i += 1 )
	  {

	    // Ecken mit Knoten identifizieren und Indexset nutzen!
	    const HostNodePointer node = (*it).template subEntity<2>(i);
	    int global_index_node = gridPart.indexSet().index( *node );
	    std :: cout << "node = " << node->geometry().corner(0) << std :: endl;
	    std :: cout << "node->geometry().corner(" << i << ") = " << it->geometry().corner(i) << std :: endl;
            std :: cout << "host_loc_value[ i ] = " << host_loc_value[ i ] << std :: endl << std :: endl;
	    if ( !(node->geometry().corner(0) == it->geometry().corner(i)) )
	    { std :: cout << "Error!" << std :: endl; }


	  }
	  
      }
#if 0
for ( ; it != solution.dend(); ++it )
 {
    std :: cout << "it.index() = " << it->index() << std :: endl;
 }
#endif  
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
