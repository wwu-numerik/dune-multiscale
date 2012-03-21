#ifndef Elliptic_MSEM_Solver_HH
#define Elliptic_MSEM_Solver_HH

#include <dune/common/fmatrix.hh>

#include <dune/subgrid/subgrid.hh>


#include <dune/fem/solver/oemsolver/oemsolver.hh>
#include <dune/fem/solver/inverseoperators.hh>
#include <dune/fem/operator/discreteoperatorimp.hh>

#include <dune/fem/operator/2order/lagrangematrixsetup.hh>
#include <dune/fem/operator/matrix/spmatrix.hh>

#include <dune/multiscale/tools/assembler/righthandside_assembler.hh>
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

   typedef typename DiscreteFunctionSpace :: GridType Grid;

   typedef typename GridType ::template Codim< 0 > :: template Partition< All_Partition > :: LevelIterator LevelEntityIteratorType;

   enum { faceCodim = 1 };

   typedef typename GridPart :: IntersectionIteratorType IntersectionIterator;

   typedef typename LagrangePointSet :: template Codim< faceCodim > 
                                     :: SubEntityIteratorType
    FaceDofIterator;

   // --------------------------- subgrid typedefs ------------------------------------

   typedef SubGrid< GridType::dimension , Grid > SubGridType; 

   typedef LeafGridPart< SubGridType > SubGridPart;

   typedef LagrangeDiscreteFunctionSpace < FunctionSpace, SubGridPart, 1 > //1=POLORDER
      SubgridDiscreteFunctionSpace;

   typedef AdaptiveDiscreteFunction < SubgridDiscreteFunctionSpace > SubgridDiscreteFunction;

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

  public:
   Elliptic_MsFEM_Solver( const DiscreteFunctionSpace &discreteFunctionSpace )
     : discreteFunctionSpace_( discreteFunctionSpace ),
       data_file_( NULL )
     {}

   Elliptic_MsFEM_Solver( const DiscreteFunctionSpace &discreteFunctionSpace, std :: ofstream& data_file )
     : discreteFunctionSpace_( discreteFunctionSpace ),
       data_file_( &data_file )
     {}

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

     // discrete elliptic MsFEM operator (corresponds with MsFEM Matrix)
     typedef DiscreteEllipticMsFEMOperator< DiscreteFunction, DiffusionOperator > EllipticMsFEMOperatorType;
     const Grid &grid = discreteFunctionSpace_.grid();

     LevelEntityIteratorType coarse_level_it = grid.template lbegin< 0 >( coarse_level );

     // create subgrid:
     SubGridType subGrid( grid );
     subGrid.createBegin();



//! 1. Anlegen des Coarse-Level discrete finction spaces mittels subgrid


#if 0
     const GridPart &gridPart = discreteFunctionSpace_.gridPart();

     //! define the right hand side assembler tool
     // (for linear and non-linear elliptic and parabolic problems, for sources f and/or G )
     RightHandSideAssembler< DiscreteFunctionType > rhsassembler;

     //! define the discrete (elliptic) operator that describes our problem
     // ( effect of the discretized differential operator on a certain discrete function )
     DiscreteEllipticOperator< DiscreteFunction, DiffusionOperator, DummyMassType > discrete_elliptic_op( discreteFunctionSpace_, diffusion_op );
     // discrete elliptic operator (corresponds with FEM Matrix)

     //! (stiffness) matrix
     FEMMatrix fem_matrix( "FEM stiffness matrix", discreteFunctionSpace_, discreteFunctionSpace_ );

     //! right hand side vector
     // right hand side for the finite element method:
     DiscreteFunction fem_rhs( "fem newton rhs", discreteFunctionSpace_ );
     fem_rhs.clear();

     std :: cout << "Solving linear problem." << std :: endl;

     if ( data_file_ )
      {
        if (data_file_->is_open())
         {
           *data_file_ << "Solving linear problem with standard FEM and resolution level " << discreteFunctionSpace_.grid().maxLevel() << "." << std :: endl;
           *data_file_ << "------------------------------------------------------------------------------" << std :: endl;
         }
      }

     // to assemble the computational time
     Dune::Timer assembleTimer;

     // assemble the stiffness matrix
     discrete_elliptic_op.assemble_matrix( fem_matrix );

     std::cout << "Time to assemble standard FEM stiffness matrix: " << assembleTimer.elapsed() << "s" << std::endl;

     if ( data_file_ )
      {
        if (data_file_->is_open())
         {
           *data_file_ << "Time to assemble standard FEM stiffness matrix: " << assembleTimer.elapsed() << "s" << std::endl;
         }
      }

     // assemble right hand side
     rhsassembler.template assemble< 2 * DiscreteFunctionSpace :: polynomialOrder + 2 >( f , fem_rhs);

     //oneLinePrint( std::cout , fem_rhs );


     // --- boundary treatment ---
     // set the dirichlet points to zero (in righ hand side of the fem problem)
     typedef typename DiscreteFunctionSpace :: IteratorType EntityIterator;
     EntityIterator endit = discreteFunctionSpace_.end();
     for( EntityIterator it = discreteFunctionSpace_.begin(); it != endit; ++it )
       {

          IntersectionIterator iit = gridPart.ibegin( *it );
          const IntersectionIterator endiit = gridPart.iend( *it );
          for( ; iit != endiit; ++iit )
            {

              if( !(*iit).boundary() )
                continue;

              LocalFunction rhsLocal = fem_rhs.localFunction( *it );
              const LagrangePointSet &lagrangePointSet
                = discreteFunctionSpace_.lagrangePointSet( *it );

              const int face = (*iit).indexInInside();

              FaceDofIterator faceIterator
                = lagrangePointSet.template beginSubEntity< faceCodim >( face );
              const FaceDofIterator faceEndIterator
                = lagrangePointSet.template endSubEntity< faceCodim >( face );
              for( ; faceIterator != faceEndIterator; ++faceIterator )
                rhsLocal[ *faceIterator ] = 0;

            }

       }
     // --- end boundary treatment ---

     InverseFEMMatrix fem_biCGStab( fem_matrix, 1e-8, 1e-8, 20000, VERBOSE );
     fem_biCGStab( fem_rhs, solution );

     if ( data_file_ )
      {
        if (data_file_->is_open())
         {
           *data_file_ << "---------------------------------------------------------------------------------" << std :: endl;
           *data_file_ << "Standard FEM problem solved in " << assembleTimer.elapsed() << "s." << std :: endl << std :: endl << std :: endl;
         }
      }

     // oneLinePrint( std::cout , solution );
#endif

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
