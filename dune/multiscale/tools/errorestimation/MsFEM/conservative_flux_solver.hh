#ifndef DiscreteEllipticMsFEMConservativeFluxSolver_HH
#define DiscreteEllipticMsFEMConservativeFluxSolver_HH

#include <vector>

#include <dune/fem/operator/matrix/spmatrix.hh>
#include <dune/common/fmatrix.hh>

#include <dune/fem/quadrature/cachingquadrature.hh>

#include <dune/grid/common/quadraturerules.hh>

#include <dune/fem/operator/common/operator.hh>


// FLUX_SOLVER_VERBOSE: 0 = false, 1 = true
#define FLUX_SOLVER_VERBOSE false

// VTK output for conservative flux solution
#define VTK_OUTPUT

// dune-subgrid include:
#include <dune/subgrid/subgrid.hh>

// dune-fem includes:
#include <dune/fem/gridpart/gridpart.hh>
#include <dune/fem/operator/2order/lagrangematrixsetup.hh>

namespace Dune
{

  // define output traits
  struct ConFluxProblemDataOutputParameters : public DataOutputParameters {

  public:

  std::string my_prefix_;
  std::string my_path_;

  void set_prefix( std::string my_prefix )
    {
      my_prefix_ = my_prefix;
      // std :: cout << "Set prefix. my_prefix_ = " << my_prefix_ << std :: endl;
    }

  void set_path( std::string my_path )
    {
      my_path_ = my_path;
    }

  // base of file name for data file
  std::string prefix() const 
    {
      if (my_prefix_ == "")
        return "solutions";
      else
        return my_prefix_;
    }

  // path where the data is stored
  std::string path() const 
    {
      if (my_path_ == "")
        return "data_output_msfem_conservative_flux";
      else
        return my_path_;

    }


  // format of output:
  int outputformat() const
    {
      //return 0; // GRAPE (lossless format)
      return 1; // VTK
      //return 2; // VTK vertex data
      //return 3; // gnuplot
    }


  };


  // Imp stands for Implementation
  template< class SubGridDiscreteFunctionType, class DiscreteFunctionType, class DiffusionOperatorType >
  class ConservativeFluxOperator
  : public Operator< typename SubGridDiscreteFunctionType::RangeFieldType,
                     typename SubGridDiscreteFunctionType::RangeFieldType,
		     SubGridDiscreteFunctionType,
		     SubGridDiscreteFunctionType >
  {
    
    typedef ConservativeFluxOperator< SubGridDiscreteFunctionType, DiscreteFunctionType, DiffusionOperatorType > This;

  public:
    typedef SubGridDiscreteFunctionType SubGridDiscreteFunction;
    typedef DiscreteFunctionType DiscreteFunction;
    typedef DiffusionOperatorType DiffusionModel;

    typedef typename SubGridDiscreteFunction::DiscreteFunctionSpaceType SubGridDiscreteFunctionSpace;
    typedef typename DiscreteFunction::DiscreteFunctionSpaceType DiscreteFunctionSpace;
    
    typedef typename DiscreteFunctionSpace::GridPartType GridPart;
    typedef typename DiscreteFunctionSpace::GridType GridType;

    typedef typename SubGridDiscreteFunctionSpace::GridPartType SubGridPart;
    typedef typename SubGridDiscreteFunctionSpace::GridType SubGridType;
    
    typedef typename DiscreteFunctionSpace::RangeFieldType RangeFieldType;
    typedef typename DiscreteFunctionSpace::DomainType DomainType;
    typedef typename DiscreteFunctionSpace::RangeType RangeType;
    typedef typename DiscreteFunctionSpace::JacobianRangeType
      JacobianRangeType;
      
  protected:
    static const int dimension = GridPart::GridType::dimension;
    static const int polynomialOrder = DiscreteFunctionSpace::polynomialOrder;

    typedef typename DiscreteFunction::LocalFunctionType LocalFunction;
    typedef typename SubGridDiscreteFunction::LocalFunctionType SubGridLocalFunction;
    
    typedef typename DiscreteFunctionSpace::BaseFunctionSetType BaseFunctionSet;
    typedef typename DiscreteFunctionSpace::LagrangePointSetType LagrangePointSet;
    typedef typename LagrangePointSet::template Codim< 1 >::SubEntityIteratorType FaceDofIterator;
    
    typedef typename DiscreteFunctionSpace::IteratorType Iterator;
    typedef typename Iterator::Entity Entity;
    typedef typename Entity :: EntityPointer EntityPointer;
    typedef typename Entity::Geometry Geometry;
    
    typedef typename GridPart::IntersectionIteratorType IntersectionIterator;
    typedef typename IntersectionIterator::Intersection Intersection;
    typedef typename Intersection::LocalCoordinate LocalCoordinate;


    typedef typename SubGridDiscreteFunctionSpace::BaseFunctionSetType SubGridBaseFunctionSet;
    typedef typename SubGridDiscreteFunctionSpace::LagrangePointSetType SubGridLagrangePointSet;
    typedef typename SubGridLagrangePointSet::template Codim< 1 >::SubEntityIteratorType SubGridFaceDofIterator;
    
    typedef typename SubGridDiscreteFunctionSpace::IteratorType SubGridIterator;
    typedef typename SubGridIterator::Entity SubGridEntity;
    typedef typename SubGridEntity::Geometry SubGridGeometry;
    
    typedef typename SubGridPart::IntersectionIteratorType SubGridIntersectionIterator;
    typedef typename SubGridIntersectionIterator::Intersection SubGridIntersection;

    typedef CachingQuadrature< GridPart, 0 > Quadrature;

    typedef QuadratureRule< double, 1 > FaceQuadratureRule;
    typedef CachingQuadrature< GridPart, 1 > FaceQuadrature;

    typedef CachingQuadrature< SubGridPart, 0 > SubGridQuadrature;

    typedef typename GridType :: template Codim< 1 > :: Geometry FaceGeometryType;

  public:

    ConservativeFluxOperator ( const SubGridDiscreteFunctionSpace &subDiscreteFunctionSpace,
                               const DiscreteFunctionSpace &discreteFunctionSpace,
                               const DiffusionModel &diffusion_op )
    : subDiscreteFunctionSpace_( subDiscreteFunctionSpace ),
      discreteFunctionSpace_( discreteFunctionSpace ),
      diffusion_operator_( diffusion_op )
    {}

  private:
    ConservativeFluxOperator ( const This & );

  public:

    // dummy operator
    virtual void
    operator() ( const SubGridDiscreteFunction &u, SubGridDiscreteFunction &w ) const;

    // assemble stiffness matrix for local problems
    template< class MatrixType >
    void assemble_matrix ( MatrixType &global_matrix ) const;

    // the right hand side assembler methods
    void assemble_RHS ( // direction 'e_i'
                        JacobianRangeType &e_i,
                        // solution of the local corrector problem
                        const SubGridDiscreteFunction &local_corrector_e_i,
                        // rhs local msfem problem:
                        SubGridDiscreteFunction &rhs_flux_problem ) const;

    void printLocalRHS( SubGridDiscreteFunction &rhs) const;

    double normRHS( SubGridDiscreteFunction &rhs) const;

  private:
    const SubGridDiscreteFunctionSpace& subDiscreteFunctionSpace_;
    const DiscreteFunctionSpace& discreteFunctionSpace_;
    const DiffusionModel& diffusion_operator_;

  };

  // dummy implementation of "operator()"
  // 'w' = effect of the discrete operator on 'u'
  template< class SubGridDiscreteFunctionImp, class DiscreteFunctionImp, class DiffusionImp >
  void ConservativeFluxOperator< SubGridDiscreteFunctionImp, DiscreteFunctionImp, DiffusionImp >
  :: operator() ( const SubGridDiscreteFunctionImp &u, SubGridDiscreteFunctionImp &w ) const 
  {

    std :: cout << "the ()-operator of the LocalProblemOperator class is not yet implemented and still a dummy." << std :: endl;
    std :: abort();

  }

#if 1
  //! assemble system matrix
  template< class SubGridDiscreteFunctionImp, class DiscreteFunctionImp, class DiffusionImp >
  template< class MatrixType >
  void ConservativeFluxOperator< SubGridDiscreteFunctionImp, DiscreteFunctionImp, DiffusionImp >::assemble_matrix ( MatrixType &global_matrix ) const
  {

    typedef typename MatrixType::LocalMatrixType LocalMatrix;

    global_matrix.reserve();
    global_matrix.clear();

    // local grid basis functions:
    std::vector< RangeType > phi( subDiscreteFunctionSpace_.mapper().maxNumDofs() );
    
    const SubGridIterator end = subDiscreteFunctionSpace_.end();
    for( SubGridIterator it = subDiscreteFunctionSpace_.begin(); it != end; ++it )
    {

      const SubGridEntity &sub_grid_entity = *it;

      EntityPointer host_entity_pointer = subDiscreteFunctionSpace_.gridPart().grid().template getHostEntity<0>( sub_grid_entity );

      const SubGridGeometry &sub_grid_geometry = sub_grid_entity.geometry();
      assert( sub_grid_entity.partitionType() == InteriorEntity );

      LocalMatrix local_matrix = global_matrix.localMatrix( sub_grid_entity, sub_grid_entity );

      const SubGridBaseFunctionSet &baseSet = local_matrix.domainBaseFunctionSet();
      const unsigned int numBaseFunctions = baseSet.numBaseFunctions();

      const IntersectionIterator iend = discreteFunctionSpace_.gridPart().iend( *host_entity_pointer );
      for( IntersectionIterator iit = discreteFunctionSpace_.gridPart().ibegin( *host_entity_pointer ); iit != iend; ++iit )
        {
          FaceQuadrature faceQuadrature( discreteFunctionSpace_.gridPart(), *iit, 2*subDiscreteFunctionSpace_.order()+ 1, FaceQuadrature::INSIDE);

          const FaceGeometryType& faceGeometry = iit->geometry();

//! spaeter loeschen:
#if 1

if ( iit->neighbor() )
  {
    EntityPointer outside_it = iit->outside();
    if ( subDiscreteFunctionSpace_.gridPart().grid().template contains<0>( *outside_it ) == true )
      { continue; }
  }
          // check if face
#endif

          const size_t numQuadraturePoints = faceQuadrature.nop();

          RangeType check_sum( 0.0 );

          for( size_t quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint )
            {

              const LocalCoordinate local_point = faceGeometry.local( faceQuadrature.point( quadraturePoint ) );

              const DomainType global_point = faceGeometry.global( local_point );

              // integration factors
              const double integrationFactor = faceGeometry.integrationElement( local_point );

              // weight
              const double quadratureWeight = faceQuadrature.weight(  quadraturePoint );

              check_sum += integrationFactor * quadratureWeight;

              for( unsigned int i = 0; i < numBaseFunctions; ++i )
                {
                  baseSet.evaluate( i, faceQuadrature[ quadraturePoint ], phi[ i ]);
                }


              for( unsigned int i = 0; i < numBaseFunctions; ++i )
                {

                 for( unsigned int j = 0; j < numBaseFunctions; ++j )
                   {

                     local_matrix.add( j, i, integrationFactor * quadratureWeight * (phi[ i ][ 0 ] * phi[ j ][ 0 ]) );

                   }

                }

            } // done loop over all quadrature points

          if ( check_sum != faceGeometry.volume())
           { std :: cout << "Error in Face Quadrature." << std :: endl; abort(); }

        }

    }


  }



  template< class SubGridDiscreteFunctionImp, class DiscreteFunctionImp, class DiffusionImp >
  void ConservativeFluxOperator< SubGridDiscreteFunctionImp, DiscreteFunctionImp, DiffusionImp >::printLocalRHS( SubGridDiscreteFunctionImp &rhs) const
    {

      typedef typename SubGridDiscreteFunctionImp::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
      typedef typename DiscreteFunctionSpaceType :: IteratorType IteratorType;
      typedef typename DiscreteFunctionImp::LocalFunctionType LocalFunctionType;

      const DiscreteFunctionSpaceType &discreteFunctionSpace    
        = rhs.space();

      const IteratorType endit = discreteFunctionSpace.end();
      for( IteratorType it = discreteFunctionSpace.begin(); it != endit; ++it )
      {

        LocalFunctionType elementOfRHS = rhs.localFunction( *it ); 

        const int numDofs = elementOfRHS.numDofs(); 
	for( int i = 0; i < numDofs; ++i )
        {
         std :: cout << "Number of Dof: " << i << " ; " << rhs.name() << " : " << elementOfRHS[ i ] << std :: endl;
        }

      }

    }  // end method


  template< class SubGridDiscreteFunctionImp, class DiscreteFunctionImp, class DiffusionImp >
  double ConservativeFluxOperator< SubGridDiscreteFunctionImp, DiscreteFunctionImp, DiffusionImp >::normRHS( SubGridDiscreteFunctionImp &rhs) const
    {

      double norm = 0.0;

      typedef typename SubGridDiscreteFunctionImp::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
      typedef typename DiscreteFunctionSpaceType :: IteratorType IteratorType;
      typedef typename IteratorType::Entity EntityType;
      typedef typename SubGridDiscreteFunctionImp::LocalFunctionType LocalFunctionType;
      typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
      typedef typename DiscreteFunctionSpaceType::GridType GridType;
      typedef typename GridType::template Codim<0>::Geometry
         EnGeometryType; 

      const DiscreteFunctionSpaceType &discreteFunctionSpace    
        = rhs.space();

      const IteratorType endit = discreteFunctionSpace.end();
      for( IteratorType it = discreteFunctionSpace.begin(); it != endit; ++it )
      {

        // entity
        const EntityType& entity = *it;

        // create quadrature for given geometry type 
        CachingQuadrature <GridPartType , 0 > quadrature(entity,2*discreteFunctionSpace.order()+2); 

        // get geoemetry of entity
        const EnGeometryType& geo = entity.geometry();

        LocalFunctionType localRHS = rhs.localFunction( *it ); 

        // integrate 
        const int quadratureNop = quadrature.nop();
        for(int quadraturePoint = 0; quadraturePoint < quadratureNop; ++quadraturePoint)
        {
          const double weight = quadrature.weight(quadraturePoint) * 
              geo.integrationElement(quadrature.point(quadraturePoint));

          RangeType value(0.0);
          localRHS.evaluate(quadrature[quadraturePoint],value);

          norm += weight * value * value;
        }

      }

     return norm;
     
    }  // end method



  // assemble the right hand side of the conservative flux problem 
  // ----------------------------------------------

  template< class SubGridDiscreteFunctionImp, class DiscreteFunctionImp, class DiffusionImp >
  //template< class MatrixType >
  void ConservativeFluxOperator< SubGridDiscreteFunctionImp, DiscreteFunctionImp, DiffusionImp >::assemble_RHS
       ( // direction 'e'
         JacobianRangeType &e_i,
         // solution of the local corrector problem
         const SubGridDiscreteFunction &local_corrector_e_i,
         // rhs flux problem:
         SubGridDiscreteFunction &rhs_flux_problem ) const
  {
    
    const SubGridDiscreteFunctionSpace &subDiscreteFunctionSpace = rhs_flux_problem.space();

    // set entries to zero:
    rhs_flux_problem.clear();

    // gradient of micro scale base function:
    std::vector< JacobianRangeType > gradient_phi( subDiscreteFunctionSpace.mapper().maxNumDofs() );

    RangeType rhs_L2_Norm = 0.0; 

    const SubGridIterator end = subDiscreteFunctionSpace.end();
    for( SubGridIterator it = subDiscreteFunctionSpace.begin(); it != end; ++it )
    {

      const SubGridEntity &local_grid_entity = *it;
      const SubGridGeometry &geometry = local_grid_entity.geometry();
      assert( local_grid_entity.partitionType() == InteriorEntity );

      SubGridLocalFunction elementOfRHS = rhs_flux_problem.localFunction( local_grid_entity );

      const SubGridBaseFunctionSet &baseSet = elementOfRHS.baseFunctionSet();
      const unsigned int numBaseFunctions = baseSet.numBaseFunctions();

      SubGridQuadrature quadrature( local_grid_entity, 2*subDiscreteFunctionSpace.order()+2 );
      const size_t numQuadraturePoints = quadrature.nop();
      for( size_t quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint )
      {

        const typename SubGridQuadrature::CoordinateType &local_point = quadrature.point( quadraturePoint );

        // remember, we are concerned with: - \int_{U(T)} (A^eps)(x) e · ∇ \phi(x)

        // global point in the subgrid
        DomainType global_point = geometry.global( local_point );

        const double weight = quadrature.weight( quadraturePoint ) * geometry.integrationElement( local_point );

        // transposed of the the inverse jacobian
        const FieldMatrix< double, dimension, dimension > &inverse_jac
          = geometry.jacobianInverseTransposed( local_point );


        // A^eps(x) ( e
        // diffusion operator evaluated in 'x' multiplied with e
        JacobianRangeType diffusion_in_e_i;
        diffusion_operator_.diffusiveFlux( global_point, e_i, diffusion_in_e_i );

        JacobianRangeType total_diffusive_flux;

        JacobianRangeType grad_corrector_e_i;
	SubGridLocalFunction localized_corrector_e_i = local_corrector_e_i.localFunction( local_grid_entity );
        localized_corrector_e_i.jacobian( quadrature[ quadraturePoint ], grad_corrector_e_i );
        diffusion_operator_.diffusiveFlux( global_point, grad_corrector_e_i, total_diffusive_flux );       

        total_diffusive_flux[ 0 ] += diffusion_in_e_i[ 0 ];

        for( unsigned int i = 0; i < numBaseFunctions; ++i )
        {
          // jacobian of the base functions, with respect to the reference element
          JacobianRangeType gradient_phi_ref_element;
          baseSet.jacobian( i, quadrature[ quadraturePoint ], gradient_phi_ref_element );

          // multiply it with transpose of jacobian inverse to obtain the jacobian with respect to the real entity
          inverse_jac.mv( gradient_phi_ref_element[ 0 ], gradient_phi[ i ][ 0 ] );
        }

        for( unsigned int i = 0; i < numBaseFunctions; ++i )
        {
          elementOfRHS[ i ] -= weight * (total_diffusive_flux[ 0 ] * gradient_phi[ i ][ 0 ]);

        }

      }

    } // end for-loop of 'Iterator'

  }

#endif

//! ------------------------------------------------------------------------------------------------
//! ------------------------------------------------------------------------------------------------



//! ------------------------------------------------------------------------------------------------



//! ------------------------------------------------------------------------------------------------
//! --------------------- the essential local msfem problem solver class ---------------------------



  template< class SubGridDiscreteFunctionType,
            class HostDiscreteFunctionType,
            class DiffusionOperatorType >
  class ConservativeFluxProblemSolver	
  {
  public:

    //! ---------------- typedefs for the HostDiscreteFunctionSpace -----------------------

    //! type of discrete function space
    typedef typename HostDiscreteFunctionType :: DiscreteFunctionSpaceType
      HostDiscreteFunctionSpaceType;

    //! type of (non-discrete )function space
    typedef typename HostDiscreteFunctionSpaceType :: FunctionSpaceType FunctionSpaceType;

    //! type of grid partition
    typedef typename HostDiscreteFunctionSpaceType :: GridPartType HostGridPartType;

    //! type of grid
    typedef typename HostDiscreteFunctionSpaceType :: GridType HostGridType;

    typedef typename HostGridType :: Traits :: LeafIndexSet HostGridLeafIndexSet;

    typedef typename HostDiscreteFunctionSpaceType :: IteratorType HostGridEntityIteratorType;

    typedef typename HostGridEntityIteratorType :: Entity HostEntityType;

    typedef typename HostEntityType :: EntityPointer HostEntityPointerType;
    
    typedef typename HostGridType ::template Codim<0> :: Geometry HostGridEntityGeometry;
        
    typedef typename HostDiscreteFunctionType :: LocalFunctionType HostLocalFunctionType;
    
    typedef typename HostGridPartType :: IntersectionIteratorType HostIntersectionIterator;


    //! ---------------- typedefs for the SubGridDiscreteFunctionSpace -----------------------
    //  ( typedefs for the local grid and the corresponding local ('sub') )discrete space ) 

    //! type of discrete function space
    typedef typename SubGridDiscreteFunctionType :: DiscreteFunctionSpaceType
      SubGridDiscreteFunctionSpaceType;

    typedef typename SubGridDiscreteFunctionSpaceType :: RangeFieldType RangeFieldType;
    typedef typename SubGridDiscreteFunctionSpaceType :: DomainType DomainType;
    typedef typename SubGridDiscreteFunctionSpaceType :: RangeType RangeType;
    typedef typename SubGridDiscreteFunctionSpaceType :: JacobianRangeType
      JacobianRangeType;

    //! type of grid partition
    typedef typename SubGridDiscreteFunctionSpaceType :: GridPartType SubGridPartType;

    //! type of grid
    typedef typename SubGridDiscreteFunctionSpaceType :: GridType SubGridType;


    typedef typename SubGridDiscreteFunctionSpaceType :: IteratorType SubGridIteratorType;
    
    typedef typename SubGridIteratorType :: Entity SubGridEntityType;
    
    typedef typename SubGridEntityType :: EntityPointer SubGridEntityPointerType;
    
    typedef typename SubGridDiscreteFunctionType :: LocalFunctionType SubGridLocalFunctionType;
    
    typedef typename SubGridDiscreteFunctionSpaceType :: LagrangePointSetType SubGridLagrangePointSetType;
    
//!-----------------------------------------------------------------------------------------



//! ------------------ Matrix Traits for the local Problems ---------------------

    enum { faceCodim = 1 };
    typedef typename SubGridLagrangePointSetType :: template Codim< faceCodim > :: SubEntityIteratorType SubGridFaceDofIteratorType;

    //! polynomial order of base functions
    enum { polynomialOrder = SubGridDiscreteFunctionSpaceType :: polynomialOrder };

    // flux problem matrix traits
    struct FluxProbMatrixTraits
     {
       typedef SubGridDiscreteFunctionSpaceType RowSpaceType;
       typedef SubGridDiscreteFunctionSpaceType ColumnSpaceType;
       typedef LagrangeMatrixSetup< false > StencilType;
       typedef ParallelScalarProduct< SubGridDiscreteFunctionSpaceType > ParallelScalarProductType;

       template< class M >
       struct Adapter
        {
          typedef LagrangeParallelMatrixAdapter< M > MatrixAdapterType;
        };
     };

    typedef SparseRowMatrixOperator< SubGridDiscreteFunctionType, SubGridDiscreteFunctionType, FluxProbMatrixTraits > FluxProbFEMMatrix;

    // OEMGMRESOp //OEMBICGSQOp // OEMBICGSTABOp
    typedef OEMBICGSTABOp< SubGridDiscreteFunctionType, FluxProbFEMMatrix > InverseFluxProbFEMMatrix;

  private:

    const DiffusionOperatorType& diffusion_;
    const HostDiscreteFunctionSpaceType& hostDiscreteFunctionSpace_;

  public:


    //! constructor - with diffusion operator A^{\epsilon}(x)
    ConservativeFluxProblemSolver( const HostDiscreteFunctionSpaceType &hostDiscreteFunctionSpace,
                                   const DiffusionOperatorType& diffusion_operator )
    : hostDiscreteFunctionSpace_( hostDiscreteFunctionSpace ),
      diffusion_( diffusion_operator )
     { }


    template < class Stream >
    void oneLinePrint( Stream& stream, const SubGridDiscreteFunctionType& func )
    {
      typedef typename SubGridDiscreteFunctionType::ConstDofIteratorType
         DofIteratorType;
      DofIteratorType it = func.dbegin();
      stream << "\n" << func.name() << ": [ ";
      for ( ; it != func.dend(); ++it )
         stream << std::setw(5) << *it << "  ";

      stream << " ] " << std::endl;
     }




    //! ----------- method: solve the local MsFEM problem ------------------------------------------

    void solve( JacobianRangeType &e_i, // direction 'e_i'
                const SubGridDiscreteFunctionType &local_corrector_e_i,
                SubGridDiscreteFunctionType &conservative_flux )
    {

      // set solution equal to zero:
      conservative_flux.clear();
      

      const SubGridDiscreteFunctionSpaceType &localDiscreteFunctionSpace = local_corrector_e_i.space();


      //! the matrix in our linear system of equations
      FluxProbFEMMatrix flux_prob_system_matrix( "Conservative Flux Problem System Matrix", localDiscreteFunctionSpace, localDiscreteFunctionSpace );

      //! define the discrete (elliptic) local MsFEM problem operator
      // ( effect of the discretized differential operator on a certain discrete function )
      // discrete elliptic operator describing the elliptic local msfem problems
      typedef ConservativeFluxOperator< SubGridDiscreteFunctionType, HostDiscreteFunctionType, DiffusionOperatorType > ConservativeFluxOperatorType;
      ConservativeFluxOperatorType cf_problem_operator( localDiscreteFunctionSpace, hostDiscreteFunctionSpace_, diffusion_ );

      const SubGridPartType &subGridPart = localDiscreteFunctionSpace.gridPart();
      const SubGridType &subGrid = localDiscreteFunctionSpace.grid();

      //! right hand side vector of the algebraic local MsFEM problem
      SubGridDiscreteFunctionType rhs( "RHS of Conservative Flux Problem", localDiscreteFunctionSpace );
      rhs.clear();

      // assemble the stiffness matrix
      cf_problem_operator.assemble_matrix( flux_prob_system_matrix );

      // assemble right hand side of algebraic local msfem problem
      cf_problem_operator.assemble_RHS( e_i, local_corrector_e_i, rhs );

      // oneLinePrint( std::cout , rhs );

      const double norm_rhs = cf_problem_operator.normRHS( rhs );

      if ( !( rhs.dofsValid() ) )
        { std :: cout << "Local Flux Problem RHS invalid." << std :: endl;
          abort(); }

      if ( norm_rhs < /*1e-06*/ 1e-30 )
        {
          conservative_flux.clear();
          std :: cout << "Local Flux Problem with solution zero." << std :: endl;
        }
      else
        {
          InverseFluxProbFEMMatrix flux_prob_biCGStab( flux_prob_system_matrix, 1e-8, 1e-8, 20000, FLUX_SOLVER_VERBOSE );
          flux_prob_biCGStab( rhs, conservative_flux );
        }

     if ( !(conservative_flux.dofsValid()) )
       {
         std::cout << "Solution of the Local Flux Problem is invalid!" << std::endl;
         std :: abort();
       }

     oneLinePrint( std::cout , conservative_flux );
    
    }

    //! ----------- end method: solve local MsFEM problem ------------------------------------------



    // create a hostgrid function from a subgridfunction
    void subgrid_to_hostrid_function( const SubGridDiscreteFunctionType &sub_func,
                                            HostDiscreteFunctionType &host_func )
    {
      
       host_func.clear();


       const SubGridDiscreteFunctionSpaceType &subDiscreteFunctionSpace = sub_func.space();
       const SubGridType &subGrid = subDiscreteFunctionSpace.grid();

       SubgridIteratorType sub_endit = subDiscreteFunctionSpace.end();
       for( SubgridIteratorType sub_it = subDiscreteFunctionSpace.begin(); sub_it != sub_endit; ++sub_it )
          {
#if 0
             const SubgridEntityType &sub_entity = *sub_it;

             HostEntityPointerType host_entity_pointer = subGrid.template getHostEntity<0>( *sub_it );
             const HostEntityType& host_entity = *host_entity_pointer;

             SubLocalFunctionType sub_loc_value = sub_func.localFunction( sub_entity );
             HostLocalFunctionType host_loc_value = host_func.localFunction( host_entity );

             const unsigned int numBaseFunctions = sub_loc_value.baseFunctionSet().numBaseFunctions();
             for( unsigned int i = 0; i < numBaseFunctions; ++i )
               {
                 host_loc_value[ i ] = sub_loc_value[ i ];
               }
#endif
          }

    }

#if 0     
    // method for solving and saving the solutions of the local msfem problems
    // for the whole set of macro-entities and for every unit vector e_i

    //! ---- method: solve and save the whole set of local msfem problems -----

    // Use the host-grid entities of Level 'computational_level' as computational domains for the subgrid computations
    void assemble_all( bool silent = true /* state information on subgrids */ )
    {
 
      enum { dimension = GridType :: dimension};
      enum { maxnumOfBaseFct = 100 }; 

      JacobianRangeType e[dimension];
      for( int i = 0; i < dimension ; ++i )
        for( int j = 0; j < dimension ; ++j )
         {
           if ( i == j )
             { e[i][0][j] = 1.0; }
           else
             { e[i][0][j] = 0.0; }
         }

      const HostGridPartType& hostGridPart = hostDiscreteFunctionSpace_.gridPart();

      HostGridType& hostGrid = hostDiscreteFunctionSpace_.gridPart().grid();

      // number of coarse grid entities (of codim 0).
      int number_of_coarse_grid_entities = specifier_.getNumOfCoarseEntities();

      std :: cout << "in method 'assemble_all': number_of_coarse_grid_entities = " << number_of_coarse_grid_entities << std :: endl;
      
      // --------------- writing data output ---------------------
      // typedefs and initialization
      #ifdef VTK_OUTPUT

      typedef Tuple<HostDiscreteFunctionType*> IOTupleType;
      typedef DataOutput<HostGridType, IOTupleType> DataOutputType;

      // general output parameters
      LocalProblemDataOutputParameters outputparam;
      outputparam.set_path( path_ );

      // sequence stamp
      std::stringstream outstring;
      // -------------------------------------------------------
      #endif

      long double starting_time = clock();

      // we want to determine minimum, average and maxiumum time for solving a local msfem problem in the current method
      double minimum_time_c_p = 1000000;
      double average_time_c_p = 0;
      double maximum_time_c_p = 0;

      HostDiscreteFunctionSpaceType& coarseSpace = specifier_.coarseSpace();

      const HostGridLeafIndexSet& coarseGridLeafIndexSet = coarseSpace.gridPart().grid().leafIndexSet();     
      
      for( HostGridEntityIteratorType coarse_it = coarseSpace.begin(); coarse_it != coarseSpace.end(); ++coarse_it )
        {

          int coarse_index = coarseGridLeafIndexSet.index( *coarse_it );

          #if 0
          std :: cout << "coarse_it->geometry().corner(0) = " << coarse_it->geometry().corner(0) << std :: endl;
          std :: cout << "coarse_it->geometry().corner(1) = " << coarse_it->geometry().corner(1) << std :: endl;
          std :: cout << "coarse_it->geometry().corner(2) = " << coarse_it->geometry().corner(2) << std :: endl;
          #endif

          bool writer_is_open = false;

	  char location_lps[50];
          sprintf( location_lps, "_localProblemSolutions_%d", coarse_index );
          std::string location_lps_s( location_lps );

          std :: string locprob_solution_location = path_ + location_lps_s;

          DiscreteFunctionWriter dfw( (locprob_solution_location).c_str() );

          writer_is_open = dfw.open();

          if ( writer_is_open )
          {

	    SubGridType& subGrid = subgrid_list_.getSubGrid( coarse_index );
            SubGridPartType subGridPart( subGrid );

            SubDiscreteFunctionSpaceType subDiscreteFunctionSpace( subGridPart );

	    char name_loc_sol[50];
            sprintf( name_loc_sol, "Local Problem Solution %d", coarse_index );
            std::string name_local_solution( name_loc_sol );

	    //! only for dimension 2!
            SubDiscreteFunctionType local_problem_solution_0( name_local_solution, subDiscreteFunctionSpace );
            local_problem_solution_0.clear();

            SubDiscreteFunctionType local_problem_solution_1( name_local_solution, subDiscreteFunctionSpace );
            local_problem_solution_1.clear();

            std :: cout << "Number of the local problem: " << dimension * coarse_index << " (of " << (dimension*number_of_coarse_grid_entities)-1 << " problems in total)" << std :: endl;

            // take time
            long double time_now = clock();

	    // solve the problems
	    solvelocalproblem( e[0], local_problem_solution_0 );

	    // min/max time
            if ( (clock()-time_now)/CLOCKS_PER_SEC > maximum_time_c_p )
               { maximum_time_c_p = (clock()-time_now)/CLOCKS_PER_SEC; }
            if ( (clock()-time_now)/CLOCKS_PER_SEC < minimum_time_c_p )
               { minimum_time_c_p = (clock()-time_now)/CLOCKS_PER_SEC; }

            std :: cout << "Number of the local problem: " << (dimension*coarse_index)+1 << " (of " << (dimension*number_of_coarse_grid_entities)-1 << " problems in total)" << std :: endl;
	    
            // take time
            time_now = clock();
	    
	    // solve the problems
	    solvelocalproblem( e[1], local_problem_solution_1 );
	    
	    // min/max time
            if ( (clock()-time_now)/CLOCKS_PER_SEC > maximum_time_c_p )
               { maximum_time_c_p = (clock()-time_now)/CLOCKS_PER_SEC; }
            if ( (clock()-time_now)/CLOCKS_PER_SEC < minimum_time_c_p )
               { minimum_time_c_p = (clock()-time_now)/CLOCKS_PER_SEC; }

            dfw.append( local_problem_solution_0);

            dfw.append( local_problem_solution_1);

            // oneLinePrint( std::cout , local_problem_solution_0 );
            // oneLinePrint( std::cout , local_problem_solution_1 );

	    HostDiscreteFunctionType host_local_solution( name_local_solution, hostDiscreteFunctionSpace_ );
            subgrid_to_hostrid_function( local_problem_solution_0, host_local_solution );
    
            // --------------- writing data output ---------------------
            // (writing)
            #ifdef VTK_OUTPUT
	    
            // --------- data output local solution --------------

            // create and initialize output class
            IOTupleType local_solution_series_0( &host_local_solution );

            char ls_name_0[50];
            sprintf( ls_name_0, "local_problem_solution_e0_%d", coarse_index );
            std::string ls_name_0_s( ls_name_0 );

            outputparam.set_prefix( ls_name_0_s );
            DataOutputType localsol_dataoutput_0( hostDiscreteFunctionSpace_.gridPart().grid(), local_solution_series_0, outputparam );

            // write data
            outstring << "local-problem-solution-0";
            localsol_dataoutput_0.writeData( 1.0 /*dummy*/, outstring.str() );
            // clear the std::stringstream:
            outstring.str(std::string());

            // -------------------------------------------------------
            #endif

            subgrid_to_hostrid_function( local_problem_solution_1, host_local_solution );

            // --------------- writing data output ---------------------
            // (writing)
            #ifdef VTK_OUTPUT

            // --------- data output local solution --------------

            // create and initialize output class
            IOTupleType local_solution_series_1( &host_local_solution );

            char ls_name_1[50];
            sprintf( ls_name_1, "local_problem_solution_e1_%d", coarse_index );
            std::string ls_name_1_s( ls_name_1 );

            outputparam.set_prefix( ls_name_1_s );
            DataOutputType localsol_dataoutput_1( hostDiscreteFunctionSpace_.gridPart().grid(), local_solution_series_1, outputparam );

            // write data
            outstring << "local-problem-solution-1";
            localsol_dataoutput_1.writeData( 1.0 /*dummy*/, outstring.str() );
            // clear the std::stringstream:
            outstring.str(std::string());

            // -------------------------------------------------------
            #endif	  	 

          }


      } // end: 'if ( writer_is_open )'


      if (data_file_)
      {
         if (data_file_->is_open())
         {
           (*data_file_) << std :: endl;
           (*data_file_) << "In method: assemble_all." << std :: endl << std :: endl;
           (*data_file_) << "MsFEM problems solved for " << number_of_coarse_grid_entities << " coarse grid entities." << std :: endl;
           (*data_file_) << dimension*number_of_coarse_grid_entities << " local MsFEM problems solved in total." << std :: endl;
           (*data_file_) << "Minimum time for solving a local problem = " << minimum_time_c_p << "s." << std :: endl;
           (*data_file_) << "Maximum time for solving a localproblem = " << maximum_time_c_p << "s." << std :: endl;
           (*data_file_) << "Average time for solving a localproblem = " << ((clock()-starting_time)/CLOCKS_PER_SEC)/(dimension*number_of_coarse_grid_entities) << "s." << std :: endl;
           (*data_file_) << "Total time for computing and saving the localproblems = " << ((clock()-starting_time)/CLOCKS_PER_SEC) << "s," << std :: endl << std :: endl;
         }
      }


    }

#endif

 }; //end class




} // end namespace Dune

#endif // #ifndef DiscreteEllipticMsFEMLocalProblem_HH
