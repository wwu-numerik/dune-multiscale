#ifndef DiscreteEllipticMsFEMLocalProblem_HH
#define DiscreteEllipticMsFEMLocalProblem_HH

#include <dune/common/fmatrix.hh>

#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/operator/common/operator.hh>

// artificical mass coefficient to guarantee uniqueness and existence of the cell problem solution
//  (should be as small as possible)
#define CELL_MASS_WEIGHT 0.0000001

// CELLSOLVER_VERBOSE: 0 = false, 1 = true
#define CELLSOLVER_VERBOSE false

#include <dune/fem/operator/2order/lagrangematrixsetup.hh>

#include <dune/multiscale/operators/disc_func_writer/discretefunctionwriter.hh>

namespace Dune
{

  // define output traits
  struct CellProblemDataOutputParameters : public DataOutputParameters {

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
        return "data_output_hmm";
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
  template< class PeriodicDiscreteFunctionImp, class DiffusionImp >
  class DiscreteCellProblemOperator
  : public Operator< typename PeriodicDiscreteFunctionImp::RangeFieldType, typename PeriodicDiscreteFunctionImp::RangeFieldType, PeriodicDiscreteFunctionImp, PeriodicDiscreteFunctionImp >
  {
    typedef DiscreteCellProblemOperator< PeriodicDiscreteFunctionImp, DiffusionImp > This;

  public:
    typedef PeriodicDiscreteFunctionImp DiscreteFunction;
    typedef DiffusionImp DiffusionModel;

    typedef typename DiscreteFunction::DiscreteFunctionSpaceType DiscreteFunctionSpace;

    typedef typename DiscreteFunctionSpace::GridPartType GridPart;
    typedef typename DiscreteFunctionSpace::GridType GridType;
    typedef typename DiscreteFunctionSpace::RangeFieldType RangeFieldType;

    typedef typename DiscreteFunctionSpace::DomainType DomainType;
    typedef typename DiscreteFunctionSpace::RangeType RangeType;
    typedef typename DiscreteFunctionSpace::JacobianRangeType
      JacobianRangeType;

  protected:
    static const int dimension = GridPart::GridType::dimension;
    static const int polynomialOrder = DiscreteFunctionSpace::polynomialOrder;

    typedef typename DiscreteFunction::LocalFunctionType LocalFunction;

    typedef typename DiscreteFunctionSpace::BaseFunctionSetType BaseFunctionSet;
    typedef typename DiscreteFunctionSpace::LagrangePointSetType LagrangePointSet;
    typedef typename LagrangePointSet::template Codim< 1 >::SubEntityIteratorType FaceDofIterator;

    typedef typename DiscreteFunctionSpace::IteratorType Iterator;
    typedef typename Iterator::Entity Entity;
    typedef typename Entity::Geometry Geometry;

    typedef typename GridPart::IntersectionIteratorType IntersectionIterator;
    typedef typename IntersectionIterator::Intersection Intersection;

    typedef CachingQuadrature< GridPart, 0 > Quadrature;

  public:
    DiscreteCellProblemOperator( const DiscreteFunctionSpace &periodicDiscreteFunctionSpace, const DiffusionModel &diffusion_op )
    : periodicDiscreteFunctionSpace_( periodicDiscreteFunctionSpace ),
      diffusion_operator_( diffusion_op )
    {}
        
  private:
    DiscreteCellProblemOperator ( const This & );

  public:

    // dummy operator
    virtual void
    operator() ( const DiscreteFunction &u, DiscreteFunction &w ) const;

    template< class MatrixType >
    void assemble_matrix ( const DomainType &x_T, MatrixType &global_matrix ) const;

    // the right hand side assembler methods
    void assembleCellRHS_linear ( //the global quadrature point in the macro grid element T
                                  const DomainType &x_T,
                                  //\nabla_x \Phi_H(x_T) (the coarse function to reconstruct):
                                  JacobianRangeType &grad_coarse_function,
                                  // rhs cell problem:
                                  DiscreteFunction &cell_problem_RHS ) const;

    void printCellRHS( DiscreteFunction &rhs) const;

    double normRHS( DiscreteFunction &rhs) const;

  private:
    const DiscreteFunctionSpace &periodicDiscreteFunctionSpace_;
    const DiffusionModel &diffusion_operator_;
  };


  // dummy implementation of "operator()"
  // 'w' = effect of the discrete operator on 'u'
  template< class DiscreteFunctionImp, class DiffusionImp >
  void DiscreteCellProblemOperator< DiscreteFunctionImp, DiffusionImp >::operator() ( const DiscreteFunction &u, DiscreteFunction &w ) const 
  {

    std :: cout << "the ()-operator of the DiscreteCellProblemOperator class is not yet implemented and still a dummy." << std :: endl;
    std :: abort();

  }


  //! stiffness matrix for a linear elliptic diffusion operator 
  // we obtain entries of the following kind
  // (cell problem for the macro grid element 'T' and for the base-function '\Phi_H',
  //  x_T denotes the barycenter of T, \delta denotes the cell size )
  //  \int_Y A_h^{\eps}(t,x_T + \delta*y) \nabla phi_h_i(y) \cdot \nabla phi_h_j(y)
  //    + CELL_MASS_WEIGHT * \int_Y phi_h_i(y) \phi_h_j(y)
  //    (the second summand yields an artificical mass term to guarantee uniqueness and existence
  //     for the problem with periodic boundary condition.
  //     This is an alternative to the 'average zero' condition.)
  template< class PeriodicDiscreteFunctionImp, class DiffusionImp >
  template< class MatrixType >
  void DiscreteCellProblemOperator< PeriodicDiscreteFunctionImp, DiffusionImp >::assemble_matrix (const DomainType &x_T, MatrixType &global_matrix ) const
  // x_T is the barycenter of the macro grid element T
  {
    typedef typename MatrixType::LocalMatrixType LocalMatrix;

    Problem::ModelProblemData model_info;
    const double delta = model_info.getDelta();

    global_matrix.reserve();
    global_matrix.clear();


    // micro scale base function:
    std::vector< RangeType > phi( periodicDiscreteFunctionSpace_.mapper().maxNumDofs() );

    // gradient of micro scale base function:
    std::vector< typename BaseFunctionSet::JacobianRangeType > gradient_phi( periodicDiscreteFunctionSpace_.mapper().maxNumDofs() );

    const Iterator end = periodicDiscreteFunctionSpace_.end();
    for( Iterator it = periodicDiscreteFunctionSpace_.begin(); it != end; ++it )
    {

      const Entity &cell_grid_entity = *it;
      const Geometry &cell_grid_geometry = cell_grid_entity.geometry();
      assert( cell_grid_entity.partitionType() == InteriorEntity );

      LocalMatrix local_matrix = global_matrix.localMatrix( cell_grid_entity, cell_grid_entity );

      const BaseFunctionSet &baseSet = local_matrix.domainBaseFunctionSet();
      const unsigned int numBaseFunctions = baseSet.numBaseFunctions();

      // for constant diffusion "2*discreteFunctionSpace_.order()" is sufficient, for the general case, it is better to use a higher order quadrature:
      Quadrature quadrature( cell_grid_entity, 2*periodicDiscreteFunctionSpace_.order()+2 );
      const size_t numQuadraturePoints = quadrature.nop();
      for( size_t quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint )
      {
        // local (barycentric) coordinates (with respect to cell grid entity)
        const typename Quadrature::CoordinateType &local_point = quadrature.point( quadraturePoint );

        // global point in the unit cell Y
        DomainType global_point = cell_grid_geometry.global( local_point );

        // x_T + (delta * global_point)
        DomainType x_T_delta_global_point;
        for( int k = 0; k < dimension; ++k )
          {
            x_T_delta_global_point[ k ] = x_T[ k ] + (delta * global_point[ k ]);
          }

        const double weight = quadrature.weight( quadraturePoint ) *
             cell_grid_geometry.integrationElement( local_point );

        // transposed of the the inverse jacobian
        const FieldMatrix< double, dimension, dimension > &inverse_jac
          = cell_grid_geometry.jacobianInverseTransposed( local_point );

        for( unsigned int i = 0; i < numBaseFunctions; ++i )
        {
          // jacobian of the base functions, with respect to the reference element
          typename BaseFunctionSet::JacobianRangeType gradient_phi_ref_element;
          baseSet.jacobian( i, quadrature[ quadraturePoint ], gradient_phi_ref_element );

          // multiply it with transpose of jacobian inverse to obtain the jacobian with respect to the real entity
          inverse_jac.mv( gradient_phi_ref_element[ 0 ], gradient_phi[ i ][ 0 ] );

          baseSet.evaluate( i, quadrature[ quadraturePoint ], phi[ i ]);

        }

        for( unsigned int i = 0; i < numBaseFunctions; ++i )
        {
          // A( x_T + \delta y, \nabla \phi )
          // diffusion operator evaluated in (x_T + \delta y , \nabla \phi)
          typename LocalFunction::JacobianRangeType diffusion_in_gradient_phi;
          diffusion_operator_.diffusiveFlux( x_T_delta_global_point , gradient_phi[ i ], diffusion_in_gradient_phi );
          for( unsigned int j = 0; j < numBaseFunctions; ++j )
            {
              // stiffness contribution
              local_matrix.add( j, i, weight * (diffusion_in_gradient_phi[ 0 ] * gradient_phi[ j ][ 0 ]) );

              // mass contribution
              local_matrix.add( j, i, CELL_MASS_WEIGHT * weight * (phi[ i ][ 0 ] * phi[ j ][ 0 ]) );

            }
        }
      }
    }

  }


#if 1
  template< class DiscreteFunctionImp, class DiffusionImp >
  void DiscreteCellProblemOperator< DiscreteFunctionImp, DiffusionImp >::printCellRHS( DiscreteFunctionImp &rhs) const
    {

      typedef typename DiscreteFunctionImp::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
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


  template< class DiscreteFunctionImp, class DiffusionImp >
  double DiscreteCellProblemOperator< DiscreteFunctionImp, DiffusionImp >::normRHS( DiscreteFunctionImp &rhs) const
    {

      double norm = 0.0;

      typedef typename DiscreteFunctionImp::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
      typedef typename DiscreteFunctionSpaceType :: IteratorType IteratorType;
      typedef typename IteratorType::Entity EntityType;
      typedef typename DiscreteFunctionImp::LocalFunctionType LocalFunctionType;
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
#endif

#if 1
  // assemble the right hand side of a cell problem
  // ----------------------------------------------

  // assemble method for the case of a linear diffusion operator
  // (in this case, no Newton method is required, which is why there is no dependency on an old fine-scale discrete function / old iteration step )

  // we compute the following entries for each fine-scale base function phi_h_i:
  // - \int_Y A^{\eps}( x_T + \delta*y ) \nabla_x PHI_H(x_T) \cdot \nabla_y phi_h_i(y)
  template< class DiscreteFunctionImp, class DiffusionImp >
  //template< class MatrixType >
  void DiscreteCellProblemOperator< DiscreteFunctionImp, DiffusionImp >::assembleCellRHS_linear
       ( // the global quadrature point in the macro grid element T
         const DomainType &x_T,
         // \nabla_x \Phi_H(x_T):
         JacobianRangeType &gradient_PHI_H,
         // rhs cell problem:
         DiscreteFunctionImp &cell_problem_RHS ) const
  {


    typedef typename DiscreteFunction::DiscreteFunctionSpaceType DiscreteFunctionSpace;
    typedef typename DiscreteFunction::LocalFunctionType LocalFunction;

    typedef typename DiscreteFunctionSpace::BaseFunctionSetType BaseFunctionSet;
    typedef typename DiscreteFunctionSpace::IteratorType Iterator;
    typedef typename Iterator::Entity Entity;
    typedef typename Entity::Geometry Geometry;

    typedef typename DiscreteFunctionSpace::GridPartType GridPart;
    typedef CachingQuadrature< GridPart, 0 > Quadrature;

    const DiscreteFunctionSpace &discreteFunctionSpace = cell_problem_RHS.space();

    // set entries to zero:
    cell_problem_RHS.clear();

    // model problem data:
    Problem::ModelProblemData problem_info;

    // get edge length of cell:
    const double delta = problem_info.getDelta();

    // gradient of micro scale base function:
    std::vector< JacobianRangeType > gradient_phi( discreteFunctionSpace.mapper().maxNumDofs() );

    RangeType rhs_L2_Norm = 0.0; 

    const Iterator end = discreteFunctionSpace.end();
    for( Iterator it = discreteFunctionSpace.begin(); it != end; ++it )
    {

      const Entity &cell_grid_entity = *it;
      const Geometry &geometry = cell_grid_entity.geometry();
      assert( cell_grid_entity.partitionType() == InteriorEntity );

      LocalFunction elementOfRHS = cell_problem_RHS.localFunction( cell_grid_entity );

      const BaseFunctionSet &baseSet = elementOfRHS.baseFunctionSet();
      const unsigned int numBaseFunctions = baseSet.numBaseFunctions();

      Quadrature quadrature( cell_grid_entity, 2*discreteFunctionSpace.order()+2 );
      const size_t numQuadraturePoints = quadrature.nop();
      for( size_t quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint )
      {

        const typename Quadrature::CoordinateType &local_point = quadrature.point( quadraturePoint );

        // global point in the unit cell Y:
        DomainType global_point = geometry.global( local_point );

        // x_T + (delta * global_point)
        DomainType x_T_delta_global_point;
        for( int k = 0; k < dimension; ++k )
          {
            x_T_delta_global_point[ k ] = x_T[ k ] + (delta * global_point[ k ]);
          }

        const double weight = quadrature.weight( quadraturePoint ) * geometry.integrationElement( local_point );

        // transposed of the the inverse jacobian
        const FieldMatrix< double, dimension, dimension > &inverse_jac
          = geometry.jacobianInverseTransposed( local_point );


        // A^{\eps}( x_T + \delta y) \nabla_x PHI_H(x_T)
        // diffusion operator evaluated in (x_T + \delta y) multiplied with \nabla_x PHI_H(x_T)
        JacobianRangeType diffusion_in_gradient_PHI_H;
        diffusion_operator_.diffusiveFlux( x_T_delta_global_point, gradient_PHI_H, diffusion_in_gradient_PHI_H );

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
          elementOfRHS[ i ] -= weight * (diffusion_in_gradient_PHI_H[ 0 ] * gradient_phi[ i ][ 0 ]);
        }

      }
    }

  }
#endif


//! ------------------------------------------------------------------------------------------------
//! ------------------------------------------------------------------------------------------------





//! ------------------------------------------------------------------------------------------------
//! ---------------- the cell problem numbering manager classes ------------------------------------

// comparison class for the CellProblemNumberingManager:
template< class GridPartType, class DomainType, class EntityPointerType >
struct classcomp {

  bool operator() (const std::pair<EntityPointerType, int>& left_entity_pair,
                   const std::pair<EntityPointerType, int>& right_entity_pair) const
  {

    // compare the barycenteres of the entities with the lexicographic order, than compare the int's (number of local base function)

    typedef CachingQuadrature< GridPartType, 0 > Quadrature;



    // ------ right element

    const typename EntityPointerType::Entity::Geometry &geometry_right = (*(right_entity_pair.first)).geometry();

    Quadrature quadrature_right( (*(right_entity_pair.first)), 0 );

    // local barycenter (with respect to entity)
    const typename Quadrature::CoordinateType &local_point_right = quadrature_right.point(0);

    DomainType barycenter_right_entity = geometry_right.global( local_point_right );


    // ------ left element

    const typename EntityPointerType::Entity::Geometry &geometry_left = (*(left_entity_pair.first)).geometry();

    Quadrature quadrature_left( (*(left_entity_pair.first)), 0 );

    // local barycenter (with respect to entity)
    const typename Quadrature::CoordinateType &local_point_left = quadrature_left.point(0);

    DomainType barycenter_left_entity = geometry_left.global( local_point_left );



    enum { dimension = GridPartType::GridType::dimension};

    int current_axis = dimension-1;

    while ( current_axis >= 0 )
     {

	if ( barycenter_left_entity[ current_axis ] < barycenter_right_entity[ current_axis ] )
	    { return true; }
	else if ( barycenter_left_entity[ current_axis ] > barycenter_right_entity[ current_axis ] )
	    { return false; }

	current_axis -= 1;

      }

    if ( left_entity_pair.second < right_entity_pair.second )
      {
        return true; }
    else
      { return false; }

    return true;
   }
};



// comparison class for the CellProblemNumberingManager (just comparison of two entities!)
template< class GridPartType, class DomainType, class EntityPointerType >
struct entity_compare {

  bool operator() ( EntityPointerType left_entity,
                    EntityPointerType right_entity ) const
  {

    // compare the barycenteres of the entities with the lexicographic order

    typedef CachingQuadrature< GridPartType, 0 > Quadrature;



    // ------ right element

    const typename EntityPointerType::Entity::Geometry &geometry_right = (*right_entity).geometry();

    Quadrature quadrature_right( *right_entity , 0 );

    // local barycenter (with respect to entity)
    const typename Quadrature::CoordinateType &local_point_right = quadrature_right.point(0);

    DomainType barycenter_right_entity = geometry_right.global( local_point_right );


    // ------ left element

    const typename EntityPointerType::Entity::Geometry &geometry_left = (*left_entity).geometry();

    Quadrature quadrature_left( *left_entity , 0 );

    // local barycenter (with respect to entity)
    const typename Quadrature::CoordinateType &local_point_left = quadrature_left.point(0);

    DomainType barycenter_left_entity = geometry_left.global( local_point_left );



    enum { dimension = GridPartType::GridType::dimension};

    int current_axis = dimension-1;

    while ( current_axis >= 0 )
     {

	if ( barycenter_left_entity[ current_axis ] < barycenter_right_entity[ current_axis ] )
	    { return true; }
	else if ( barycenter_left_entity[ current_axis ] > barycenter_right_entity[ current_axis ] )
	    { return false; }

	current_axis -= 1;

      }

    return false;
   }
};



// only for the combination entity + number of local base function on entity
template< class DiscreteFunctionSpaceType >
class CellProblemNumberingManager{

public:

  typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;
  typedef typename GridPartType :: GridType GridType;

  typedef typename GridType :: template Codim<0> :: Entity EntityType; 
  typedef typename GridType :: template  Codim<0> :: EntityPointer EntityPointerType; 

  typedef typename DiscreteFunctionSpaceType :: BaseFunctionSetType BaseFunctionSetType;
  typedef typename DiscreteFunctionSpaceType :: IteratorType IteratorType;

  typedef typename DiscreteFunctionSpaceType :: DomainType DomainType;

  typedef classcomp< GridPartType, DomainType, EntityPointerType > CompClass;

  typedef std::map< std::pair<EntityPointerType, int> , int, CompClass > CellNumMapType;

  // for the comparison of two entities:
  typedef entity_compare< GridPartType, DomainType, EntityPointerType > CompEntityClass;

  typedef std::map< EntityPointerType , int, CompEntityClass > CellNumMapNLType;

  CellNumMapType *cell_numbering_map_;
  CellNumMapNLType *cell_numbering_map_NL_;

  // simpliefied: in general we need CellNumMapType for the cell problem numering in the linear setting (entity and local number of base function) and in the nonlinear case we need CellNumMapNLType (NL stands for nonlinear).
  // CellNumMapType is also required in the nonlinear case if we use the standard MsFEM formulation (no PGF)

  inline explicit CellProblemNumberingManager ( DiscreteFunctionSpaceType &discreteFunctionSpace)
    {

       cell_numbering_map_ = new CellNumMapType;
       cell_numbering_map_NL_ = new CellNumMapNLType;

       int counter = 0;
       int number_of_entity = 0;

       IteratorType endit = discreteFunctionSpace.end();
       for(IteratorType it = discreteFunctionSpace.begin(); it != endit ; ++it )
         {

           cell_numbering_map_NL_->insert( std::make_pair( it , number_of_entity ) );

           const BaseFunctionSetType baseSet
              = discreteFunctionSpace.baseFunctionSet( *it );

           // number of base functions on entity
           const int numBaseFunctions = baseSet.numBaseFunctions();

           for( int i = 0; i < numBaseFunctions; ++i )
             {
               std::pair<EntityPointerType, int>  idPair( it , i );
               cell_numbering_map_->insert( std::make_pair( idPair , counter ) );
               counter++;
             }

           number_of_entity++;

         }
    }

  // use 'cp_num_manager.get_number_of_cell_problem( it, i )'
  inline int get_number_of_cell_problem ( EntityPointerType &ent, const int &numOfBaseFunction ) const
   {
     std::pair<EntityPointerType, int>  idPair( ent , numOfBaseFunction );
     return (*cell_numbering_map_)[idPair];
   }

  // use 'cp_num_manager.get_number_of_cell_problem( it )'
  // Note: 'get_number_of_cell_problem( it )' is NOT equal to 'get_number_of_cell_problem( it , 0 )'!
  inline int get_number_of_cell_problem ( EntityPointerType &ent ) const
   {
     return (*cell_numbering_map_NL_)[ent];
   }

};

//! ------------------ end of the cell problem numbering manager classes ---------------------------
//! ------------------------------------------------------------------------------------------------



//! ------------------------------------------------------------------------------------------------



//! ------------------------------------------------------------------------------------------------
//! --------------------- the essential cell problem solver class ----------------------------------

  template< class PeriodicDiscreteFunctionImp, class DiffusionOperatorImp >
  class CellProblemSolver	
  {
  public:

    //! type of discrete functions
    typedef PeriodicDiscreteFunctionImp PeriodicDiscreteFunctionType;

    //! type of discrete function space
    typedef typename PeriodicDiscreteFunctionType :: DiscreteFunctionSpaceType
      PeriodicDiscreteFunctionSpaceType;

    //! type of grid partition
    typedef typename PeriodicDiscreteFunctionSpaceType :: GridPartType PeriodicGridPartType;

    //! type of grid
    typedef typename PeriodicDiscreteFunctionSpaceType :: GridType PeriodicGridType;

    //! type of range vectors
    typedef typename PeriodicDiscreteFunctionSpaceType :: RangeType RangeType;

    //! type of range vectors
    typedef typename PeriodicDiscreteFunctionSpaceType :: DomainType DomainType;

    //! polynomial order of base functions
    enum { polynomialOrder = PeriodicDiscreteFunctionSpaceType :: polynomialOrder };

    //! type of the (possibly non-linear) diffusion operator
    typedef DiffusionOperatorImp DiffusionType;

    struct CellMatrixTraits
     {
       typedef PeriodicDiscreteFunctionSpaceType RowSpaceType;
       typedef PeriodicDiscreteFunctionSpaceType ColumnSpaceType;
       typedef LagrangeMatrixSetup< false > StencilType;
       typedef ParallelScalarProduct< PeriodicDiscreteFunctionSpaceType > ParallelScalarProductType;

       template< class M >
       struct Adapter
        {
          typedef LagrangeParallelMatrixAdapter< M > MatrixAdapterType;
        };
     };

    typedef SparseRowMatrixOperator< PeriodicDiscreteFunctionType, PeriodicDiscreteFunctionType, CellMatrixTraits > CellFEMMatrix;

    // OEMGMRESOp //OEMBICGSQOp // OEMBICGSTABOp
    typedef OEMBICGSTABOp< PeriodicDiscreteFunctionType, CellFEMMatrix > InverseCellFEMMatrix;

    // discrete elliptic operator describing the elliptic cell problems
    typedef DiscreteCellProblemOperator< PeriodicDiscreteFunctionType, DiffusionType > CellProblemOperatorType;


  private:

    const PeriodicDiscreteFunctionSpaceType &periodicDiscreteFunctionSpace_; //Referenz &, wenn & verwendet, dann unten:
    DiffusionType &diffusion_;

    std :: ofstream *data_file_;

  public:

    //! constructor - with diffusion operator A^{\epsilon}(x)
    CellProblemSolver( const PeriodicDiscreteFunctionSpaceType &periodicDiscreteFunctionSpace,
                             DiffusionType &diffusion_operator )
    : periodicDiscreteFunctionSpace_( periodicDiscreteFunctionSpace ),
      diffusion_( diffusion_operator ),
      data_file_( NULL )
    {
    }

    //! constructor - with diffusion operator A^{\epsilon}(x)
    CellProblemSolver( const PeriodicDiscreteFunctionSpaceType &periodicDiscreteFunctionSpace,
                             DiffusionType &diffusion_operator,
                             std :: ofstream &data_file )
    : periodicDiscreteFunctionSpace_( periodicDiscreteFunctionSpace ),
      diffusion_( diffusion_operator ),
      data_file_( &data_file )
    {
    }

    //! ----------- method: solve cell problem ------------------------------------------

    template < class JacobianRangeImp >
    void solvecellproblem( JacobianRangeImp &gradient_PHI_H,
                           // the barycenter x_T of a macro grid element 'T'
                           const DomainType &globalQuadPoint,
                                 PeriodicDiscreteFunctionType &cell_problem_solution)
    {

      // set solution equal to zero:
      cell_problem_solution.clear();

      //! the matrix in our linear system of equations
      // in the non-linear case, it is the matrix for each iteration step
      CellFEMMatrix cell_system_matrix( "Cell Problem System Matrix", periodicDiscreteFunctionSpace_, periodicDiscreteFunctionSpace_ );

      //! define the discrete (elliptic) cell problem operator
      // ( effect of the discretized differential operator on a certain discrete function )
      CellProblemOperatorType cell_problem_op( periodicDiscreteFunctionSpace_, diffusion_);

      //! right hand side vector of the algebraic cell problem
      // (in the non-linear setting it changes for every iteration step)
      PeriodicDiscreteFunctionType cell_problem_rhs( "rhs of cell problem", periodicDiscreteFunctionSpace_ );
      cell_problem_rhs.clear();

      // NOTE:
      // is the right hand side of the cell problem equal to zero or almost identical to zero?
      // if yes, the solution of the cell problem is also identical to zero. The solver is getting a problem with this situation, which is why we do not solve cell problems for zero-right-hand-side, since we already know the result.

      // assemble the stiffness matrix
      cell_problem_op.assemble_matrix( globalQuadPoint, cell_system_matrix );

      // assemble right hand side of algebraic cell problem
      cell_problem_op.assembleCellRHS_linear( globalQuadPoint, gradient_PHI_H, cell_problem_rhs );

      const double norm_rhs = cell_problem_op.normRHS( cell_problem_rhs );

      if ( !( cell_problem_rhs.dofsValid() ) )
        { std :: cout << "Cell Problem RHS invalid." << std :: endl;
          abort(); }

      if ( norm_rhs < /*1e-06*/ 1e-10 )
        {
          cell_problem_solution.clear();
          //std :: cout << "Cell problem with solution zero." << std :: endl;
        }
      else
        {
          InverseCellFEMMatrix cell_fem_biCGStab( cell_system_matrix, 1e-8, 1e-8, 20000, CELLSOLVER_VERBOSE );
          cell_fem_biCGStab( cell_problem_rhs, cell_problem_solution );
        }

     if ( !(cell_problem_solution.dofsValid()) )
       {
         std::cout << "Current solution of the cell problem invalid!" << std::endl;
         std :: abort();
       }

    }

    //! ----------- end method: solve cell problem ------------------------------------------






    // method for solving and saving the solutions of the cell problems
    // for the whole set of macroscopic base function

    //! ---- method: solve and save the cell problems for the set of macroscopic base functions -----

    // here we need a 'cell problem numbering manager' to determine the number of the cell problem
    // (a combination of number of entity and number of local base function) 
    // Structure:
    // Struktur der Indizierung fuer das Abspeichern der Loesungen der Zellprobleme:
    // wir loesen Zellprobleme fuer jede Entity des Makro-Grids und jede Basisfunktion, die einen nichtleeren support auf dieser Entity besitzt, also schematisch:
    // Sei n=0,...,N einer Durchnumerierung der Entitys und i=0,...I_n eine zu einer festen Entity gehoerende Nummerierung der Basisfunktionen mit nicht-leeren support.
    // Die Durchnummerierung der Loesungen der Zellprobleme k=0,...,K ist dann gegeben durch: k(n,i_n) = ( sum_(l=0)^(n-1) ( I_l + 1) ) + i_n
    // NOTE: es verhaelt sich NICHT wie die vorhandene Methode mapToGlobal(entity,i) ! (die gibt die globale Nummer der Basisfunktion zurueck, es gibt aber  deutlich mehr Zellprobleme zum Loesen!
    // (das wird aber alles im Hintergrund vom 'cell problem numbering manager')

    // compute and save solutions of the cell problems for the base function set of the 'discreteFunctionSpace'
    // requires cell problem numbering manager
    template < class DiscreteFunctionImp, class CellProblemNumberingManagerImp >
    void saveTheSolutions_baseSet(
          const typename DiscreteFunctionImp::DiscreteFunctionSpaceType &discreteFunctionSpace,
          const CellProblemNumberingManagerImp &cp_num_manager, // just to check, if we use the correct numeration
          const std :: string &filename )
    {

      typedef DiscreteFunctionImp DiscreteFunctionType;

      typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType
         DiscreteFunctionSpaceType;

      typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;

      typedef typename DiscreteFunctionSpaceType :: GridType GridType;

      typedef typename DiscreteFunctionSpaceType :: JacobianRangeType
        JacobianRangeType;

      typedef typename DiscreteFunctionSpaceType :: BaseFunctionSetType
        BaseFunctionSetType;

      typedef typename DiscreteFunctionSpaceType :: IteratorType IteratorType;

      typedef typename DiscreteFunctionType :: LocalFunctionType
         LocalFunctionType;

      typedef typename GridType :: template Codim< 0 > :: Entity EntityType;
      typedef typename EntityType :: Geometry EntityGeometryType;

      typedef CachingQuadrature< GridPartType, 0 > EntityQuadratureType;

      enum { dimension = GridType :: dimension};
      enum { maxnumOfBaseFct = 100 }; 

      bool writer_is_open = false;

      std :: string cell_solution_location = "data/MsFEM/"+filename+"_cellSolutions_baseSet";
      DiscreteFunctionWriter dfw( (cell_solution_location).c_str() );

      writer_is_open = dfw.open();

      long double starting_time = clock();

      // we want to determine minimum, average and maxiumum time for solving a cell problem in the current method
      double minimum_time_c_p = 1000000;
      double average_time_c_p = 0;
      double maximum_time_c_p = 0;

      int number_of_cell_problem = 0;

      if ( writer_is_open )
      {

        IteratorType endit = discreteFunctionSpace.end();
        for(IteratorType it = discreteFunctionSpace.begin(); it != endit ; ++it)
          {

           // gradients of the macroscopic base functions
           JacobianRangeType gradientPhi[ maxnumOfBaseFct ];

           // entity
           const EntityType& entity = *it;

           const BaseFunctionSetType baseSet
              = discreteFunctionSpace.baseFunctionSet( entity );

           const EntityGeometryType &geometry = entity.geometry();

           EntityQuadratureType quadrature( entity, 0 );

           DomainType barycenter_of_entity = geometry.global( quadrature.point(0) );

           // number of base functions on entity
           const int numBaseFunctions = baseSet.numBaseFunctions();

           // calc Jacobian inverse before volume is evaluated
           const FieldMatrix< double, dimension, dimension > &inv
                    = geometry.jacobianInverseTransposed( quadrature.point( 0 /*=quadraturePoint*/ ) );

           PeriodicDiscreteFunctionType correctorPhi_i( "corrector Phi_i" , periodicDiscreteFunctionSpace_ );


           for( int i = 0; i < numBaseFunctions; ++i )
            {
              baseSet.jacobian( i, quadrature[0 /*=quadraturePoint*/], gradientPhi[ i ] );
              // multiply with transpose of jacobian inverse
              gradientPhi[ i ][ 0 ] = FMatrixHelp :: mult( inv, gradientPhi[ i ][ 0 ] );
            }

           for( int i = 0; i < numBaseFunctions; ++i )
            {

             correctorPhi_i.clear();

             // take time
             long double time_now = clock();

             solvecellproblem<JacobianRangeType>
                  ( gradientPhi[ i ], barycenter_of_entity, correctorPhi_i );

             // min/max time
             if ( (clock()-time_now)/CLOCKS_PER_SEC > maximum_time_c_p )
               { maximum_time_c_p = (clock()-time_now)/CLOCKS_PER_SEC; }
             if ( (clock()-time_now)/CLOCKS_PER_SEC < minimum_time_c_p )
               { minimum_time_c_p = (clock()-time_now)/CLOCKS_PER_SEC; }

             dfw.append( correctorPhi_i );

             // check if we use a correct numeration of the cell problems:
             if ( !(cp_num_manager.get_number_of_cell_problem( it, i ) == number_of_cell_problem) )
               { std :: cout << "Numeration of cell problems incorrect." << std :: endl;
                 std :: abort(); }

             number_of_cell_problem++;

            }

         } // end: for-loop: IteratorType it

      } // end: 'if ( writer_is_open )'


      if (data_file_)
      {
         if (data_file_->is_open())
         {
           (*data_file_) << std :: endl;
           (*data_file_) << "In method: saveTheSolutions_baseSet." << std :: endl << std :: endl;
           (*data_file_) << "Cell problems solved for " << discreteFunctionSpace.grid().size(0) << " leaf entities." << std :: endl;
           (*data_file_) << "Minimum time for solving a cell problem = " << minimum_time_c_p << "s." << std :: endl;
           (*data_file_) << "Maximum time for solving a cell problem = " << maximum_time_c_p << "s." << std :: endl;
           (*data_file_) << "Average time for solving a cell problem = " << ((clock()-starting_time)/CLOCKS_PER_SEC)/number_of_cell_problem << "s." << std :: endl;
           (*data_file_) << "Total time for computing and saving the cell problems = " << ((clock()-starting_time)/CLOCKS_PER_SEC) << "s," << std :: endl << std :: endl;
         }
      }

    }


 }; //end class




}

#endif // #ifndef DiscreteElliptic_HH
