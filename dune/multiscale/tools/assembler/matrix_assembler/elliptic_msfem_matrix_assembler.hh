#ifndef DiscreteEllipticMSFEMOperator_HH
#define DiscreteEllipticMSFEMOperator_HH

#include <dune/common/fmatrix.hh>

#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/multiscale/tools/solver/MsFEM/msfem_localproblems/localproblemsolver.hh>

#include <dune/fem/operator/2order/lagrangematrixsetup.hh>

namespace Dune
{

  template< class HostDiscreteFunctionImp, class SubGridImp >
  class SubGridList
  {
  public:
    //! ---------------- typedefs for the HostDiscreteFunctionSpace -----------------------

    typedef HostDiscreteFunctionImp HostDiscreteFunctionType;
    
    //! type of discrete function space
    typedef typename HostDiscreteFunctionType :: DiscreteFunctionSpaceType
      HostDiscreteFunctionSpaceType;

    //! type of (non-discrete )function space
    typedef typename HostDiscreteFunctionSpaceType :: FunctionSpaceType FunctionSpaceType;

    //! type of grid partition
    typedef typename HostDiscreteFunctionSpaceType :: GridPartType HostGridPartType;

    //! type of grid
    typedef typename HostDiscreteFunctionSpaceType :: GridType HostGridType;
    
    typedef typename HostDiscreteFunctionSpaceType :: IteratorType MaxLevelHostIteratorType;
    
    typedef typename MaxLevelHostIteratorType :: Entity HostEntityType;

    typedef typename HostEntityType :: EntityPointer HostEntityPointerType;
    
    typedef typename HostGridType :: template Codim< 0 > :: template Partition< All_Partition > :: LevelIterator HostgridLevelEntityIteratorType;
   
    typedef typename HostGridType :: Traits :: LevelIndexSet HostGridLevelIndexSet;

    typedef typename HostEntityType :: template Codim< 2 > :: EntityPointer HostNodePointer;
    
    typedef typename HostGridPartType :: IntersectionIteratorType HostIntersectionIterator;

    //! ---------------- typedefs for the SubgridDiscreteFunctionSpace -----------------------
    //  ( typedefs for the local grid and the corresponding local ('sub') )discrete space ) 

    //! type of grid
    typedef SubGridImp SubGridType; 

    //! type of grid part
    typedef LeafGridPart< SubGridType > SubGridPartType;
    
    template< typename EntityPointerCollectionType >
    void enrichment( HostEntityPointerType& hit,
                     const HostGridPartType& hostGridPart,
                     SubGridType& subGrid,
                     EntityPointerCollectionType& entities_sharing_same_node,
                     int& layer,
                     int &Schritt )
    {

      layer -= 1;
      // loop over the nodes of the enity
      for ( int i = 0; i < (*hit).template count<2>(); i += 1 )
	{

	   const HostNodePointer node = (*hit).template subEntity<2>(i);
	   int global_index_node = hostGridPart.indexSet().index( *node );

	   for( int j = 0; j < entities_sharing_same_node[global_index_node].size(); j += 1 )
	      {
                 Schritt += 1;
		 if ( !( subGrid.template contains <0>( *entities_sharing_same_node[ global_index_node ][ j ] ) ) )
		   { subGrid.insertPartial( *entities_sharing_same_node[ global_index_node ][ j ] ); }
		 
		 if ( layer > 0 )
		  {
		    enrichment( entities_sharing_same_node[ global_index_node ][ j ],
			        hostGridPart, subGrid, entities_sharing_same_node, layer, Schritt );
		    layer += 1;
		  }
	      }

	}

   }
    
    SubGridList( const HostDiscreteFunctionSpaceType& hostSpace, const int& computational_level, bool silent = true )
    : hostSpace_( hostSpace ),
      computational_level_( computational_level ),
      silent_( silent )
    {
      
      const HostGridPartType& hostGridPart = hostSpace_.gridPart();

      HostGridType& hostGrid = hostSpace_.gridPart().grid();

      int number_of_nodes = hostGrid.size( hostGrid.maxLevel(), 2 /*codim*/ );
      
      // -------- identify the entities that share a certain node -------

      std :: vector< std :: vector < HostEntityPointerType > > entities_sharing_same_node( number_of_nodes );
      
      for( MaxLevelHostIteratorType it = hostSpace_.begin(); it != hostSpace_.end(); ++it )
        {
	  int number_of_nodes_in_entity = (*it).template count<2>();
	  for ( int i = 0; i < number_of_nodes_in_entity; i += 1 )
	    {
	      const HostNodePointer node = (*it).template subEntity<2>(i);
	      int global_index_node = hostGridPart.indexSet().index( *node );
	      
	      entities_sharing_same_node[ global_index_node ].push_back( it );
	    }
        }

#if 0
      for ( int i = 0; i < number_of_nodes ; i+=1 )
        {
          std :: cout << "Knoten " << i << " wird von " << entities_sharing_same_node[i].size() << " Entities geteilt." << std :: endl;
        }
#endif

      // ---------------------------------------------------------------- 
     
      
      // maximum level defined in this grid. Levels are numbered 0 ... maxLevel with 0 the coarsest level.
      int maxLevel = hostGrid.maxLevel();
      int level_difference = maxLevel - computational_level_;

      // number of grid entities of a given codim on a given level in this process.
      int number_of_level_host_entities = hostGrid.size( computational_level_, 0 /*codim*/ );


      std :: cout << "number_of_level_host_entities = " << number_of_level_host_entities << std :: endl;

      
// number of layers per coarse grid entity T:  U(T) is created by enrichting T with n(T)-layers.
std :: vector < int > number_of_layers( number_of_level_host_entities );
for ( int i = 0; i < number_of_level_host_entities; i+=1 ) { number_of_layers[i] = 2; }

      subGrid = new SubGridType* [ number_of_level_host_entities ];
      
      //! ----------- create subgrids --------------------

      const HostGridLevelIndexSet& hostGridLevelIndexSet = hostGrid.levelIndexSet(computational_level_);

      //SubGridType* subGrid[ number_of_level_host_entities ];

      HostgridLevelEntityIteratorType level_iterator_end = hostGrid.template lend< 0 >( computational_level_ );
      HostgridLevelEntityIteratorType level_iterator_begin = hostGrid.template lbegin< 0 >( computational_level_ );

      int lev_index = 0;
      for ( HostgridLevelEntityIteratorType lit = level_iterator_begin;
         lit != level_iterator_end ; ++lit )
        {

          subGrid[ lev_index ] = new SubGridType( hostGrid );
          subGrid[ lev_index ]->createBegin();

          //level_entity_collection.push_back( lit );
          lev_index += 1;

        }

      for( MaxLevelHostIteratorType it = hostSpace_.begin(); it != hostSpace_.end(); ++it )
        {
	  int number_of_nodes_in_entity = (*it).template count<2>();
	  for ( int i = 0; i < number_of_nodes_in_entity; i += 1 )
	    {
	      const HostNodePointer node = (*it).template subEntity<2>(i);
	      int global_index_node = hostGridPart.indexSet().index( *node );
	      
	      entities_sharing_same_node[ global_index_node ].push_back( it );
	    }
        }





int Schritt = 0;



      // a maxlevel iterator for the codim 0 hostgrid entities:
      MaxLevelHostIteratorType host_endit = hostSpace_.end();
      for( MaxLevelHostIteratorType host_it = hostSpace_.begin();
           host_it != host_endit ;
           ++host_it )
        {

           const HostEntityType& host_entity = *host_it;
	   
           int number_of_nodes_in_entity = (*host_it).template count<2>();
	   

           // get the level 'computational_level'-father of host_entity (which is a maxlevel entity)
           HostEntityPointerType level_father_entity = host_it;
           for (int lev = 0; lev < level_difference; ++lev)
             level_father_entity = level_father_entity->father();

           int father_index = hostGridLevelIndexSet.index( *level_father_entity );
           // std :: cout << "father_index = " << father_index << std :: endl;

           Schritt += 1;
           if ( !( subGrid[ father_index ]->template contains <0>(host_entity) ) )
            { subGrid[ father_index ]->insertPartial( host_entity ); }


           // check the neighbor entities and look if they belong to the same father
           // if yes, continue
           // if not, enrichement with 'n(T)'-layers
           bool all_neighbors_have_same_father = true;
           const HostIntersectionIterator iend = hostGridPart.iend( host_entity );
           for( HostIntersectionIterator iit = hostGridPart.ibegin( host_entity ); iit != iend; ++iit )
             {

                if ( iit->neighbor() ) //if there is a neighbor entity
                  {
                    // check if the neighbor entity is in the subgrid
                   const HostEntityPointerType neighborHostEntityPointer = iit->outside();
                   const HostEntityType& neighborHostEntity = *neighborHostEntityPointer;
		   
                   HostEntityPointerType level_father_neighbor_entity = neighborHostEntityPointer;
                   for (int lev = 0; lev < level_difference; ++lev)
                       level_father_neighbor_entity = level_father_neighbor_entity->father();
	   
	   
                   if ( !(level_father_neighbor_entity == level_father_entity) )
                    {
                      all_neighbors_have_same_father = false;
                    }

                  }
                else
		  {
		    all_neighbors_have_same_father = false;
		  }
		
              }
           if ( all_neighbors_have_same_father == true )
	      { continue; }

	   int layers = number_of_layers[ father_index ];
	   enrichment( host_it, hostGridPart, *subGrid[ father_index ], entities_sharing_same_node, layers, Schritt );

#if 0
           // enrichment:
           //! TODO ENRICHMENT for layers > 1!
           HostEntityPointerType hit = host_it;
           for ( int layer = 0; layer < number_of_layers[ father_index ]; layer += 1 )
	    {
	      for ( int i = 0; i < (*hit).template count<2>(); i += 1 )
	        {

	          const HostNodePointer node = (*hit).template subEntity<2>(i);
	          int global_index_node = hostGridPart.indexSet().index( *node );

	          for( int j = 0; j < entities_sharing_same_node[global_index_node].size(); j += 1 )
	            {
		      if ( !( subGrid[ father_index ]->template contains <0>( *entities_sharing_same_node[ global_index_node ][ j ] ) ) )
		       { subGrid[ father_index ]->insertPartial( *entities_sharing_same_node[ global_index_node ][ j ] ); }
		    }
		
	        }
	     }
#endif
            
        }


std :: cout << "Schritt = " << Schritt << std :: endl; abort();

      for ( int i = 0; i < number_of_level_host_entities; ++i )
       {
          subGrid[ i ]->createEnd();
          if ( !silent_ )
            { subGrid[ i ]->report(); }
       }

      //! ----------- end create subgrids --------------------

    }
    
  SubGridType& getSubGrid( int i )
  {
    int size = hostSpace_.gridPart().grid().size( computational_level_, 0 /*codim*/ );
    
    if ( i < size )    
     { return *(subGrid[ i ]); }
    else
     { std :: cout << "Error. Subgrid-Index too large." << std :: endl; }
  }
    
  private:
    
    const HostDiscreteFunctionSpaceType &hostSpace_;
    const int computational_level_;
    
    bool silent_;
    
    SubGridType** subGrid;
    
  public:
    
    //Destruktor (Dynamisch angeforderter Speicher wird wieder freigegeben)
    ~SubGridList()
     {
        int size = hostSpace_.gridPart().grid().size( computational_level_, 0 /*codim*/ );
        
        for (unsigned int i = 0;i<size; ++i)
         delete subGrid[i];
        delete[] subGrid;

     }

  };

  
  // Imp stands for Implementation
  template< class CoarseDiscreteFunctionImp, class FineDiscreteFunctionImp, class DiffusionImp >
  class DiscreteEllipticMsFEMOperator
  : public Operator< typename CoarseDiscreteFunctionImp::RangeFieldType,
                     typename CoarseDiscreteFunctionImp::RangeFieldType,
		      CoarseDiscreteFunctionImp, CoarseDiscreteFunctionImp >
  {
    typedef DiscreteEllipticMsFEMOperator< CoarseDiscreteFunctionImp, FineDiscreteFunctionImp, DiffusionImp > This;

  public:
    
    typedef CoarseDiscreteFunctionImp CoarseDiscreteFunction;
    typedef FineDiscreteFunctionImp FineDiscreteFunction;
   
    typedef DiffusionImp DiffusionModel;
    

    typedef typename CoarseDiscreteFunction::DiscreteFunctionSpaceType CoarseDiscreteFunctionSpace;    
    typedef typename FineDiscreteFunction::DiscreteFunctionSpaceType FineDiscreteFunctionSpace;

    typedef typename FineDiscreteFunctionSpace::GridPartType FineGridPart;
    typedef typename FineDiscreteFunctionSpace::GridType FineGrid;

    typedef typename FineDiscreteFunctionSpace::RangeFieldType RangeFieldType;
    typedef typename FineDiscreteFunctionSpace::DomainType DomainType;
    typedef typename FineDiscreteFunctionSpace::RangeType RangeType;
    typedef typename FineDiscreteFunctionSpace::JacobianRangeType
      JacobianRangeType;
      
    typedef SubGrid< WORLDDIM , FineGrid > SubGridType;
    typedef SubGridList< FineDiscreteFunction, SubGridType > SubGridListType;
    
    typedef MsFEMLocalProblemSolver< FineDiscreteFunction, SubGridListType, DiffusionModel > MsFEMLocalProblemSolverType;

  protected:
    static const int dimension = FineGridPart::GridType::dimension;
    static const int polynomialOrder = FineDiscreteFunctionSpace::polynomialOrder;

    typedef typename FineDiscreteFunction::LocalFunctionType FineLocalFunction;

    typedef typename FineDiscreteFunctionSpace::BaseFunctionSetType FineBaseFunctionSet;
    typedef typename FineDiscreteFunctionSpace::LagrangePointSetType FineLagrangePointSet;
    typedef typename FineLagrangePointSet::template Codim< 1 >::SubEntityIteratorType FineFaceDofIterator;

    typedef typename FineDiscreteFunctionSpace::IteratorType FineIterator;
    typedef typename FineIterator::Entity FineEntity;
    typedef typename FineEntity::Geometry FineGeometry;

    typedef typename FineGridPart::IntersectionIteratorType FineIntersectionIterator;
    typedef typename FineIntersectionIterator::Intersection FineIntersection;

    typedef typename FineGrid :: template Codim< 0 > :: template  Partition< All_Partition > :: LevelIterator FineLevelEntityIterator;

    typedef CachingQuadrature< FineGridPart, 0 > FineQuadrature;


  public:

    typedef typename CoarseDiscreteFunctionSpace::GridPartType CoarseGridPart;
    typedef typename CoarseDiscreteFunctionSpace::GridType CoarseGrid;

  protected:

    typedef typename CoarseDiscreteFunction::LocalFunctionType CoarseLocalFunction;

    typedef typename CoarseDiscreteFunctionSpace::BaseFunctionSetType CoarseBaseFunctionSet;
    typedef typename CoarseDiscreteFunctionSpace::LagrangePointSetType CoarseLagrangePointSet;
    typedef typename CoarseLagrangePointSet::template Codim< 1 >::SubEntityIteratorType CoarseFaceDofIterator;

    typedef typename CoarseDiscreteFunctionSpace::IteratorType CoarseIterator;
    typedef typename CoarseIterator::Entity CoarseEntity;
    typedef typename CoarseEntity::Geometry CoarseGeometry;

    typedef typename CoarseGridPart::IntersectionIteratorType CoarseIntersectionIterator;
    typedef typename CoarseIntersectionIterator::Intersection CoarseIntersection;

    typedef typename CoarseGrid :: template Codim< 0 > :: template  Partition< All_Partition > :: LevelIterator CoarseLevelEntityIterator;

    typedef CachingQuadrature< CoarseGridPart, 0 > CoarseQuadrature;
    
  public:
    
    DiscreteEllipticMsFEMOperator( const CoarseDiscreteFunctionSpace &coarseDiscreteFunctionSpace,
                                   const FineDiscreteFunctionSpace &fineDiscreteFunctionSpace,
                                   const DiffusionModel &diffusion_op,
                                   std :: ofstream& data_file,
                                   std :: string path = ""  )
    : coarseDiscreteFunctionSpace_( coarseDiscreteFunctionSpace ),
      fineDiscreteFunctionSpace_( fineDiscreteFunctionSpace ),
      diffusion_operator_( diffusion_op ),
      data_file_( &data_file ),
      path_( path )
    {
      bool silence = false;
      
      const int coarse_level = coarseDiscreteFunctionSpace_.gridPart().grid().maxLevel();
      
      SubGridListType subgrid_list( fineDiscreteFunctionSpace_ , coarse_level , silence );
      
      MsFEMLocalProblemSolverType loc_prob_solver( fineDiscreteFunctionSpace_, subgrid_list, diffusion_operator_, data_file, path_ );
      loc_prob_solver.assemble_all( coarse_level, silence );
    }


  private:
    DiscreteEllipticMsFEMOperator ( const This & );

  public:

    // dummy operator
    virtual void
    operator() ( const CoarseDiscreteFunction &u, CoarseDiscreteFunction &w ) const;


    template< class MatrixType >
    void assemble_matrix ( MatrixType &global_matrix ) const;

  private:
    
    const FineDiscreteFunctionSpace &fineDiscreteFunctionSpace_;
    const CoarseDiscreteFunctionSpace &coarseDiscreteFunctionSpace_;
    const DiffusionModel &diffusion_operator_;

    // data file for saving information
    std :: ofstream *data_file_;

    // path where to save the data output
    std :: string path_;

  };



  // dummy implementation of "operator()"
  // 'w' = effect of the discrete operator on 'u'
  template< class CoarseDiscreteFunctionImp, class FineDiscreteFunctionImp, class DiffusionImp >
  void DiscreteEllipticMsFEMOperator< CoarseDiscreteFunctionImp, 
                                      FineDiscreteFunctionImp, 
				       DiffusionImp > :: operator() ( const CoarseDiscreteFunction &u, CoarseDiscreteFunction &w ) const 
  {

    std :: cout << "the ()-operator of the DiscreteEllipticMsFEMOperator class is not yet implemented and still a dummy." << std :: endl;
    std :: abort();

  }


  template< class CoarseDiscreteFunctionImp, class FineDiscreteFunctionImp, class DiffusionImp >
  template< class MatrixType >
  void DiscreteEllipticMsFEMOperator< CoarseDiscreteFunctionImp, 
                                      FineDiscreteFunctionImp, 
				       DiffusionImp > :: assemble_matrix ( MatrixType &global_matrix ) const
  {
#if 0
    
    //! der braucht einen macro-space (der wird auch als subgrid-space generiert) und den total globalen Feinskalen-Raum
    
   
    // the local problem:
    // Let 'T' denote a coarse grid element and
    // let 'U(T)' denote the environment of 'T' that corresponds with the subgrid.
    
    // if Petrov-Galerkin-MsFEM
    #ifdef PGF
    std :: cout << "Assembling Petrov-Galerkin-MsFEM Matrix." << std :: endl;
    #else
    std :: cout << "Assembling MsFEM Matrix." << std :: endl;
    #endif

    std :: string local_solution_location;

    // the file/place, where we saved the solutions of the cell problems
    local_solution_location = "data/MsFEM/"+(*filename_)+"/local_problems/_localProblemSolutions_baseSet";


    bool reader_is_open = false;
    // reader for the cell problem data file:
    DiscreteFunctionReader discrete_function_reader( (local_solution_location).c_str() );
    reader_is_open = discrete_function_reader.open();

    typedef typename MatrixType::LocalMatrixType LocalMatrix;

    global_matrix.reserve();
    global_matrix.clear();

    std::vector< typename BaseFunctionSet::JacobianRangeType > gradient_Phi( discreteFunctionSpace_.mapper().maxNumDofs() );
    
//! assemble local grids, load local solutions
    
    //! LevelEntityIterator
    
    const Iterator macro_grid_end = discreteFunctionSpace_.end();
    for( Iterator macro_grid_it = discreteFunctionSpace_.begin(); macro_grid_it != macro_grid_end; ++macro_grid_it )
    {

      // the coarse grid element T:
      const Entity &macro_grid_entity = *macro_grid_it;
      const Geometry &macro_grid_geometry = macro_grid_entity.geometry();
      assert( macro_grid_entity.partitionType() == InteriorEntity );

      LocalMatrix local_matrix = global_matrix.localMatrix( macro_grid_entity, macro_grid_entity );

      const BaseFunctionSet &macro_grid_baseSet = local_matrix.domainBaseFunctionSet();
      const unsigned int numMacroBaseFunctions = macro_grid_baseSet.numBaseFunctions();

#if 0
      
      // 1 point quadrature!! That is how we compute and save the cell problems.
      // If you want to use a higher order quadrature, you also need to change the computation of the cell problems!
      Quadrature one_point_quadrature( macro_grid_entity, 0 );

      // the barycenter of the macro_grid_entity
      const typename Quadrature::CoordinateType &local_macro_point = one_point_quadrature.point( 0 /*=quadraturePoint*/ );
      DomainType macro_entity_barycenter = macro_grid_geometry.global( local_macro_point );

      const double macro_entity_volume = one_point_quadrature.weight( 0 /*=quadraturePoint*/ ) * macro_grid_geometry.integrationElement( local_macro_point );

      // transposed of the the inverse jacobian
      const FieldMatrix< double, dimension, dimension > &inverse_jac
          = macro_grid_geometry.jacobianInverseTransposed( local_macro_point );

      int cell_problem_id [ numMacroBaseFunctions ];

      DiscreteFunction* corrector_Phi[discreteFunctionSpace_.mapper().maxNumDofs()];

      for( unsigned int i = 0; i < numMacroBaseFunctions; ++i )
        {

          // get number of cell problem from entity and number of base function
          cell_problem_id[i] = lp_num_manager_.get_number_of_local_problem( macro_grid_it, i );

          // jacobian of the base functions, with respect to the reference element
          typename BaseFunctionSet::JacobianRangeType gradient_Phi_ref_element;
          macro_grid_baseSet.jacobian( i, one_point_quadrature[ 0 ], gradient_Phi_ref_element );

          // multiply it with transpose of jacobian inverse to obtain the jacobian with respect to the real entity
          inverse_jac.mv( gradient_Phi_ref_element[ 0 ], gradient_Phi[ i ][ 0 ] );

          // ( Q^eps(\Phi_i) ○ F ):
          corrector_Phi[i] = new DiscreteFunction( "Corrector Function of Phi", localDiscreteFunctionSpace_ );
          corrector_Phi[i]->clear();
          #ifdef AD_HOC_COMPUTATION
          MsFEMLocalProblemSolverType cell_problem_solver( localDiscreteFunctionSpace_, diffusion_operator_ );
          cell_problem_solver.template solvelocalproblem<JacobianRangeType>
                ( gradient_Phi[ i ], macro_entity_barycenter, *(corrector_Phi[ i ]) );
          #else
          if (reader_is_open)
            { discrete_function_reader.read( cell_problem_id[ i ], *(corrector_Phi[ i ]) ); }
          #endif

        }

      for( unsigned int i = 0; i < numMacroBaseFunctions; ++i )
        {

          for( unsigned int j = 0; j < numMacroBaseFunctions; ++j )
           {

            RangeType local_integral = 0.0;


#if 1
DomainType A_0_MSFEM(0.0);
#endif

            // iterator for the micro grid ( grid for the reference element T_0 )
            const Iterator micro_grid_end = localDiscreteFunctionSpace_.end();
            for( Iterator micro_grid_it = localDiscreteFunctionSpace_.begin(); micro_grid_it != micro_grid_end; ++micro_grid_it )
              {

                // remember:
                // |det(A)| \int_{T_0} (A^eps ○ F)(x) ( ∇\Phi_i(x_T) + (A^{-1})^T ∇( Q^eps(\Phi_i) ○ F )(x)) · ( ∇\Phi_j(x_T) + (A^{-1})^T ∇( Q^eps(\Phi_j) ○ F )(x))

#if 1
gradient_Phi[ i ][ 0 ][ 0 ] = 1.0;
gradient_Phi[ i ][ 0 ][ 1 ] = 0.0;
#endif

                const Entity &micro_grid_entity = *micro_grid_it;
                const Geometry &micro_grid_geometry = micro_grid_entity.geometry();
                assert( micro_grid_entity.partitionType() == InteriorEntity );

                // ( Q^eps(\Phi_i) ○ F ):
                typename DiscreteFunction::LocalFunctionType localized_corrector_i = corrector_Phi[ i ]->localFunction( micro_grid_entity );
                // ( Q^eps(\Phi_j) ○ F ):
                typename DiscreteFunction::LocalFunctionType localized_corrector_j = corrector_Phi[ j ]->localFunction( micro_grid_entity );

                // higher order quadrature, since A^{\epsilon} is highly variable
                Quadrature micro_grid_quadrature( micro_grid_entity, 2*localDiscreteFunctionSpace_.order()+2 );
                const size_t numQuadraturePoints = micro_grid_quadrature.nop();

                for( size_t microQuadraturePoint = 0; microQuadraturePoint < numQuadraturePoints; ++microQuadraturePoint )
                  {

                    // local (barycentric) coordinates (with respect to entity)
                    const typename Quadrature::CoordinateType &local_micro_point = micro_grid_quadrature.point( microQuadraturePoint );

                    DomainType global_point_in_T_0 = micro_grid_geometry.global( local_micro_point );

                    const double weight_micro_quadrature = micro_grid_quadrature.weight(  microQuadraturePoint ) * micro_grid_geometry.integrationElement( local_micro_point );

                    // ∇( Q^eps(\Phi_i) ○ F )   and   ∇( Q^eps(\Phi_j) ○ F )
                    JacobianRangeType grad_corrector_i, grad_corrector_j;
                    localized_corrector_i.jacobian( micro_grid_quadrature[ microQuadraturePoint ], grad_corrector_i );
                    localized_corrector_j.jacobian( micro_grid_quadrature[ microQuadraturePoint ], grad_corrector_j );

#if 1
                    // global point in the reference element T_0
                    DomainType global_point = micro_grid_geometry.global( local_micro_point );

                    // 'F(x)', i.e. F ( global point in the reference element T_0 )
                    // (the transformation of the global point in T_0 to its position in T)
                    DomainType global_point_transformed(0.0);

                    for( int k = 0; k < dimension; ++k )
                      for( int l = 0; l < dimension; ++l )
                        global_point_transformed[ k ] += ( val_A[ k ][ l ] * global_point_in_T_0[ l ] );

                    global_point_transformed += corner_0_of_T;

                    // F(x) = Ax + a_0, F : T_0 -> T is given by
                    // A_11 = a_1(1) - a_0(1)     A_12 = a_2(1) - a_0(1)
                    // A_21 = a_1(2) - a_0(2)     A_22 = a_2(2) - a_0(2)

#endif
                    // direction of the diffusion: ( ∇\Phi_i(x_T) + (A^{-1})^T ∇( Q^eps(\Phi_i) ○ F )(x))
                    JacobianRangeType direction_of_diffusion( 0.0 );
                    for( int k = 0; k < dimension; ++k )
                      {
                        for( int l = 0; l < dimension; ++l )
                          {
                           direction_of_diffusion[ 0 ][ k ] += val_A_inverse_transposed[ k ][ l ] * grad_corrector_i[ 0 ][ l ];
                          }

                        direction_of_diffusion[ 0 ][ k ] += gradient_Phi[ i ][ 0 ][ k ];

                      }

                    JacobianRangeType diffusion_in_gradient_Phi_reconstructed;
                    diffusion_operator_.diffusiveFlux( global_point_transformed,
                                                       direction_of_diffusion, diffusion_in_gradient_Phi_reconstructed );


                    // if test function reconstruction
                    #ifndef PGF


                    // ( ∇\Phi_j(x_T) + (A^{-1})^T ∇( Q^eps(\Phi_j) ○ F )(x)):
                    JacobianRangeType grad_reconstruction_Phi_j( 0.0 );
                    for( int k = 0; k < dimension; ++k )
                      {
                        for( int l = 0; l < dimension; ++l )
                          { grad_reconstruction_Phi_j[ 0 ][ k ] += val_A_inverse_transposed[ k ][ l ] * grad_corrector_j[ 0 ][ l ]; }
                        grad_reconstruction_Phi_j[ 0 ][ k ] += gradient_Phi[ j ][ 0 ][ k ];
                      }

                    local_integral += weight_micro_quadrature * ( diffusion_in_gradient_Phi_reconstructed[ 0 ] * grad_reconstruction_Phi_j[ 0 ]);
                    #else
                    local_integral += weight_micro_quadrature * ( diffusion_in_gradient_Phi_reconstructed[ 0 ] * gradient_Phi[ j ][ 0 ]);

#if 1
A_0_MSFEM[ 0 ] +=  2.0 * weight_micro_quadrature * diffusion_in_gradient_Phi_reconstructed[ 0 ][ 0 ];
A_0_MSFEM[ 1 ] +=  2.0 * weight_micro_quadrature * diffusion_in_gradient_Phi_reconstructed[ 0 ][ 1 ];
#endif

                    #endif

                  }
              }

#if 1
std :: cout << "A_0_MSFEM[ 0 ] = " << A_0_MSFEM[ 0 ] << std :: endl;
std :: cout << "A_0_MSFEM[ 1 ] = " << A_0_MSFEM[ 1 ] << std :: endl;
#endif

            // add |det(A)|*\int_{T_0} ...
            local_matrix.add( j, i, abs_det_A * local_integral );
           }

        }

      // delete?
      //delete[] corrector_Phi;
#endif
    }

    //discrete_function_reader.close();

    // boundary treatment
    const GridPart &gridPart = discreteFunctionSpace_.gridPart();
    for( Iterator it = discreteFunctionSpace_.begin(); it != macro_grid_end; ++it )
    {
      const Entity &entity = *it;
      if( !entity.hasBoundaryIntersections() )
        continue;

      LocalMatrix local_matrix = global_matrix.localMatrix( entity, entity );

      const LagrangePointSet &lagrangePointSet = discreteFunctionSpace_.lagrangePointSet( entity );

      const IntersectionIterator iend = gridPart.iend( entity );
      for( IntersectionIterator iit = gridPart.ibegin( entity ); iit != iend; ++iit )
      {
        const Intersection &intersection = *iit;
        if( !intersection.boundary() )
          continue;

        const int face = intersection.indexInInside();
        const FaceDofIterator fdend = lagrangePointSet.template endSubEntity< 1 >( face );
        for( FaceDofIterator fdit = lagrangePointSet.template beginSubEntity< 1 >( face ); fdit != fdend; ++fdit )
          local_matrix.unitRow( *fdit );
      }
    }
#endif
  }


//! ------------------------------------------------------------------------------------------------
//! ------------------------------------------------------------------------------------------------


}

#endif // #ifndef DiscreteElliptic_HH
