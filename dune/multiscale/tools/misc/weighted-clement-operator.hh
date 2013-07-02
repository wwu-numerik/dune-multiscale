// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_DIVERGENCE_HH
#define DUNE_DIVERGENCE_HH

#ifdef HAVE_CMAKE_CONFIG
  #include "cmake_config.h"
#else
  #include "config.h"
#endif // ifdef HAVE_CMAKE_CONFIG

#include <unordered_map>
#include <unordered_set>

#include <dune/common/fmatrix.hh>
#include <dune/common/timer.hh>
#include <dune/fem/storage/array.hh>
#include <dune/fem/quadrature/quadrature.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/operator/matrix/spmatrix.hh>
#include <dune/fem/solver/oemsolver.hh>
#include <dune/fem/operator/2order/lagrangematrixsetup.hh>
#include <dune/stuff/fem/localmatrix_proxy.hh>
#include <dune/stuff/la/container/pattern.hh>

#include <dune/multiscale/common/traits.hh>
#include <dune/multiscale/msfem/msfem_traits.hh>
#include <dune/multiscale/msfem/msfem_grid_specifier.hh>
#include <dune/multiscale/msfem/localproblems/subgrid-list.hh>
#include <dune/multiscale/tools/misc.hh>
#include <dune/multiscale/tools/misc/linear-lagrange-interpolation.hh>


namespace Dune {
namespace Multiscale {
namespace MsFEM {

template < class FunctionSpaceTraits >
std::vector<int> mapEach(const Dune::Fem::DiscreteFunctionSpaceInterface<FunctionSpaceTraits>& space,
                         const typename Dune::Fem::DiscreteFunctionSpaceInterface<FunctionSpaceTraits>::EntityType& entity)
{
  const auto& mapper = space.blockMapper();
  std::vector<int> ret(mapper.numDofs(entity));
  auto add = [&](int localDof, int globalDof){ ret[localDof] = globalDof; };
  mapper.mapEach(entity, add);
  return ret;
}

template <class IndexSetType>
struct EntityPointerHash {
  const IndexSetType& index_set_;

  EntityPointerHash(const IndexSetType& index_set)
    :index_set_(index_set)
  {}

  template < class GridImp, class IteratorImp >
  std::size_t operator() (const Dune::EntityPointer< GridImp, IteratorImp >& ptr) const
  {
    return index_set_.index(*ptr);
  }
};

template<class DomainSpace, class RangeSpace>
class ClemementPattern : public DSL::SparsityPatternDefault {
  typedef DSL::SparsityPatternDefault BaseType;
  typedef typename DomainSpace::EntityType::EntityPointer DomainEntityPointerType;
  typedef typename RangeSpace::EntityType::EntityPointer RangeEntityPointerType;
  typedef std::unordered_set<RangeEntityPointerType,
                             EntityPointerHash<typename RangeSpace::IndexSetType> >
    FineEntitySetType;
  typedef std::unordered_map<DomainEntityPointerType,
                             FineEntitySetType,
                             EntityPointerHash<typename DomainSpace::IndexSetType> >
    SupportMapType;
  SupportMapType support_map_;
public:
  //! creates an entity/neighbor pattern with domainSpace.size() == #rows sets
  ClemementPattern(const DomainSpace& domainSpace,
                   const RangeSpace& rangeSpace,
                   const MacroMicroGridSpecifier& specifier)
    : BaseType(domainSpace.size())
    , support_map_(domainSpace.gridPart().grid().size(0),
                   typename SupportMapType::hasher(domainSpace.indexSet()))
  {

    for (const auto& domain_entity : domainSpace) {
      const auto globalI_vec = mapEach(domainSpace, domain_entity);
      FineEntitySetType range_set(specifier.getLevelDifference()*3,
                                 typename FineEntitySetType::hasher(rangeSpace.indexSet()));
      const auto father_of_loc_grid_ent =
        Stuff::Grid::make_father(rangeSpace.gridPart().grid().leafIndexSet(),
                                 domainSpace.grid().template getHostEntity< 0 >(domain_entity),
                                 specifier.getLevelDifference());
      for(const auto& range_entity : rangeSpace)
      {
        if (!Stuff::Grid::entities_identical(range_entity, *father_of_loc_grid_ent))
          continue;
        range_set.insert(RangeEntityPointerType(range_entity));
        for (const auto i : globalI_vec) {
          const auto globalJ_vec = mapEach(rangeSpace,range_entity);
          auto& columns = BaseType::inner(i);
          for (const auto j : globalJ_vec) {
            columns.insert(j);
          }
        }
      }
      support_map_.insert(std::make_pair(DomainEntityPointerType(domain_entity), range_set));
    }
  }

  const SupportMapType& support() const {
    return support_map_;
  }
};

class WeightedClementOperator
: public Operator< typename SubGridList::SubGridDiscreteFunction::RangeFieldType,
                   typename CommonTraits::DiscreteFunctionType::RangeFieldType,
                   typename SubGridList::SubGridDiscreteFunction,
                   typename CommonTraits::DiscreteFunctionType>,
  public OEMSolver::PreconditionInterface
{
private:
  typedef std::vector<std::shared_ptr<CommonTraits::DiscreteFunctionType>> CoarseBasisFunctionList;
  typedef CommonTraits::DiscreteFunctionType CoarseDiscreteFunction;

  //! type of discrete functions
  typedef SubGridList::SubGridDiscreteFunction DiscreteFunctionType;

  typedef CoarseDiscreteFunction CoarseDiscreteFunctionType;

  //! type of this LaplaceFEOp
  typedef WeightedClementOperator WeightedClementOpType;

  //! needs to be friend for conversion check
  friend class Conversion<WeightedClementOpType,OEMSolver::PreconditionInterface>;

  typedef WeightedClementOpType ThisType;

  //! type of discrete function space
  typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType
    DiscreteFunctionSpaceType;

  typedef typename CoarseDiscreteFunction :: DiscreteFunctionSpaceType
    CoarseDiscreteFunctionSpaceType;

  //! field type of range
  typedef typename DiscreteFunctionSpaceType :: RangeType
    RangeType;

  typedef typename DiscreteFunctionSpaceType :: DomainType
    DomainType;

  //! field type of range
  typedef typename DiscreteFunctionSpaceType :: RangeFieldType
    RangeFieldType;

  //! type of grid partition
  typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;
  //! type of grid
  typedef typename DiscreteFunctionSpaceType :: GridType GridType;

  typedef typename CoarseDiscreteFunctionSpaceType :: GridPartType CoarseGridPartType;

  //! type of jacobian
  typedef typename DiscreteFunctionSpaceType :: JacobianRangeType
    JacobianRangeType;
  //! type of the base function set
  typedef typename DiscreteFunctionSpaceType :: BasisFunctionSetType
    BasisFunctionSetType;

  typedef typename CoarseDiscreteFunctionSpaceType :: BasisFunctionSetType
    CoarseBasisFunctionSetType;

  enum{dimRange = GridType::dimension};

  //! polynomial order of base functions
  enum { polynomialOrder = DiscreteFunctionSpaceType :: polynomialOrder };

  //! The grid's dimension
  enum { dimension = GridType :: dimension };

  //! type of quadrature to be used
  typedef Fem::CachingQuadrature< GridPartType, 0 > QuadratureType;
  typedef Fem::CachingQuadrature< CoarseGridPartType, 0 > CoarseQuadratureType;

  typedef Dune::Fem::SparseRowMatrixTraits < typename SubGridList::SubGridDiscreteFunctionSpace,
                                       typename SubGridList::HostDiscreteFunctionSpaceType > WeightedClementMatrixObjectTraits;

  typedef typename WeightedClementMatrixObjectTraits
    :: MatrixObject< Dune::LagrangeMatrixTraits< WeightedClementMatrixObjectTraits > >
    :: MatrixObjectType LinearOperatorType;

  //! get important types from the MatrixObject
  typedef typename LinearOperatorType :: LocalMatrixType LocalMatrixType;
  typedef typename LinearOperatorType :: PreconditionMatrixType PreconditionMatrixType;
  typedef typename LinearOperatorType :: MatrixType MatrixType;
  typedef typename DiscreteFunctionSpaceType :: MapperType MapperType;
  typedef std::vector<DomainType> CoarseNodeVectorType;

  // type of DofManager
  typedef Fem::DofManager< GridType > DofManagerType;


public:
  //! constructor
  explicit WeightedClementOperator( const DiscreteFunctionSpaceType& space,
                              const CoarseDiscreteFunctionSpaceType& coarse_space,
                              const CoarseNodeVectorType& coarse_nodes,
                              const CoarseBasisFunctionList& coarse_basis,
                              const std::map<int,int>& global_id_to_internal_id,
                              const MacroMicroGridSpecifier& specifier )
  : discreteFunctionSpace_( space ),
    coarse_space_( coarse_space ),
    coarse_nodes_( coarse_nodes ),
    coarse_basis_( coarse_basis ),
    global_id_to_internal_id_( global_id_to_internal_id ),
    dofManager_( DofManagerType :: instance( space.grid() ) ),
    specifier_( specifier )
    , sparsity_pattern_(discreteFunctionSpace_, coarse_space_, specifier_)
    , linearOperator_( discreteFunctionSpace_, coarse_space_ ),
    sequence_( -1 ),
    gradCache_( discreteFunctionSpace_.mapper().maxNumDofs() ),
    values_( discreteFunctionSpace_.mapper().maxNumDofs() )
  {
  }

private:
  // prohibit copying
  WeightedClementOperator ( const ThisType & ) = delete;

public:                                                           /*@LST0S@*/
  //! apply the operator
  virtual void operator() ( const DiscreteFunctionType &u,
                            CoarseDiscreteFunctionType &w ) const
  {
    systemMatrix().apply( u, w );                                 /*@\label{sto:matrixEval}@*/
  }                                                               /*@LST0E@*/

  //! return reference to preconditioning matrix, used by OEM-Solver
  const PreconditionMatrixType &preconditionMatrix () const
  {
    return systemMatrix().preconditionMatrix();
  }
                                                                  /*@LST0S@*/
  virtual void applyTransposed ( const CoarseDiscreteFunctionType &u,
                                 DiscreteFunctionType &w) const
  {
    systemMatrix().apply_t(u,w);                /*@\label{sto:applytransposed}@*/
  }                                                             /*@LST0E@*/

  //! return true if preconditioning is enabled
  bool hasPreconditionMatrix () const
  {
    return linearOperator_.hasPreconditionMatrix();
  }

  //! print the system matrix into a stream
  void print ( std :: ostream & out = std :: cout ) const
  {
    systemMatrix().matrix().print( out );
  }

  //! return reference to discreteFunctionSpace
  const DiscreteFunctionSpaceType &discreteFunctionSpace () const
  {
    return discreteFunctionSpace_;
  }

  /*! \brief obtain a reference to the system matrix
   *
   *  The assembled matrix is returned. If the system matrix has not been
   *  assembled, yet, the assembly is performed.
   *
   *  \returns a reference to the system matrix
   */
  LinearOperatorType &systemMatrix () const                        /*@LST0S@*/
  {
    // if stored sequence number it not equal to the one of the
    // dofManager (or space) then the grid has been changed
    // and matrix has to be assembled new
    if( sequence_ != dofManager_.sequence() )                     /*@\label{sto:sequence}@*/
      assemble();

    return linearOperator_;
  }

  /** \brief perform a grid walkthrough and assemble the global matrix */
  // the coarse basis functions that belong to nodes on the boundary of the subgrid are included,
  // excluded are only basis functions that belong to nodes on the boundary of Omega
  void assemble ()  const
  {
    const DiscreteFunctionSpaceType &space = discreteFunctionSpace();

    // reserve memory for matrix
    linearOperator_.reserve();                                   /*@LST0E@*/

    // create timer (also stops time)
    Timer timer;

    // clear matrix                                             /*@LST0S@*/
    linearOperator_.clear();

    typedef typename DiscreteFunctionSpaceType :: IteratorType IteratorType;
    typedef typename CoarseDiscreteFunctionSpaceType :: IteratorType CoarseIteratorType;

    typedef typename IteratorType :: Entity EntityType;
    typedef typename CoarseIteratorType :: Entity CoarseEntityType;

    typedef typename EntityType :: Geometry GeometryType;
    typedef typename CoarseEntityType :: Geometry CoarseGeometryType;

    // coefficients in the matrix that describes the weighted Clement interpolation, i.e. coff[c] = (\int_{\Omega} \Phi_j)^{-1}
    std::vector<double> coff( coarse_space_.size(), 0.0 );

    CoarseIteratorType coarse_end = coarse_space_.end();
    for(CoarseIteratorType it = coarse_space_.begin(); it != coarse_end; ++it)
    {
      CoarseEntityType& entity = *it;

      std::vector< std::size_t > indices;
      coarse_space_.mapper().map(entity, indices);

      // cache geometry of entity
      const CoarseGeometryType coarse_geometry = entity.geometry();

      assert(entity.partitionType() == InteriorEntity);

      std::vector< RangeType > phi( coarse_space_.mapper().maxNumDofs() );

      // get base function set
      const CoarseBasisFunctionSetType &coarse_baseSet = coarse_space_.basisFunctionSet( entity );
      const auto numBaseFunctions = coarse_baseSet.size();

      // create quadrature of appropriate order
      CoarseQuadratureType quadrature( entity, 2 * polynomialOrder + 2 );

      // loop over all quadrature points
      const size_t numQuadraturePoints = quadrature.nop();
      for (size_t quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint)
      {
        const typename CoarseQuadratureType::CoordinateType& local_point = quadrature.point(quadraturePoint);

        const double weight = quadrature.weight(quadraturePoint) * coarse_geometry.integrationElement(local_point);

        coarse_baseSet.evaluateAll( quadrature[quadraturePoint], phi );
        for (unsigned int i = 0; i < numBaseFunctions; ++i)
         {
          coff[ indices[i] ] += weight * phi[i];
         }
      }
    }

    for ( size_t c = 0; c < coff.size(); ++c )
    {
      if ( coff[c] != 0.0 )
        coff[c] = 1.0 / coff[c];
    }

    for(const auto& sp_it : sparsity_pattern_.support())
    {
      const auto& entity = *sp_it.first;

       for(const auto& coarse_entity_ptr : sp_it.second)
       {

         const auto& coarse_entity = *coarse_entity_ptr;
          DSFe::LocalMatrixProxy<LinearOperatorType> localMatrix(linearOperator_, entity, coarse_entity, 1e-12);

          const CoarseGeometryType coarse_geometry = coarse_entity.geometry();

          // get base function set
          const CoarseBasisFunctionSetType &coarse_baseSet = coarse_space_.basisFunctionSet( coarse_entity );
          const auto coarse_numBaseFunctions = coarse_baseSet.size();

          const auto& coarse_lagrangepoint_set = specifier_.coarseSpace().lagrangePointSet(coarse_entity);

          // only implemented for 3 Lagrange Points, i.e. piecewise linear functions
         //! @todo Attention: 2D simplex only
          assert( coarse_numBaseFunctions == 3 );
          std::vector< RangeType > coarse_phi_corner_0( coarse_numBaseFunctions );
          std::vector< RangeType > coarse_phi_corner_1( coarse_numBaseFunctions );
          std::vector< RangeType > coarse_phi_corner_2( coarse_numBaseFunctions );

          std::vector< DomainType > coarse_corners( coarse_numBaseFunctions );
          std::vector< std::size_t > coarse_global_dof_number;
          coarse_space_.mapper().map(coarse_entity, coarse_global_dof_number);

          // coarse_corner_phi_j[i] = coarse basis function i evaluated in corner j
          coarse_baseSet.evaluateAll( coarse_lagrangepoint_set.point( 0 ) , coarse_phi_corner_0 );
          coarse_baseSet.evaluateAll( coarse_lagrangepoint_set.point( 1 ) , coarse_phi_corner_1 );
          coarse_baseSet.evaluateAll( coarse_lagrangepoint_set.point( 2 ) , coarse_phi_corner_2 );

          for(size_t loc_point = 0; loc_point < coarse_numBaseFunctions ; ++loc_point ) {
             coarse_corners[ loc_point ] = coarse_geometry.global(coarse_lagrangepoint_set.point( loc_point ) );
          }

          // LinearLagrangeInterpolation2D should be eventually replaced by LinearLagrangeFunction2D
          LinearLagrangeInterpolation2D< DiscreteFunctionSpaceType > coarse_basis_interpolation_0
             ( coarse_corners[0], coarse_phi_corner_0[0], coarse_corners[1], coarse_phi_corner_1[0], coarse_corners[2], coarse_phi_corner_2[0] );

          LinearLagrangeInterpolation2D< DiscreteFunctionSpaceType > coarse_basis_interpolation_1
             ( coarse_corners[0], coarse_phi_corner_0[1], coarse_corners[1], coarse_phi_corner_1[1], coarse_corners[2], coarse_phi_corner_2[1] );

          LinearLagrangeInterpolation2D< DiscreteFunctionSpaceType > coarse_basis_interpolation_2
             ( coarse_corners[0], coarse_phi_corner_0[2], coarse_corners[1], coarse_phi_corner_1[2], coarse_corners[2], coarse_phi_corner_2[2] );

          // cache geometry of entity
          const GeometryType geometry = entity.geometry();

          std::vector< RangeType > fine_phi( space.mapper().maxNumDofs() );

          // get base function set
          const BasisFunctionSetType &baseSet = space.basisFunctionSet( entity );
          const auto numBaseFunctions = baseSet.size();

          // create quadrature of appropriate order
          const QuadratureType quadrature( entity, 2 * polynomialOrder + 2 );

          // loop over all quadrature points
          const size_t numQuadraturePoints = quadrature.nop();
          for (size_t quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint)
          {
             const auto& local_point = quadrature.point(quadraturePoint);
             DomainType global_point = geometry.global( quadrature.point(quadraturePoint) );

             const double weight = quadrature.weight(quadraturePoint) * geometry.integrationElement(local_point);

             baseSet.evaluateAll( quadrature[quadraturePoint], fine_phi );

             for (unsigned int j = 0; j < numBaseFunctions; ++j)
             {
               for (unsigned int i = 0; i < coarse_numBaseFunctions; ++i)
               {
                  if ( specifier_.is_coarse_boundary_node( coarse_global_dof_number[i] ) )
                    { continue; }

                  RangeType coarse_phi_i( 0.0 );

                  if (i == 0)
                    coarse_basis_interpolation_0.evaluate(global_point, coarse_phi_i);
                  if (i == 1)
                    coarse_basis_interpolation_1.evaluate(global_point, coarse_phi_i);
                  if (i == 2)
                    coarse_basis_interpolation_2.evaluate(global_point, coarse_phi_i);

                  localMatrix.add( i, j, weight * coff[coarse_global_dof_number[i]] * coarse_phi_i * fine_phi[j] );
               }
             }
          }
       }
    }

    // get elapsed time
    const double assemblyTime = timer.elapsed();
    // in verbose mode print times
    if ( Fem::Parameter :: verbose () )
      std :: cout << "Time to assemble weighted clement operator: " << assemblyTime << "s" << std :: endl;

    // get grid sequence number from space (for adaptive runs)    /*@LST0S@*/
    sequence_ = dofManager_.sequence();

  }

  /** \attention This imp. differs from the one before the refactor in that it iterates over both domain and range entities
   *   in order to determine which local matrices to touch. I have no idea why this compiled before, since
   *   type(domain_entity) != type(range_entity) and therefore calling "localMatrix(entity, entity)"
   *   should never have been possible.
   **/
  void boundaryTreatment () const
  {
    for(const auto& entity : discreteFunctionSpace_)
    {
      for(const auto& coarse_entity : coarse_space_)
      {
        if (!Stuff::Grid::entities_identical(entity, coarse_entity))
          continue;

        // if entity has boundary intersections
        if( entity.hasBoundaryIntersections() )
        {
          // get local matrix from matrix object
          LocalMatrixType localMatrix = linearOperator_.localMatrix(entity, coarse_entity);

          const auto& lagrangePointSet = discreteFunctionSpace_.lagrangePointSet(entity);

          const auto endiit = discreteFunctionSpace_.gridPart().iend( entity );
          for( auto iit = discreteFunctionSpace_.gridPart().ibegin( entity );
               iit != endiit ; ++iit )
          {
            if ( iit->neighbor() ) // if there is a neighbor entity
              continue;

            const int face = (*iit).indexInInside();
            const auto fdend = lagrangePointSet.endSubEntity< 1 >(face);
            for (auto fdit = lagrangePointSet.beginSubEntity< 1 >(face); fdit != fdend; ++fdit)
              localMatrix.unitRow(*fdit);
          }
        }
      }
    }
  }

protected:
  const DiscreteFunctionSpaceType &discreteFunctionSpace_;
  const CoarseDiscreteFunctionSpaceType& coarse_space_;
  const DofManagerType &dofManager_;

  const MacroMicroGridSpecifier& specifier_;
  ClemementPattern<DiscreteFunctionSpaceType,CoarseDiscreteFunctionSpaceType> sparsity_pattern_;
  mutable LinearOperatorType linearOperator_;

  const CoarseNodeVectorType& coarse_nodes_;
  const CoarseBasisFunctionList& coarse_basis_;
  const std::map<int,int>& global_id_to_internal_id_;

  //! flag indicating whether the system matrix has been assembled
  mutable int sequence_;

  mutable Fem :: DynamicArray< JacobianRangeType > gradCache_;
  mutable Fem :: DynamicArray< RangeType > values_;
  mutable RangeFieldType weight_;
};

} //namespace MsFEM {
} //namespace Multiscale {
} //namespace Dune {

#endif
