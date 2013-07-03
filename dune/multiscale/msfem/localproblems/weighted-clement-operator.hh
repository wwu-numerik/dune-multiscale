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

#include <dune/common/fmatrix.hh>
#include <dune/fem/quadrature/quadrature.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/operator/matrix/spmatrix.hh>
#include <dune/fem/operator/2order/lagrangematrixsetup.hh>

#include <dune/multiscale/common/traits.hh>
#include <dune/multiscale/msfem/msfem_traits.hh>
#include <dune/multiscale/msfem/localproblems/subgrid-list.hh>
#include <dune/multiscale/msfem/localproblems/clement_pattern.hh>


namespace Dune {
namespace Multiscale {
namespace MsFEM {


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
                              const MacroMicroGridSpecifier& specifier );

private:
  // prohibit copying
  WeightedClementOperator ( const ThisType & ) = delete;

public:
  //! apply the operator
  virtual void operator() ( const DiscreteFunctionType &u,
                            CoarseDiscreteFunctionType &w ) const;

  //! return reference to preconditioning matrix, used by OEM-Solver
  const PreconditionMatrixType &preconditionMatrix () const;

  virtual void applyTransposed ( const CoarseDiscreteFunctionType &u,
                                 DiscreteFunctionType &w) const;

  //! return true if preconditioning is enabled
  bool hasPreconditionMatrix () const;

  //! print the system matrix into a stream
  void print ( std :: ostream & out = std :: cout ) const;

  //! return reference to discreteFunctionSpace
  const DiscreteFunctionSpaceType &discreteFunctionSpace () const;

  /*! \brief obtain a reference to the system matrix
   *
   *  The assembled matrix is returned. If the system matrix has not been
   *  assembled, yet, the assembly is performed.
   *
   *  \returns a reference to the system matrix
   */
  const LinearOperatorType &systemMatrix () const;

  /** \brief perform a grid walkthrough and assemble the global matrix */
  // the coarse basis functions that belong to nodes on the boundary of the subgrid are included,
  // excluded are only basis functions that belong to nodes on the boundary of Omega
  void assemble ()  const;

  /** \attention This imp. differs from the one before the refactor in that it iterates over both domain and range entities
   *   in order to determine which local matrices to touch. I have no idea why this compiled before, since
   *   type(domain_entity) != type(range_entity) and therefore calling "localMatrix(entity, entity)"
   *   should never have been possible.
   **/
  void boundaryTreatment () const;

protected:
  const DiscreteFunctionSpaceType& discreteFunctionSpace_;
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
