#ifndef RB_DUNE_FEM_REDUCEDBASISSPACE_BASEFUNCTIONSET_HH
#define RB_DUNE_FEM_REDUCEDBASISSPACE_BASEFUNCTIONSET_HH

#include <iostream>

#include <dune/fem/storage/array.hh>
#include <dune/fem/space/basefunctions/basefunctioninterface.hh>
#include <dune/fem/space/basefunctions/basefunctionsetinterface.hh>
#include "../discfunclist/discfunclist_xdr.hh"
#include "../discfunclist/discfunclist_mem.hh"

using namespace Dune;

namespace Dune {
namespace Multiscale {
template< class BaseFunctionImp >
class ReducedBasisBaseFunctionSet;

template< class BaseFunctionImp >
class ReducedBasisBaseFunctionSetTraits
{
public:
  typedef BaseFunctionImp BaseFunctionType;

private:
  typedef ReducedBasisBaseFunctionSetTraits< BaseFunctionType > ThisType;

public:
  typedef ReducedBasisBaseFunctionSet< BaseFunctionImp > BaseFunctionSetType;

  typedef typename BaseFunctionType::DiscreteFunctionSpaceType BaseFunctionSpaceType;

/*    typedef BaseFunctionSpaceType                                    DiscreteFunctionSpaceType;*/

  typedef BaseFunctionType DiscreteFunctionType;

  typedef typename BaseFunctionSpaceType::BaseFunctionSetType
    ::FunctionSpaceType FunctionSpaceType;
  typedef double AttributeType;
};

/** \class ReducedBasisBaseFunctionSet
   *  \brief The ReducedBasisBaseFunctionSet class provides
   *
   *  This class is needed to build the space and provides the functionality of
   *  the space for example the jacobian method is implemented here
   */
template< class BaseFunctionImp >
class ReducedBasisBaseFunctionSet
  : public BaseFunctionSetDefault< ReducedBasisBaseFunctionSetTraits< BaseFunctionImp > >
{
private:
  typedef ReducedBasisBaseFunctionSet< BaseFunctionImp > ThisType;

public:
  typedef ThisType BaseFunctionSetType;

  typedef ReducedBasisBaseFunctionSetTraits< BaseFunctionImp > Traits;

  typedef BaseFunctionSetDefault< Traits > BaseType;

  typedef typename Traits::BaseFunctionType            BaseFunctionType;
  typedef typename Traits::BaseFunctionSpaceType       BaseFunctionSpaceType;
  typedef typename BaseFunctionType::LocalFunctionType LocalBaseFunctionType;

  typedef typename BaseFunctionSpaceType::GridPartType GridPartType;

  typedef typename Traits::FunctionSpaceType FunctionSpaceType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType  RangeFieldType;

  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType  RangeType;

  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
  typedef typename FunctionSpaceType::HessianRangeType  HessianRangeType;

  static const int dimDomain = FunctionSpaceType::dimDomain;
  static const int dimRange = FunctionSpaceType::dimRange;

  typedef DiscreteFunctionList_xdr< Traits > BaseFunctionListType;
  // typedef DynamicArray< BaseFunctionType* >                        BaseFunctionListType;

private:
  typedef typename GridPartType::GridType
    ::template Codim< 0 >::Entity EntityCodim0Type;

public:
  /** \brief default constructor */
  inline ReducedBasisBaseFunctionSet()
    : baseFunctionList_(NULL)
      , entity_(NULL)
  {}

  /** constructor
     *
     *  This constructor initializes the base function set, but does not bind
     *  it to an entity.
     *
     *  \param[in]  baseFunctionList  array containing the discrete functions
     *                                to be used as base functions
     */
  inline explicit ReducedBasisBaseFunctionSet(const BaseFunctionListType& baseFunctionList)
    : baseFunctionList_(&baseFunctionList)
      , entity_(NULL)
  {}

  /** constructor
     *
     *  This constructor initializes the base function set and binds it to an
     *  entity.
     *
     *  \param[in]  baseFunctionList  array containing the discrete functions
     *                                to be used as base functions
     *  \param[in]  entity            entity (of codim 0) to bind the base
     *                                function set to
     */
  inline ReducedBasisBaseFunctionSet(const BaseFunctionListType& baseFunctionList,
                                     const EntityCodim0Type& entity)
    : baseFunctionList_(&baseFunctionList)
      , entity_(&entity)
  {}

  /** \brief copy constructor
     *
     *  \param[in]  other  base function set to copy
     */
  inline ReducedBasisBaseFunctionSet(const ThisType& other)
    : baseFunctionList_(other.baseFunctionList_)
      , entity_(other.entity_)
  {}

  /** \brief copy another ReducedBasisBaseFunctionSet
     *
     *  \param[in]  other  base function set to copy
     */
  inline ThisType& operator=(const ThisType& other) {
    baseFunctionList_ = other.baseFunctionList_;
    entity_ = other.entity_;
    return *this;
  }

  /* \copydoc Dune::RB::BaseFunctionSetInterface::evaluate(const int,const FieldVector<deriType,diffOrd> &,const
     *PointType &,RangeType &)const */
  template< int diffOrd, class deriType, class PointType >
  inline void evaluate(const int baseFunction,
                       const FieldVector< deriType, diffOrd >& diffVariable,
                       const PointType& x,
                       RangeType& phi) const;

// !NEW

  /* \copydoc Dune::RB::BaseFunctionSetInterface::evaluate(const int,const FieldVector<deriType,diffOrd> &,const
     *PointType &,RangeType &)const */
  template< class PointType >
  inline void evaluate(const int baseFunction,
                       const PointType& x,
                       RangeType& phi) const;

  /* \copydoc Dune::RB::BaseFunctionSetInterface::evaluateSingle(const int baseFunction,const PointType &x,const
     *RangeType &psi) const */
  template< class PointType >
  inline RangeFieldType evaluateSingle(const int baseFunction,
                                       const PointType& x,
                                       RangeType& psi) const;

  /** \brief obtain the entity, this base function set belongs to */
  inline const EntityCodim0Type& entity() const {
    assert(entity_ != NULL);
    return *entity_;
  }

  inline GeometryType geometryType() const {
    return entity().geometry().type();
  }

  /* \copydoc Dune::RB::BaseFunctionSetInterface::numBaseFunctions() */
  inline int numBaseFunctions() const {
    assert(baseFunctionList_ != NULL);
    return baseFunctionList_->size();
  }

protected:
  const BaseFunctionListType* baseFunctionList_;
  const EntityCodim0Type* entity_;

  int numBaseFunctions_;
};
}
}

#include "basefunctionset_inline.hh"

#endif // ifndef RB_DUNE_FEM_REDUCEDBASISSPACE_BASEFUNCTIONSET_HH
