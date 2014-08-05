#ifndef DUNE_MULTISCALE_MSFEM_DIFFUSION_EVALUATION_HH
#define DUNE_MULTISCALE_MSFEM_DIFFUSION_EVALUATION_HH

#include <dune/multiscale/msfem/msfem_traits.hh>
#include <dune/multiscale/problems/base.hh>
#include <dune/gdt/functionals/l2.hh>

namespace Dune {
namespace Multiscale {
namespace MsFEM {

// forward, to be used in the traits
template< class LocalizableFunctionImp >
class CoarseBasisProduct;


/**
 *  \brief Traits for the Product evaluation.
 */
template< class LocalizableFunctionImp >
class CoarseBasisProductTraits
{
public:
  typedef CoarseBasisProduct< LocalizableFunctionImp > derived_type;
  typedef LocalizableFunctionImp            LocalizableFunctionType;
  static_assert(std::is_base_of< Dune::Stuff::IsLocalizableFunction, LocalizableFunctionImp >::value,
                "LocalizableFunctionImp has to be derived from Stuff::IsLocalizableFunction.");
};


template< class LocalizableFunctionImp >
class CoarseBasisProduct
  : public GDT::LocalEvaluation::Codim0Interface< CoarseBasisProductTraits< LocalizableFunctionImp >, 1 >
{
public:
  typedef CoarseBasisProductTraits< LocalizableFunctionImp >   Traits;
  typedef typename Traits::LocalizableFunctionType  LocalizableFunctionType;

  CoarseBasisProduct(const Multiscale::CommonTraits::BaseFunctionSetType& coarse_base,
                     const LocalizableFunctionType& inducingFunction,
                     const std::size_t coarseBaseFunc)
    : inducingFunction_(inducingFunction)
    , coarse_base_set_(coarse_base)
    , coarseBaseFunc_(coarseBaseFunc)
  {}

  template< class EntityType >
  class LocalfunctionTuple
  {
    typedef typename LocalizableFunctionType::LocalfunctionType LocalfunctionType;
  public:
    typedef std::tuple< std::shared_ptr< LocalfunctionType > > Type;
  };

  template< class EntityType >
  typename LocalfunctionTuple< EntityType >::Type localFunctions(const EntityType& entity) const
  {
    return std::make_tuple(inducingFunction_.local_function(entity));
  }

  /**
   * \brief extracts the local functions and calls the correct order() method
   */
  template< class E, class D, int d, class R, int rT, int rCT >
  size_t order(const typename LocalfunctionTuple< E >::Type& localFuncs,
               const Stuff::LocalfunctionSetInterface< E, D, d, R, rT, rCT >& testBase) const
  {
    const auto localFunction = std::get< 0 >(localFuncs);
    return order(*localFunction, testBase);
  }

  /**
   *  \todo add copydoc
   *  \return localFunction.order() + testBase.order()
   */
  template< class E, class D, int d, class R, int rL, int rCL, int rT, int rCT >
  size_t order(const Stuff::LocalfunctionInterface< E, D, d, R, rL, rCL >& localFunction,
               const Stuff::LocalfunctionSetInterface< E, D, d, R, rT, rCT >& testBase) const
  {
    return localFunction.order() + testBase.order();
  } // int order(...)

  /**
   * \brief extracts the local functions and calls the correct evaluate() method
   */
  template< class E, class D, int d, class R, int rT, int rCT >
   void evaluate(const typename LocalfunctionTuple< E >::Type& localFuncs,
                 const Stuff::LocalfunctionSetInterface< E, D, d, R, rT, rCT >& testBase,
                 const Dune::FieldVector< D, d >& localPoint,
                 Dune::DynamicVector< R >& ret) const
  {
    const auto& localFunction = std::get< 0 >(localFuncs);
    evaluate(*localFunction, testBase, localPoint, ret);
  }

  template< class E, class D, int d, class R, int rL, int rCL, int rT, int rCT >
  void evaluate(const Stuff::LocalfunctionInterface< E, D, d, R, rL, rCL >& localFunction,
                const Stuff::LocalfunctionSetInterface< E, D, d, R, rT, rCT >& testBase,
                const Dune::FieldVector< D, d >& localPoint,
                Dune::DynamicVector< R >& ret) const
  {
    typedef Dune::FieldVector< R, 1 > RangeType;
    // evaluate local function
    const auto& entity = testBase.entity();
    const auto global_point = entity.geometry().global(localPoint);
    const auto& coarse_entity = coarse_base_set_.entity();
    const auto quadInCoarseLocal = coarse_entity.geometry().local(global_point);
    const auto coarseBaseFuncJacs = coarse_base_set_.jacobian(quadInCoarseLocal);
    const auto direction = coarseBaseFuncJacs[coarseBaseFunc_];
    DMP::DiffusionBase::RangeType functionValue;
    localFunction.evaluate(localPoint, functionValue);
    functionValue[0][0] = functionValue[0][0] * direction[0][0];
    functionValue[0][1] = functionValue[1][1] * direction[0][1];
    // evaluate test base
    const std::size_t size = testBase.size();
    std::vector< RangeType > testValues = testBase.evaluate(localPoint);
    // compute product
    assert(ret.size() >= size);

    for (size_t ii = 0; ii < size; ++ii) {
      auto fo = functionValue[0];
      fo *=testValues[ii][0];
      ret[ii] =  fo[ii];
    }
  }


private:
  const LocalizableFunctionType& inducingFunction_;
  const Multiscale::CommonTraits::BaseFunctionSetType& coarse_base_set_;
  const std::size_t coarseBaseFunc_;
}; // class Product

} // namespace MsFEM {
} // namespace Multiscale {
} // namespace Dune {

#endif // DUNE_MULTISCALE_MSFEM_DIFFUSION_EVALUATION_HH
