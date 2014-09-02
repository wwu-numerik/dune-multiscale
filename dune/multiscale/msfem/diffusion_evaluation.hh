#ifndef DUNE_MULTISCALE_MSFEM_DIFFUSION_EVALUATION_HH
#define DUNE_MULTISCALE_MSFEM_DIFFUSION_EVALUATION_HH

#include <dune/multiscale/msfem/msfem_traits.hh>
#include <dune/multiscale/problems/base.hh>
#include <dune/multiscale/problems/selector.hh>
#include <dune/gdt/functionals/l2.hh>

namespace Dune {
namespace Multiscale {
namespace MsFEM {

// forward, to be used in the traits
class CoarseBasisProduct;

/**
 *  \brief Traits for the Product evaluation.
 */
class CoarseBasisProductTraits {
public:
  typedef CoarseBasisProduct derived_type;
  typedef Problem::LocalDiffusionType LocalizableFunctionType;
  static_assert(std::is_base_of<Dune::Stuff::IsLocalizableFunction, LocalizableFunctionType>::value,
                "LocalizableFunctionImp has to be derived from Stuff::IsLocalizableFunction.");
};

class CoarseBasisProduct : public GDT::LocalEvaluation::Codim0Interface<CoarseBasisProductTraits, 1> {
  typedef MsFEMTraits::LocalEntityType EntityType;

public:
  typedef CoarseBasisProductTraits Traits;
  typedef typename Traits::LocalizableFunctionType LocalizableFunctionType;

  CoarseBasisProduct(const Multiscale::CommonTraits::BaseFunctionSetType& coarse_base,
                     const LocalizableFunctionType& inducingFunction, const std::size_t coarseBaseFunc)
    : inducingFunction_(inducingFunction)
    , coarse_base_set_(coarse_base)
    , coarseBaseFunc_(coarseBaseFunc) {}

  class LocalfunctionTuple {
    typedef typename LocalizableFunctionType::LocalfunctionType LocalfunctionType;

  public:
    typedef std::tuple<std::shared_ptr<LocalfunctionType>> Type;
  };

  typename LocalfunctionTuple::Type localFunctions(const EntityType& entity) const {
    return std::make_tuple(inducingFunction_.local_function(entity));
  }

  /**
   * \brief extracts the local functions and calls the correct order() method
   */
  template <class E, class D, int d, class R, int rT, int rCT>
  size_t order(const typename LocalfunctionTuple::Type& localFuncs,
               const Stuff::LocalfunctionSetInterface<E, D, d, R, rT, rCT>& testBase) const {
    const auto localFunction = std::get<0>(localFuncs);
    return order(*localFunction, testBase);
  }

  /**
   *  \todo add copydoc
   *  \return localFunction.order() + testBase.order()
   */
  template <class E, class D, int d, class R, int rL, int rCL, int rT, int rCT>
  size_t order(const Stuff::LocalfunctionInterface<E, D, d, R, rL, rCL>& localFunction,
               const Stuff::LocalfunctionSetInterface<E, D, d, R, rT, rCT>& testBase) const {
    return localFunction.order() + testBase.order();
  } // int order(...)

  /**
   * \brief extracts the local functions and calls the correct evaluate() method
   */
  template <class D, int d, class R, int rT, int rCT>
  void evaluate(const typename LocalfunctionTuple::Type& localFuncs,
                const Stuff::LocalfunctionSetInterface<EntityType, D, d, R, rT, rCT>& testBase,
                const Dune::FieldVector<D, d>& localPoint, Dune::DynamicVector<R>& ret) const {
    const auto& localFunction = std::get<0>(localFuncs);
    evaluate(*localFunction, testBase, localPoint, ret);
  }

  template <class E, class D, int d, class R, int rL, int rCL, int rT, int rCT>
  void evaluate(const Stuff::LocalfunctionInterface<E, D, d, R, rL, rCL>& /*localFunction*/,
                const Stuff::LocalfunctionSetInterface<E, D, d, R, rT, rCT>& testBase,
                const Dune::FieldVector<D, d>& localPoint, Dune::DynamicVector<R>& ret) const {
    // evaluate local function
    const auto& entity = testBase.entity();
    const auto global_point = entity.geometry().global(localPoint);
    const auto& coarse_entity = coarse_base_set_.entity();
    const auto quadInCoarseLocal = coarse_entity.geometry().local(global_point);
    const auto coarseBaseFuncJacs = coarse_base_set_.jacobian(quadInCoarseLocal);
    const auto direction = coarseBaseFuncJacs[coarseBaseFunc_];

    DMP::JacobianRangeType flux;
    DMP::getDiffusion()->diffusiveFlux(global_point, direction, flux);
    // evaluate test base
    const std::size_t size = testBase.size();
    const auto transformed_gradients = testBase.jacobian(localPoint);
    // compute product
    assert(ret.size() >= size);
    assert(transformed_gradients.size() >= size);

    //! \TODO WTF muss hier eigentlich hin
    for (size_t ii = 0; ii < size; ++ii) {
      // transformed_gradients[ii] is FieldMatrix<double, 1, 2> --> grad_phi_s[ii][0] is FieldVector<double,2>
      ret[ii] = -1 * (flux[0] * transformed_gradients[ii][0]);
      //      DSC_LOG_DEBUG << "DIFF " << global_point << " | " << flux  << " | " << grad_phi_s[ii][0]<< std::endl;
    }
  }

private:
  const LocalizableFunctionType& inducingFunction_;
  const Multiscale::CommonTraits::BaseFunctionSetType& coarse_base_set_;
  const std::size_t coarseBaseFunc_;
}; // class CoarseBasisProduct

// forward, to be used in the traits
class DirichletProduct;

/**
 *  \brief Traits for the Product evaluation.
 */
class DirichletProductTraits {
public:
  typedef DirichletProduct derived_type;
  typedef Problem::LocalDiffusionType LocalizableFunctionType;
  static_assert(std::is_base_of<Dune::Stuff::IsLocalizableFunction, LocalizableFunctionType>::value,
                "LocalizableFunctionImp has to be derived from Stuff::IsLocalizableFunction.");
};

class DirichletProduct : public GDT::LocalEvaluation::Codim0Interface<DirichletProductTraits, 1> {
public:
  typedef DirichletProductTraits Traits;
  typedef typename Traits::LocalizableFunctionType LocalizableFunctionType;

  DirichletProduct(const MsFEMTraits::LocalGridDiscreteFunctionType& dirichlet_extension,
                   const LocalizableFunctionType& inducingFunction)
    : inducingFunction_(inducingFunction)
    , dirichlet_extension_(dirichlet_extension) {}

  template <class EntityType>
  class LocalfunctionTuple {
    typedef typename LocalizableFunctionType::LocalfunctionType LocalfunctionType;

  public:
    typedef std::tuple<std::shared_ptr<LocalfunctionType>> Type;
  };

  template <class EntityType>
  typename LocalfunctionTuple<EntityType>::Type localFunctions(const EntityType& entity) const {
    return std::make_tuple(inducingFunction_.local_function(entity));
  }

  /**
   * \brief extracts the local functions and calls the correct order() method
   */
  template <class E, class D, int d, class R, int rT, int rCT>
  size_t order(const typename LocalfunctionTuple<E>::Type& localFuncs,
               const Stuff::LocalfunctionSetInterface<E, D, d, R, rT, rCT>& testBase) const {
    const auto localFunction = std::get<0>(localFuncs);
    return order(*localFunction, testBase);
  }

  /**
   *  \todo add copydoc
   *  \return localFunction.order() + testBase.order()
   */
  template <class E, class D, int d, class R, int rL, int rCL, int rT, int rCT>
  size_t order(const Stuff::LocalfunctionInterface<E, D, d, R, rL, rCL>& localFunction,
               const Stuff::LocalfunctionSetInterface<E, D, d, R, rT, rCT>& testBase) const {
    return localFunction.order() + testBase.order();
  } // int order(...)

  /**
   * \brief extracts the local functions and calls the correct evaluate() method
   */
  template <class E, class D, int d, class R, int rT, int rCT>
  void evaluate(const typename LocalfunctionTuple<E>::Type& localFuncs,
                const Stuff::LocalfunctionSetInterface<E, D, d, R, rT, rCT>& testBase,
                const Dune::FieldVector<D, d>& localPoint, Dune::DynamicVector<R>& ret) const {
    const auto& localFunction = std::get<0>(localFuncs);
    evaluate(*localFunction, testBase, localPoint, ret);
  }

  template <class E, class D, int d, class R, int rL, int rCL, int rT, int rCT>
  void evaluate(const Stuff::LocalfunctionInterface<E, D, d, R, rL, rCL>& /*localFunction*/,
                const Stuff::LocalfunctionSetInterface<E, D, d, R, rT, rCT>& testBase,
                const Dune::FieldVector<D, d>& localPoint, Dune::DynamicVector<R>& ret) const {
    // evaluate local function
    const auto& entity = testBase.entity();
    const auto dirichlet_lf = dirichlet_extension_.local_function(entity);
    const auto direction = dirichlet_lf->jacobian(localPoint);

    DMP::JacobianRangeType flux;
    const auto global_point = entity.geometry().global(localPoint);
    DMP::getDiffusion()->diffusiveFlux(global_point, direction, flux);
    // evaluate test base
    const std::size_t size = testBase.size();
    typedef typename Stuff::LocalfunctionSetInterface<E, D, d, R, rT, rCT>::JacobianRangeType JR;
    const std::vector<JR> grad_phi_s = testBase.jacobian(localPoint);

    // compute product
    assert(ret.size() >= size);
    assert(grad_phi_s.size() >= size);

    //! \TODO WTF muss hier eigentlich hin
    for (size_t ii = 0; ii < size; ++ii) {
      // grad_phi_s[ii] is FieldMatrix<double, 1, 2> --> grad_phi_s[ii][0] is FieldVector<double,2>
      ret[ii] = -1 * (flux[0] * grad_phi_s[ii][0]);
    }
  }

private:
  const LocalizableFunctionType& inducingFunction_;
  const MsFEMTraits::LocalGridDiscreteFunctionType& dirichlet_extension_;
}; // class Product

} // namespace MsFEM {
} // namespace Multiscale {
} // namespace Dune {

#endif // DUNE_MULTISCALE_MSFEM_DIFFUSION_EVALUATION_HH
