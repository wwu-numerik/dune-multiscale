#ifndef DUNE_MULTISCALE_LA_BACKEND_HH
#define DUNE_MULTISCALE_LA_BACKEND_HH

#include <dune/stuff/la/container.hh>
#include <dune/stuff/la/solver.hh>
#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/spaces/continuouslagrange.hh>

namespace Dune {
namespace Multiscale {

template <class DiscreteFunctionSpaceType>
struct BackendChooser {

  //! \todo gdt base
//  static_assert(std::is_base_of<GDT::Spaces::ContinuousLagrangeBase<
//                , GridPartImp::dimension, RangeFieldImp, 1, 1 >
//                                DiscreteFunctionSpaceType>::value, "");

  typedef DS::LocalizableFunctionInterface< typename DiscreteFunctionSpaceType::EntityType,
                                                  typename DiscreteFunctionSpaceType::DomainFieldType, DiscreteFunctionSpaceType::dimDomain,
                                                  typename DiscreteFunctionSpaceType::RangeFieldType, DiscreteFunctionSpaceType::dimRange,
     DiscreteFunctionSpaceType::dimRangeCols > DiscreteFunctionBaseType;

  typedef typename Stuff::LA::Container< typename DiscreteFunctionSpaceType::RangeFieldType,
    Stuff::LA::ChooseBackend::istl_sparse >::VectorType DiscreteFunctionDataType;
  typedef DiscreteFunctionDataType GdtVectorType;
  typedef typename Stuff::LA::Container< typename DiscreteFunctionSpaceType::RangeFieldType,
    Stuff::LA::ChooseBackend::istl_sparse >::MatrixType LinearOperatorType;
  typedef GDT::StoredDiscreteFunction< DiscreteFunctionSpaceType, DiscreteFunctionDataType >      DiscreteFunctionType;
  typedef GDT::ConstDiscreteFunction< DiscreteFunctionSpaceType, DiscreteFunctionDataType > ConstDiscreteFunctionType;
  typedef Stuff::LA::Solver< LinearOperatorType > InverseOperatorType;
};

} // namespace Multiscale
} // namespace Dune

#endif // DUNE_MULTISCALE_LA_BACKEND_HH
