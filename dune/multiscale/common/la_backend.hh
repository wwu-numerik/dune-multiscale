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

#if DUNE_MULTISCALE_USE_ISTL
  static constexpr auto backend_type = Stuff::LA::ChooseBackend::istl_sparse;
#else
  static constexpr auto backend_type = Stuff::LA::ChooseBackend::eigen_sparse;
#endif

  typedef typename Stuff::LA::Container<typename DiscreteFunctionSpaceType::RangeFieldType, backend_type>::VectorType
  DiscreteFunctionDataType;
  typedef DiscreteFunctionDataType GdtVectorType;
  typedef typename Stuff::LA::Container<typename DiscreteFunctionSpaceType::RangeFieldType, backend_type>::MatrixType
  LinearOperatorType;
  typedef GDT::DiscreteFunction<DiscreteFunctionSpaceType, DiscreteFunctionDataType> DiscreteFunctionType;
  typedef GDT::ConstDiscreteFunction<DiscreteFunctionSpaceType, DiscreteFunctionDataType> ConstDiscreteFunctionType;
  typedef Stuff::LA::Solver<LinearOperatorType, typename DiscreteFunctionSpaceType::CommunicatorType>
  InverseOperatorType;
};

} // namespace Multiscale
} // namespace Dune

#endif // DUNE_MULTISCALE_LA_BACKEND_HH
