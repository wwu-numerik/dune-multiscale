#ifndef DUNE_MULTISCALE_LA_BACKEND_HH
#define DUNE_MULTISCALE_LA_BACKEND_HH

#include <dune/xt/la/container.hh>
#include <dune/xt/la/solver.hh>
#include <dune/xt/la/solver/istl.hh>
#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/spaces/cg.hh>

namespace Dune {
namespace Multiscale {

template <class DiscreteFunctionSpaceType>
struct BackendChooser
{

#if DUNE_MULTISCALE_USE_ISTL
  static constexpr auto backend_type = XT::LA::Backends::istl_sparse;
#else
  static constexpr auto backend_type = XT::LA::Backends::eigen_sparse;
#endif

  typedef typename XT::LA::Container<typename DiscreteFunctionSpaceType::RangeFieldType, backend_type>::VectorType
      DiscreteFunctionDataType;
  typedef DiscreteFunctionDataType GdtVectorType;
  typedef typename XT::LA::Container<typename DiscreteFunctionSpaceType::RangeFieldType, backend_type>::MatrixType
      LinearOperatorType;
  typedef GDT::DiscreteFunction<DiscreteFunctionSpaceType, DiscreteFunctionDataType> DiscreteFunctionType;
  typedef GDT::ConstDiscreteFunction<DiscreteFunctionSpaceType, DiscreteFunctionDataType> ConstDiscreteFunctionType;
  typedef XT::LA::Solver<LinearOperatorType, typename DiscreteFunctionSpaceType::DofCommunicatorType>
      InverseOperatorType;
};

} // namespace Multiscale
} // namespace Dune

#endif // DUNE_MULTISCALE_LA_BACKEND_HH
