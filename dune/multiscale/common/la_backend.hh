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

  typedef typename Stuff::LA::Container< typename DiscreteFunctionSpaceType::RangeFieldType,
    Stuff::LA::ChooseBackend::istl_sparse >::VectorType DiscreteFunctionDataType;
  typedef DiscreteFunctionDataType GdtVectorType;
  typedef typename Stuff::LA::Container< typename DiscreteFunctionSpaceType::RangeFieldType,
    Stuff::LA::ChooseBackend::istl_sparse >::MatrixType LinearOperatorType;
  typedef GDT::DiscreteFunction< DiscreteFunctionSpaceType, DiscreteFunctionDataType >      DiscreteFunctionType;
  typedef GDT::ConstDiscreteFunction< DiscreteFunctionSpaceType, DiscreteFunctionDataType > ConstDiscreteFunctionType;
  typedef Stuff::LA::Solver< LinearOperatorType, typename DiscreteFunctionSpaceType::CommunicatorType > InverseOperatorType;
};

} // namespace Multiscale
} // namespace Dune

#endif // DUNE_MULTISCALE_LA_BACKEND_HH
