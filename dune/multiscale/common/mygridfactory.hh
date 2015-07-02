#ifndef DUNE_MULTISCALE_COMMON_MYGRIDFACTORY_HH
#define DUNE_MULTISCALE_COMMON_MYGRIDFACTORY_HH

#include <dune/multiscale/common/traits.hh>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/parallel/mpicollectivecommunication.hh>
#include <dune/grid/common/gridenums.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/spgrid.hh>
#include <dune/stuff/grid/structuredgridfactory.hh>

namespace Dune {
namespace Multiscale {

//! Helper struct to make overlap for SPGrid possible
template <class GridType>
class MyGridFactory {
  static_assert(Dune::AlwaysFalse<GridType>::value, "Only Yasp and SP grid are supported");
};

template <class GridType>
struct LocalGridFactory {
  typedef typename GridType::ctype ctype;

  static std::shared_ptr<GridType> createLocalGrid(const Dune::FieldVector<ctype, CommonTraits::world_dim>& lowerLeft,
                                                  const Dune::FieldVector<ctype, CommonTraits::world_dim>& upperRight,
                                                  const Dune::array<unsigned int, CommonTraits::world_dim>& elements) {

    Dune::array<unsigned int, CommonTraits::world_dim> overlap;
    overlap.fill(0u);
    return MyGridFactory<GridType>::createCubeGrid(lowerLeft, upperRight, elements, overlap,
                                                                 Dune::MPIHelper::getLocalCommunicator());
  }
};

template <class ct, int dim, Dune::SPRefinementStrategy strategy, class Comm>
class MyGridFactory<Dune::SPGrid<ct, dim, strategy, Comm>> : public LocalGridFactory<Dune::SPGrid<ct, dim, strategy, Comm>>{
  typedef Dune::SPGrid<ct, dim, strategy, Comm> GridType;
  typedef typename GridType::ctype ctype;

public:
  static std::shared_ptr<GridType> createCubeGrid(const Dune::FieldVector<ctype, CommonTraits::world_dim>& lowerLeft,
                                                  const Dune::FieldVector<ctype, CommonTraits::world_dim>& upperRight,
                                                  const Dune::array<unsigned int, dim>& elements,
                                                  const Dune::array<unsigned int, dim>& overlap,
                                                  Dune::MPIHelper::MPICommunicator communicator) {
    DUNE_THROW(InvalidStateException, "");
    return Dune::StructuredGridFactory<GridType>::createCubeGrid(lowerLeft, upperRight, elements, overlap,
                                                                 communicator);
  }

};

template <int dim>
class MyGridFactory<Dune::YaspGrid<dim>> : public LocalGridFactory<Dune::YaspGrid<dim>> {
  typedef Dune::YaspGrid<dim> GridType;
  typedef typename GridType::ctype ctype;

public:
  static std::shared_ptr<GridType> createCubeGrid(const Dune::FieldVector<ctype, CommonTraits::world_dim>& lowerLeft,
                                                  const Dune::FieldVector<ctype, CommonTraits::world_dim>& upperRight,
                                                  const Dune::array<unsigned int, dim>& elements_in,
                                                  const Dune::array<unsigned int, dim>& overlap,
                                                  Dune::MPIHelper::MPICommunicator communicator) {
    const auto no_periodic_direction = std::bitset<dim>();
    if(DSC::FloatCmp::ne(lowerLeft, Dune::FieldVector<ctype, CommonTraits::world_dim>(0.0)))
      DUNE_THROW(Dune::InvalidStateException, "YaspGrid + Origin != 0.0 is still a no-go");
    Dune::array<int, dim> elements;
    std::copy(elements_in.begin(), elements_in.end(), elements.begin());
    auto overlap_check = overlap;
    overlap_check.fill(overlap[0]);
    for(auto i : DSC::valueRange(1,dim))
      if(overlap[i] != overlap[0])
        DUNE_THROW(Dune::InvalidStateException, "YaspGrid only supports uniform overlap");
    return std::make_shared<GridType>(communicator, upperRight, elements,
                      no_periodic_direction, overlap[0]);
  }
};

} // namespace Multiscale {
} // namespace Dune {

#endif // DUNE_MULTISCALE_COMMON_MYGRIDFACTORY_HH
