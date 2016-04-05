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

template <class ct, int dim, template <int> class Refinement, class Comm>
class MyGridFactory<Dune::SPGrid<ct, dim, Refinement, Comm>>
    : public LocalGridFactory<Dune::SPGrid<ct, dim, Refinement, Comm>> {
  typedef Dune::SPGrid<ct, dim, Refinement, Comm> GridType;
  typedef typename GridType::ctype ctype;

public:
  static std::shared_ptr<GridType> createCubeGrid(const Dune::FieldVector<ctype, CommonTraits::world_dim>& lowerLeft,
                                                  const Dune::FieldVector<ctype, CommonTraits::world_dim>& upperRight,
                                                  const Dune::array<unsigned int, dim>& elements,
                                                  const Dune::array<unsigned int, dim>& overlap,
                                                  Dune::MPIHelper::MPICommunicator communicator) {
    return Dune::StructuredGridFactory<GridType>::createCubeGrid(lowerLeft, upperRight, elements, overlap,
                                                                 communicator);
  }
};

template <int dim, class Coord>
class MyGridFactory<Dune::YaspGrid<dim, Coord>> : public LocalGridFactory<Dune::YaspGrid<dim, Coord>> {
  typedef Dune::YaspGrid<dim, Coord> GridType;
  typedef typename GridType::ctype ctype;

public:
  static std::shared_ptr<GridType> createCubeGrid(const Dune::FieldVector<ctype, CommonTraits::world_dim>& lowerLeft,
                                                  const Dune::FieldVector<ctype, CommonTraits::world_dim>& upperRight,
                                                  const Dune::array<unsigned int, dim>& elements_in,
                                                  const Dune::array<unsigned int, dim>& overlap,
                                                  Dune::MPIHelper::MPICommunicator communicator) {
    const auto no_periodic_direction = std::bitset<dim>(false);
    Dune::array<int, dim> elements;
    std::copy(elements_in.begin(), elements_in.end(), elements.begin());
    auto overlap_check = overlap;
    overlap_check.fill(overlap[0]);
    for (auto i : Dune::XT::Common::value_range(1, dim))
      if (overlap[i] != overlap[0])
        DUNE_THROW(Dune::InvalidStateException, "YaspGrid only supports uniform overlap");
    return std::make_shared<GridType>(lowerLeft, upperRight, elements, no_periodic_direction, overlap[0], communicator);
  }
};

} // namespace Multiscale {
} // namespace Dune {

#endif // DUNE_MULTISCALE_COMMON_MYGRIDFACTORY_HH
