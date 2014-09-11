#include <config.h>
#include "grid_creation.hh"

#include <dune/grid/common/gridenums.hh>

#include <dune/stuff/common/ranges.hh>
#include <dune/stuff/grid/structuredgridfactory.hh>
#include <dune/stuff/grid/information.hh>
#include <dune/multiscale/problems/selector.hh>
#include <dune/common/parallel/mpihelper.hh>

// Helper struct to make overlap for SPGrid possible
// Declared in unnamed namespace to avoid naming conflicts
namespace {
template <class GridType>
class MyGridFactory {
  typedef typename GridType::ctype ctype;
  static const int dimworld = GridType::dimensionworld;
  static const int dim = GridType::dimension;

public:
  static std::shared_ptr<GridType> createCubeGrid(const Dune::FieldVector<ctype, dimworld>& lowerLeft,
                                                  const Dune::FieldVector<ctype, dimworld>& upperRight,
                                                  const Dune::array<unsigned int, dim>& elements,
                                                  const Dune::array<unsigned int, dim>& /*overlap*/) {
    // structured grid factory allows overlap only for SPGrid at the moment, hence the following check
    BOOST_ASSERT_MSG((DSC_CONFIG_GET("msfem.oversampling_layers", 0) == 0),
                     "Oversampling may only be used in combination with SPGrid!");
    return Dune::StructuredGridFactory<GridType>::createCubeGrid(lowerLeft, upperRight, elements);
  }
};

template <class ct, int dim, Dune::SPRefinementStrategy strategy, class Comm>
class MyGridFactory<Dune::SPGrid<ct, dim, strategy, Comm>> {
  typedef Dune::SPGrid<ct, dim, strategy, Comm> GridType;
  typedef typename GridType::ctype ctype;
  static const int dimworld = GridType::dimensionworld;

public:
  static std::shared_ptr<GridType> createCubeGrid(const Dune::FieldVector<ctype, dimworld>& lowerLeft,
                                                  const Dune::FieldVector<ctype, dimworld>& upperRight,
                                                  const Dune::array<unsigned int, dim>& elements,
                                                  const Dune::array<unsigned int, dim>& overlap) {
    return Dune::StructuredGridFactory<GridType>::createCubeGrid(lowerLeft, upperRight, elements, overlap);
  }
};
}

using namespace Dune::Multiscale;
using namespace std;

typedef tuple<CommonTraits::DomainType, CommonTraits::DomainType, array<unsigned int, CommonTraits::world_dim>,
              array<unsigned int, CommonTraits::world_dim>,
              array<unsigned int, CommonTraits::world_dim>> SetupReturnType;

SetupReturnType setup() {
  BOOST_ASSERT_MSG(DSC_CONFIG.has_sub("grids"), "Parameter tree needs to have 'grids' subtree!");

  const auto world_dim = CommonTraits::world_dim;
  typedef CommonTraits::DomainType CoordType;
  const auto& gridCorners = Problem::getModelData()->gridCorners();
  CoordType lowerLeft = gridCorners.first;
  CoordType upperRight = gridCorners.second;

  const auto oversamplingLayers = DSC_CONFIG_GET("msfem.oversampling_layers", 0);
  const auto microPerMacro = DSC_CONFIG.get<CoordType>("grids.micro_cells_per_macrocell_dim", CoordType(8), world_dim);
  const auto coarse_cells = DSC_CONFIG.get<CoordType>("grids.macro_cells_per_dim", CoordType(8), world_dim);

  array<unsigned int, world_dim> elements, overCoarse, overFine;

  for (const auto i : DSC::valueRange(world_dim)) {
    elements[i] = coarse_cells[i];
    overCoarse[i] = std::ceil(double(oversamplingLayers) / double(microPerMacro[i]));
    overFine[i] = DSC_CONFIG_GET("grids.overlap", 1);
  }
  return std::make_tuple(lowerLeft, upperRight, elements, overCoarse, overFine);
}

std::shared_ptr<CommonTraits::GridType> Dune::Multiscale::make_coarse_grid() {
  CommonTraits::DomainType lowerLeft, upperRight;
  array<unsigned int, CommonTraits::world_dim> elements, overCoarse;
  std::tie(lowerLeft, upperRight, elements, overCoarse, std::ignore) = setup();
  auto coarse_gridptr =
      MyGridFactory<CommonTraits::GridType>::createCubeGrid(lowerLeft, upperRight, elements, overCoarse);
  const auto expected_elements = std::accumulate(elements.begin(), elements.end(), 1u, std::multiplies<unsigned int>());
  auto actual_elements = coarse_gridptr->size(0);
  if (int(expected_elements) != coarse_gridptr->comm().sum(actual_elements))
    DUNE_THROW(InvalidStateException, "Wonky grid distribution");
  return coarse_gridptr;
}

pair<shared_ptr<CommonTraits::GridType>, shared_ptr<CommonTraits::GridType>>
Dune::Multiscale::make_grids(const bool check_partitioning) {
  auto coarse_grid = make_coarse_grid();
  return {coarse_grid, make_fine_grid(coarse_grid, check_partitioning)};
}

std::shared_ptr<Dune::Multiscale::CommonTraits::GridType>
Dune::Multiscale::make_fine_grid(std::shared_ptr<Dune::Multiscale::CommonTraits::GridType> coarse_gridptr,
                                 const bool check_partitioning) {
  const auto world_dim = CommonTraits::world_dim;
  CommonTraits::DomainType lowerLeft, upperRight;
  array<unsigned int, world_dim> elements, overFine;
  std::tie(lowerLeft, upperRight, elements, std::ignore, overFine) = setup();
  const auto coarse_cells =
      DSC_CONFIG.get<CommonTraits::DomainType>("grids.macro_cells_per_dim", CommonTraits::DomainType(8), world_dim);
  const auto microPerMacro = DSC_CONFIG.get<CommonTraits::DomainType>("grids.micro_cells_per_macrocell_dim",
                                                                      CommonTraits::DomainType(8), world_dim);

  for (const auto i : DSC::valueRange(CommonTraits::world_dim)) {
    elements[i] = coarse_cells[i] * microPerMacro[i];
  }
  auto fine_gridptr =
      StructuredGridFactory<CommonTraits::GridType>::createCubeGrid(lowerLeft, upperRight, elements, overFine);

  // check whether grids match (may not match after load balancing if different refinements in different
  // spatial directions are used)
  DSC_LOG_DEBUG << boost::format("Rank %d has %d coarse codim-0 elements and %d fine ones\n") %
                       coarse_gridptr->comm().rank() % coarse_gridptr->size(0) % fine_gridptr->size(0) << std::endl;
  if (coarse_gridptr && check_partitioning && Dune::MPIHelper::getCollectiveCommunication().size() > 1) {
    const auto coarse_dimensions =
        DSG::dimensions(coarse_gridptr->leafGridView<PartitionIteratorType::Interior_Partition>());
    const auto fine_dimensions =
        DSG::dimensions(fine_gridptr->leafGridView<PartitionIteratorType::Interior_Partition>());
    for (const auto i : DSC::valueRange(world_dim)) {
      const bool match =
          DSC::FloatCmp::eq(coarse_dimensions.coord_limits[i].min(), fine_dimensions.coord_limits[i].min()) &&
          DSC::FloatCmp::eq(coarse_dimensions.coord_limits[i].max(), fine_dimensions.coord_limits[i].max());
      if (!match)
        DUNE_THROW(InvalidStateException, "Coarse and fine mesh do not match after load balancing, do \
                       you use different refinements in different spatial dimensions?");
    }
  }
  return fine_gridptr;
}
