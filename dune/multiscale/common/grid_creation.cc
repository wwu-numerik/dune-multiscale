#include <config.h>

#include "grid_creation.hh"

#include <dune/grid/common/gridenums.hh>
#include <dune/xt/common/ranges.hh>
#include <dune/stuff/grid/structuredgridfactory.hh>
#include <dune/stuff/grid/information.hh>
#include <dune/multiscale/problems/selector.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/multiscale/common/mygridfactory.hh>

using namespace Dune::Multiscale;
using namespace std;

typedef tuple<CommonTraits::DomainType, CommonTraits::DomainType, array<unsigned int, CommonTraits::world_dim>,
              array<unsigned int, CommonTraits::world_dim>,
              array<unsigned int, CommonTraits::world_dim>> SetupReturnType;

SetupReturnType setup(const DMP::ProblemContainer& problem) {
  BOOST_ASSERT_MSG(DXTC_CONFIG.has_sub("grids"), "Parameter tree needs to have 'grids' subtree!");

  const auto world_dim = CommonTraits::world_dim;
  typedef CommonTraits::DomainType CoordType;
  const auto& gridCorners = problem.getModelData().gridCorners();
  CoordType lowerLeft = gridCorners.first;
  CoordType upperRight = gridCorners.second;

  const auto oversamplingLayers = problem.config().get("msfem.oversampling_layers", 0);
  const auto validator = Dune::XT::Common::ValidateGreater<CommonTraits::DomainType>(CommonTraits::DomainType(1));
  const auto coarse_cells = problem.config().get<CommonTraits::DomainType>("grids.macro_cells_per_dim", world_dim, 0, validator);
  const auto microPerMacro = problem.config().get<CommonTraits::DomainType>("grids.micro_cells_per_macrocell_dim", world_dim, 0, validator);


  array<unsigned int, world_dim> elements, overCoarse, overFine;

  for (const auto i : Dune::XT::Common::value_range(world_dim)) {
    elements[i] = coarse_cells[i];
    overCoarse[i] = std::ceil(double(oversamplingLayers) / double(microPerMacro[i]));
    overFine[i] = problem.config().get("grids.overlap", 1);
  }
  return std::make_tuple(lowerLeft, upperRight, elements, overCoarse, overFine);
}

std::shared_ptr<CommonTraits::GridType>
Dune::Multiscale::make_coarse_grid(const DMP::ProblemContainer& problem,
                                   Dune::MPIHelper::MPICommunicator communicator) {
  CommonTraits::DomainType lowerLeft, upperRight;
  array<unsigned int, CommonTraits::world_dim> elements, overCoarse;
  std::tie(lowerLeft, upperRight, elements, overCoarse, std::ignore) = setup(problem);
  auto coarse_gridptr =
      MyGridFactory<CommonTraits::GridType>::createCubeGrid(lowerLeft, upperRight, elements, overCoarse, communicator);
  const auto expected_elements = std::accumulate(elements.begin(), elements.end(), 1u, std::multiplies<unsigned int>());
  auto actual_elements = coarse_gridptr->size(0) - coarse_gridptr->overlapSize(0);
  const auto sum_elements = coarse_gridptr->comm().sum(actual_elements);
  if (int(expected_elements) != sum_elements)
    DUNE_THROW(InvalidStateException, "Wonky grid distribution");
  if ((coarse_gridptr->comm().size() > 1) && (actual_elements == int(expected_elements)))
    DUNE_THROW(InvalidStateException, "Rank 0 fail");
  return coarse_gridptr;
}

pair<shared_ptr<CommonTraits::GridType>, shared_ptr<CommonTraits::GridType>>
Dune::Multiscale::make_grids(const DMP::ProblemContainer& problem, const bool check_partitioning,
                             Dune::MPIHelper::MPICommunicator communicator) {
  auto coarse_grid = make_coarse_grid(problem, communicator);
  return {coarse_grid, make_fine_grid(problem, coarse_grid, check_partitioning)};
}

std::shared_ptr<Dune::Multiscale::CommonTraits::GridType>
Dune::Multiscale::make_fine_grid(const DMP::ProblemContainer& problem,
                                 std::shared_ptr<Dune::Multiscale::CommonTraits::GridType> coarse_gridptr,
                                 const bool check_partitioning, Dune::MPIHelper::MPICommunicator communicator) {
  const auto world_dim = CommonTraits::world_dim;
  CommonTraits::DomainType lowerLeft, upperRight;
  array<unsigned int, world_dim> elements, overFine;
  std::tie(lowerLeft, upperRight, elements, std::ignore, overFine) = setup(problem);
  const auto validator = Dune::XT::Common::ValidateGreater<CommonTraits::DomainType>(CommonTraits::DomainType(1));
  const auto coarse_cells = problem.config().get<CommonTraits::DomainType>("grids.macro_cells_per_dim", world_dim, 0, validator);
  const auto microPerMacro = problem.config().get<CommonTraits::DomainType>("grids.micro_cells_per_macrocell_dim", world_dim, 0, validator);

  for (const auto i : Dune::XT::Common::value_range(CommonTraits::world_dim)) {
    elements[i] = coarse_cells[i] * microPerMacro[i];
  }
  auto fine_gridptr =
      MyGridFactory<CommonTraits::GridType>::createCubeGrid(lowerLeft, upperRight, elements, overFine, communicator);

  if (coarse_gridptr && check_partitioning && Dune::MPIHelper::getCollectiveCommunication().size() > 1) {
    // check whether grids match (may not match after load balancing if different refinements in different
    // spatial directions are used)
    MS_LOG_DEBUG << boost::format("Rank %d has %d coarse codim-0 elements and %d fine ones\n") %
                        coarse_gridptr->comm().rank() % coarse_gridptr->size(0) % fine_gridptr->size(0)
                 << std::endl;
    const auto fine_view = fine_gridptr->leafGridView<CommonTraits::InteriorPartition>();
    const auto coarse_view = coarse_gridptr->leafGridView<CommonTraits::InteriorPartition>();
    //    if(coarse_view.size(0) != std::pow(coarse_cells[0], CommonTraits::world_dim))
    //      DUNE_THROW(InvalidStateException, "snafu " << std::pow(coarse_cells[0], CommonTraits::world_dim)
    //          << " | " << coarse_view.size(0) << '\n');
    const auto coarse_dimensions = DSG::dimensions(coarse_view);
    const auto fine_dimensions = DSG::dimensions(fine_view);
    for (const auto i : Dune::XT::Common::value_range(world_dim)) {
      const bool match =
          Dune::XT::Common::FloatCmp::eq(coarse_dimensions.coord_limits[i].min(), fine_dimensions.coord_limits[i].min()) &&
          Dune::XT::Common::FloatCmp::eq(coarse_dimensions.coord_limits[i].max(), fine_dimensions.coord_limits[i].max());
      if (!match)
        DUNE_THROW(InvalidStateException, "Coarse and fine mesh do not match after load balancing, do \
                       you use different refinements in different spatial dimensions?");
    }
  }
  return fine_gridptr;
}
