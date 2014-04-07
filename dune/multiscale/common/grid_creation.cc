#include <config.h>
#include "grid_creation.hh"

#include <dune/stuff/common/ranges.hh>
#include <dune/stuff/grid/structuredgridfactory.hh>
#include <dune/stuff/grid/information.hh>

std::pair<std::shared_ptr<Dune::Multiscale::CommonTraits::GridType>,
          std::shared_ptr<Dune::Multiscale::CommonTraits::GridType>>
Dune::Multiscale::make_grids() {
  const int dim_world = CommonTraits::GridType::dimensionworld;
  typedef FieldVector<typename CommonTraits::GridType::ctype, dim_world> CoordType;
  CoordType lowerLeft(0.0);
  CoordType upperRight(1.0);

  const auto coarse_cells = DSC_CONFIG_GET("global.macro_cells_per_dim", 8);
  array<unsigned int, dim_world> elements;
  for (const auto i : DSC::valueRange(dim_world))
    elements[i] = coarse_cells;
  auto coarse_gridptr =
      StructuredGridFactory<CommonTraits::GridType>::createCubeGrid(lowerLeft, upperRight, elements, overCoarse);

  for (const auto i : DSC::valueRange(dim_world)) {
    elements[i] = coarse_cells * DSC_CONFIG_GET("global.micro_cells_per_macrocell_dim", 8);
  }
  auto fine_gridptr = StructuredGridFactory<CommonTraits::GridType>::createCubeGrid(lowerLeft, upperRight, elements);
  return {coarse_gridptr, fine_gridptr};
}
