#include <config.h>
#include "grid_creation.hh"

#include <dune/stuff/common/ranges.hh>
#include <dune/stuff/grid/structuredgridfactory.hh>
#include <dune/stuff/grid/information.hh>

std::pair<std::shared_ptr<Dune::Multiscale::CommonTraits::GridType>, std::shared_ptr<Dune::Multiscale::CommonTraits::GridType> > Dune::Multiscale::make_grids() {
  const int dim_world = CommonTraits::GridType::dimensionworld;
  typedef FieldVector<typename CommonTraits::GridType::ctype, dim_world> CoordType;
  const auto coarse_cells = DSC_CONFIG_GET("global.macro_cells_per_dim", 8);
  array<unsigned int,dim_world> elements;
  for(const auto i : DSC::valueRange(dim_world))
    elements[i] = coarse_cells;
  auto coarse_gridptr = StructuredGridFactory<CommonTraits::GridType>
                        ::createCubeGrid(CoordType(0.0), CoordType(1.0), elements);

  const auto coarse_dimensions = DSG::dimensions(*coarse_gridptr);
  CoordType lowerLeft(0);
  CoordType upperRight(0);
  for(const auto i : DSC::valueRange(dim_world))
  {
    elements[i] = coarse_cells * DSC_CONFIG_GET("global.micro_cells_per_macrocell_dim", 8);
    lowerLeft[i] = coarse_dimensions.coord_limits[i].min();
    upperRight[i] = coarse_dimensions.coord_limits[i].max();
  }
  auto fine_gridptr = StructuredGridFactory<CommonTraits::GridType>::createCubeGrid(lowerLeft, upperRight, elements);
  return {coarse_gridptr,fine_gridptr};
}
