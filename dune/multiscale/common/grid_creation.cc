#include <config.h>
#include "grid_creation.hh"

#include <dune/stuff/common/ranges.hh>
#include <dune/stuff/grid/structuredgridfactory.hh>
#include <dune/stuff/grid/information.hh>
#include <dune/multiscale/problems/selector.hh>

// Helper struct to make overlap for SPGrid possible
// Declared in unnamed namespace to avoid naming conflicts
namespace {
template< class GridType >
class MyGridFactory {
  typedef typename GridType::ctype ctype;
  static const int dimworld = GridType::dimensionworld;
  static const int dim = GridType::dimension;

public:
  static std::shared_ptr<GridType>
  createCubeGrid(const Dune::FieldVector<ctype,dimworld>& lowerLeft,
                 const Dune::FieldVector<ctype,dimworld>& upperRight,
                 const Dune::array<unsigned int,dim>& elements,
                 const Dune::array<unsigned int,dim>& /*overlap*/) {
    // structured grid factory allows overlap only for SPGrid at the moment, hence the following check
    BOOST_ASSERT_MSG((DSC_CONFIG_GET("msfem.oversampling_layers", 0)==0),
                     "Oversampling may only be used in combination with SPGrid!");
    return Dune::StructuredGridFactory<GridType>::createCubeGrid(lowerLeft, upperRight, elements);
  }
};

template< class ct, int dim, Dune::SPRefinementStrategy strategy, class Comm >
class MyGridFactory< Dune::SPGrid< ct, dim, strategy, Comm > > {
  typedef Dune::SPGrid< ct, dim, strategy, Comm > GridType;
  typedef typename GridType::ctype ctype;
  static const int dimworld = GridType::dimensionworld;

public:
  static std::shared_ptr<GridType>
  createCubeGrid(const Dune::FieldVector<ctype,dimworld>& lowerLeft,
                 const Dune::FieldVector<ctype,dimworld>& upperRight,
                 const Dune::array<unsigned int,dim>& elements,
                 const Dune::array<unsigned int,dim>& overlap) {
    return Dune::StructuredGridFactory<GridType>::createCubeGrid(lowerLeft, upperRight, elements, overlap);
  }
};
}

std::pair<std::shared_ptr<Dune::Multiscale::CommonTraits::GridType>,
          std::shared_ptr<Dune::Multiscale::CommonTraits::GridType>>
Dune::Multiscale::make_grids() {
  const int dim_world = CommonTraits::GridType::dimensionworld;
  typedef FieldVector<typename CommonTraits::GridType::ctype, dim_world> CoordType;
  const auto& gridCorners = Problem::getModelData()->gridCorners();
  CoordType lowerLeft = gridCorners.first;
  CoordType upperRight = gridCorners.second;

  const auto oversamplingLayers = DSC_CONFIG_GET("msfem.oversampling_layers", 0);
  const auto microPerMacro = DSC_CONFIG_GET("msfem.micro_cells_per_macrocell_dim", 8);
  const int overlapLayers = std::ceil(double(oversamplingLayers)/double(microPerMacro));

  const auto coarse_cells = DSC_CONFIG_GET("global.macro_cells_per_dim", 8);
  array<unsigned int, dim_world> elements;
  array<unsigned int, dim_world> overCoarse;
  for (const auto i : DSC::valueRange(dim_world)) {
    elements[i] = coarse_cells;
    overCoarse[i] = overlapLayers;
  }
  auto coarse_gridptr =
      MyGridFactory<CommonTraits::GridType>::createCubeGrid(lowerLeft, upperRight, elements, overCoarse);

  for (const auto i : DSC::valueRange(dim_world)) {
    elements[i] = coarse_cells * microPerMacro;
  }
  auto fine_gridptr = StructuredGridFactory<CommonTraits::GridType>::createCubeGrid(lowerLeft, upperRight, elements);
  return {coarse_gridptr, fine_gridptr};
}
