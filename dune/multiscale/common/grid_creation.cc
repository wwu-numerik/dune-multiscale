#include <config.h>
#include "grid_creation.hh"

#include <dune/stuff/common/ranges.hh>
#include <dune/stuff/grid/structuredgridfactory.hh>
#include <dune/stuff/grid/information.hh>
#include <dune/stuff/common/parameter/tree.hh>
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
  BOOST_ASSERT_MSG(DSC_CONFIG.hasSub("grids"), "Parameter tree needs to have 'grids' subtree!");
  
  const DSC::ExtendedParameterTree gridParameterTree(DSC_CONFIG.sub("grids"));
  const int dim_world = CommonTraits::GridType::dimensionworld;
  typedef FieldVector<typename CommonTraits::GridType::ctype, dim_world> CoordType;
  const auto& gridCorners = Problem::getModelData()->gridCorners();
  CoordType lowerLeft = gridCorners.first;
  CoordType upperRight = gridCorners.second;

  const auto oversamplingLayers = DSC_CONFIG_GET("msfem.oversampling_layers", 0);
  std::vector<int> microPerMacro = gridParameterTree.getVector("micro_cells_per_macrocell_dim", 8, dim_world);

  const std::vector<int> coarse_cells = gridParameterTree.getVector("macro_cells_per_dim", 8, dim_world);
  array<unsigned int, dim_world> elements;
  array<unsigned int, dim_world> overCoarse;
  for (const auto i : DSC::valueRange(dim_world)) {
    elements[i] = coarse_cells[i];
    overCoarse[i] = std::ceil(double(oversamplingLayers)/double(microPerMacro[i]));
  }
  auto coarse_gridptr =
      MyGridFactory<CommonTraits::GridType>::createCubeGrid(lowerLeft, upperRight, elements, overCoarse);


  for (const auto i : DSC::valueRange(dim_world)) {
    elements[i] = coarse_cells[i] * microPerMacro[i];
  }
  auto fine_gridptr = StructuredGridFactory<CommonTraits::GridType>::createCubeGrid(lowerLeft, upperRight, elements);

  // check whether grids match (may not match after load balancing if different refinements in different
  // spatial directions are used)
  if (Dune::Fem::MPIManager::size()>1) {
    const auto coarse_dimensions = DSG::dimensions<CommonTraits::GridType>(*coarse_gridptr);
    const auto fine_dimensions = DSG::dimensions<CommonTraits::GridType>(*fine_gridptr);
    const auto eps = std::numeric_limits<decltype(coarse_dimensions.coord_limits[0].max())>::epsilon();
    for (const auto i : DSC::valueRange(dim_world)) {
      const bool match = (std::abs(coarse_dimensions.coord_limits[i].min()-fine_dimensions.coord_limits[i].min())
                          < eps)
                          && (std::abs(coarse_dimensions.coord_limits[i].max()-fine_dimensions.coord_limits[i].max())
                              < eps);
      if (!match)
        DUNE_THROW(InvalidStateException, "Coarse and fine mesh do not match after load balancing, do \
                       you use different refinements in different spatial dimensions?");
    }
  }

  return {coarse_gridptr, fine_gridptr};
}
