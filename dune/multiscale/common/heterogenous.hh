#ifndef DUNE_MULTISCALE_COMMON_HETEROGENOUS_HH
#define DUNE_MULTISCALE_COMMON_HETEROGENOUS_HH

#include <dune/common/deprecated.hh>
#include <dune/common/fvector.hh>

#include <dune/grid/common/backuprestore.hh>
#include <dune/grid/common/grid.hh>
#include <dune/grid/common/entity.hh>
#include <dune/grid/io/file/vtk/function.hh>

#include <dune/xt/common/ranges.hh>
#include <dune/stuff/aliases.hh>
#include <dune/xt/grid/search.hh>

#include <dune/gdt/spaces/cg/interface.hh>
#include <dune/multiscale/common/traits.hh>

namespace Dune {
namespace Multiscale {

class LocalsolutionProxy;
class LocalGridSearch;

template <class ImpTraits, size_t domainDim, size_t rangeDim>
std::vector<typename GDT::Spaces::CGInterface<ImpTraits, domainDim, rangeDim, 1>::DomainType> global_evaluation_points(
    const GDT::Spaces::CGInterface<ImpTraits, domainDim, rangeDim, 1>& space,
    const typename GDT::Spaces::CGInterface<ImpTraits, domainDim, rangeDim, 1>::EntityType& target_entity)
{
  const auto& target_lagrangepoint_set = space.lagrange_points(target_entity);
  const auto& target_geometry = target_entity.geometry();
  const auto quadNop = target_lagrangepoint_set.size();
  std::vector<typename GDT::Spaces::CGInterface<ImpTraits, domainDim, rangeDim, 1>::DomainType> points(quadNop);
  for (size_t qP = 0; qP < quadNop; ++qP) {
    points[qP] = target_geometry.global(target_lagrangepoint_set[qP]);
  }
  return points;
}

class MsFEMProjection
{
public:
  //! signature for non-default SearchStrategy constructions
  static void project(LocalsolutionProxy& source, CommonTraits::DiscreteFunctionType& target);

protected:
  static void preprocess(CommonTraits::DiscreteFunctionType& func);
  static void postprocess(CommonTraits::DiscreteFunctionType& func);
  static void identifySharedNodes(const CommonTraits::GridViewType& gridPart, std::vector<int>& map);
};

} // namespace Multiscale
} // namespace Dune

#endif // DUNE_MULTISCALE_COMMON_HETEROGENOUS_HH
