#ifndef DUNE_MULTISCALE_PROXYGRIDVIEW_HH
#define DUNE_MULTISCALE_PROXYGRIDVIEW_HH

#include <dune/multiscale/msfem/msfem_traits.hh>
#include <dune/grid/common/gridview.hh>

namespace Dune {
namespace Multiscale {

namespace MsFEM {
class LocalGridList;
}

struct ProxyGridviewTraits : public DefaultLeafGridViewTraits<MsFEM::MsFEMTraits::LocalGridType, All_Partition> {};

class ProxyGridview : public GridView<ProxyGridviewTraits> {
  typedef GridView<ProxyGridviewTraits> BaseType;

public:
  typedef MsFEM::MsFEMTraits::LocalGridType::ctype ctype;
  static auto constexpr dimension = CommonTraits::world_dim;
  typedef MsFEM::MsFEMTraits::LocalEntityType EntityType;

  template <int cd>
  using Codim = MsFEM::MsFEMTraits::LocalGridType::Codim<cd>;

  ProxyGridview(const MsFEM::LocalGridList& localGrids);

private:
  const MsFEM::LocalGridList& localGrids_;
};

} // namespace Dune
} // namespace Multiscale

#endif // DUNE_MULTISCALE_PROXYGRIDVIEW_HH
