#include <config.h>

#include "proxygridview.hh"
#include <dune/multiscale/msfem/localproblems/localgridlist.hh>

namespace Dune {
namespace Multiscale {

ProxyGridview::ProxyGridview(const LocalGridList& localGrids)
  : BaseType(localGrids.getSubGrid(0))
  , localGrids_(localGrids)
{
  localGrids_.size();
}

} // namespace Dune
} // namespace Multiscale
