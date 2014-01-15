#include <config.h>
#include <dune/stuff/common/parameter/configcontainer.hh>
#include <sstream>

#include "dune/multiscale/msfem/msfem_traits.hh"
#include "elliptic_msfem_matrix_assembler.hh"

namespace Dune {
namespace Multiscale {
namespace MsFEM {

DiscreteEllipticMsFEMOperator::DiscreteEllipticMsFEMOperator(
    MacroMicroGridSpecifierType& specifier, const CoarseDiscreteFunctionSpace& coarseDiscreteFunctionSpace,
    // number of layers per coarse grid entity T:  U(T) is created by enrichting T with
    // n(T)-layers:
    MsFEMTraits::LocalGridListType& subgrid_list, const DiffusionModel& diffusion_op)
  : specifier_(specifier)
  , coarseDiscreteFunctionSpace_(coarseDiscreteFunctionSpace)
  , subgrid_list_(subgrid_list)
  , diffusion_operator_(diffusion_op)
  , petrovGalerkin_(false) {}

} // namespace MsFEM {
} // namespace Multiscale {
} // namespace Dune {
