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
    MsFEMTraits::SubGridListType& subgrid_list, const DiffusionModel& diffusion_op)
  : specifier_(specifier)
  , coarseDiscreteFunctionSpace_(coarseDiscreteFunctionSpace)
  , subgrid_list_(subgrid_list)
  , diffusion_operator_(diffusion_op)
  , petrovGalerkin_(DSC_CONFIG_GET("msfem.petrov_galerkin", true)) {
  // coarseDiscreteFunctionSpace_ = specifier_.coarseSpace();
  // fineDiscreteFunctionSpace_ = specifier_.fineSpace();
  MsFEMLocalProblemSolverType(specifier_.fineSpace(), specifier_, subgrid_list_, diffusion_operator_).solve_all();
}

} // namespace MsFEM {
} // namespace Multiscale {
} // namespace Dune {
