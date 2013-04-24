#include "elliptic_msfem_matrix_assembler.hh"

namespace Dune {
namespace Multiscale {
namespace MsFEM {

DiscreteEllipticMsFEMOperator::DiscreteEllipticMsFEMOperator(MacroMicroGridSpecifierType& specifier,
                              const CoarseDiscreteFunctionSpace& coarseDiscreteFunctionSpace,
                              // number of layers per coarse grid entity T:  U(T) is created by enrichting T with
                              // n(T)-layers:
                              MsFEMTraits::SubGridListType& subgrid_list,
                              const DiffusionModel& diffusion_op)
  : specifier_(specifier)
    , coarseDiscreteFunctionSpace_(coarseDiscreteFunctionSpace)
    , subgrid_list_(subgrid_list)
    , diffusion_operator_(diffusion_op) {
  bool silence = false;

  // coarseDiscreteFunctionSpace_ = specifier_.coarseSpace();
  // fineDiscreteFunctionSpace_ = specifier_.fineSpace();
  MsFEMLocalProblemSolverType loc_prob_solver(
    specifier_.fineSpace(), specifier_, subgrid_list_, diffusion_operator_);

  loc_prob_solver.assemble_all(silence);
}

void DiscreteEllipticMsFEMOperator::subgrid_to_hostrid_function(const LocalDiscreteFunction& sub_func,
                                                                                FineDiscreteFunction& host_func) {
  if ( sub_func.space().gridPart().grid().maxLevel() != host_func.space().gridPart().grid().maxLevel() )
  {
    DSC_LOG_ERROR
    << "Error in method 'subgrid_to_hostrid_function': MaxLevel of SubGrid not identical to MaxLevel of FineGrid."
    << std::endl;
  }

  host_func.clear();

  const LocalDiscreteFunctionSpace& subDiscreteFunctionSpace = sub_func.space();
  const auto& subGrid = subDiscreteFunctionSpace.grid();

  LocalGridIterator sub_endit = subDiscreteFunctionSpace.end();
  for (LocalGridIterator sub_it = subDiscreteFunctionSpace.begin(); sub_it != sub_endit; ++sub_it)
  {
    const LocalGridEntity& sub_entity = *sub_it;

    auto host_entity_pointer = subGrid.getHostEntity< 0 >(*sub_it);
    const FineEntity& host_entity = *host_entity_pointer;

    LocalGridLocalFunction sub_loc_value = sub_func.localFunction(sub_entity);
    FineLocalFunction host_loc_value = host_func.localFunction(host_entity);

    const auto numBaseFunctions = sub_loc_value.baseFunctionSet().size();
    for (unsigned int i = 0; i < numBaseFunctions; ++i)
    {
      host_loc_value[i] = sub_loc_value[i];
    }
  }
} // subgrid_to_hostrid_function



} //namespace MsFEM {
} //namespace Multiscale {
} //namespace Dune {

