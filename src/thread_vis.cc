// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include <dune/multiscale/common/main_init.hh>
#include <dune/multiscale/common/traits.hh>
#include <dune/multiscale/msfem/msfem_traits.hh>
#include <dune/multiscale/problems/selector.hh>
#include <dune/multiscale/msfem/msfem_grid_specifier.hh>
#include <dune/multiscale/msfem/localproblems/subgrid-list.hh>
#include <dune/multiscale/tools/misc/outputparameter.hh>

#include <dune/fem/space/finitevolume.hh>
#include <dune/fem/space/lagrange.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/misc/threads/threadmanager.hh>
#include <dune/fem/misc/threads/domainthreaditerator.hh>
#include <dune/fem/io/file/vtkio.hh>

#include <dune/stuff/common/ranges.hh>

namespace Dune {
namespace Multiscale {

template <class FunctionType>
void output_all(std::vector<std::unique_ptr<FunctionType>>& functions, CommonTraits::GridPartType& gridpart, std::string name )
{
  Dune::Fem::VTKIO<CommonTraits::GridPartType> vtkio(gridpart);

  for(auto& f : functions)
  {
    vtkio.addCellData(*f,f->name());
  }
  boost::filesystem::path out_filename(DSC_CONFIG_GET("global.datadir", "data"));
  out_filename /= name;
  vtkio.write(out_filename.string());
}

void partition_vis_single(std::string macroGridName, std::string function_name, int level, int threadnum = 8)
{
  Dune::Fem::ThreadManager::setMaxNumberThreads(threadnum);
  //Dune::Fem::Parameter::replace(std::string("fem.threads.partitioningmethod"), std::string("kway"));
  typedef Dune::Fem::FiniteVolumeSpace<CommonTraits::FunctionSpaceType, CommonTraits::GridPartType, 0> FVSpace;
  typedef Dune::Fem::AdaptiveDiscreteFunction<FVSpace> FVFunc;

  CommonTraits::GridPointerType macro_grid_pointer(macroGridName);
  // refine the grid 'starting_refinement_level' times:
  Dune::Fem::GlobalRefine::apply(*macro_grid_pointer, level);
  CommonTraits::GridType& grid = *macro_grid_pointer;
  CommonTraits::GridPartType gridPart(grid);
  FVSpace fv_space(gridPart);

  std::vector<std::unique_ptr<FVFunc>> functions(threadnum+1);

  Dune::Fem::DomainDecomposedIteratorStorage< CommonTraits::GridPartType > iterators(gridPart);
  iterators.update();
  for (auto thread : DSC::valueRange(threadnum))
  {
    const auto fn = function_name + "_thread_" + DSC::toString(thread);
    functions[thread] = DSC::make_unique<FVFunc>(fn, fv_space);
  }
  auto& combined = functions.back();
  combined = DSC::make_unique<FVFunc>(function_name + "_combined", fv_space);
  combined->clear();

  #ifdef _OPENMP
  #pragma omp parallel
  #endif
  {
    const auto thread =  Dune::Fem::ThreadManager::thread();
    auto& function = functions[thread];
    function->clear();
    for(const auto& entity : iterators)
    {
      auto local_function = function->localFunction(entity);
      auto local_combined = combined->localFunction(entity);
      for (const auto idx : DSC::valueRange(local_function.size())) {
        local_function[idx] = thread+1;
        local_combined[idx] = thread+1;
      }
    }
  }

  output_all(functions, gridPart, function_name+"all");
}

void subgrid_vis(const std::string& macroGridName, int total_refinement_level_,
                 int coarse_grid_level_, int number_of_layers_ )
{
  CommonTraits::GridPointerType macro_grid_pointer(macroGridName);
  // refine the grid 'starting_refinement_level' times:
  Dune::Fem::GlobalRefine::apply(*macro_grid_pointer, coarse_grid_level_);

  CommonTraits::GridType& grid = *macro_grid_pointer;
  CommonTraits::GridPartType gridPart(grid);
  // coarse grid
  CommonTraits::GridPointerType macro_grid_pointer_coarse(macroGridName);
  CommonTraits::GridType& grid_coarse = *macro_grid_pointer_coarse;
  Dune::Fem::GlobalRefine::apply(grid_coarse, coarse_grid_level_);
  CommonTraits::GridPartType gridPart_coarse(grid_coarse);

  Dune::Fem::GlobalRefine::apply(grid, total_refinement_level_ - coarse_grid_level_);

  //! ------------------------- discrete function spaces -----------------------------------
  // the global-problem function space:
  CommonTraits::DiscreteFunctionSpaceType discreteFunctionSpace(gridPart);
  CommonTraits::DiscreteFunctionSpaceType discreteFunctionSpace_coarse(gridPart_coarse);
  const auto number_of_level_host_entities = grid_coarse.size(0 /*codim*/);

  // number of layers per coarse grid entity T:  U(T) is created by enrichting T with n(T)-layers.
  MsFEM::MsFEMTraits::MacroMicroGridSpecifierType specifier(discreteFunctionSpace_coarse, discreteFunctionSpace);
  for (int i = 0; i < number_of_level_host_entities; ++i) {
    specifier.setNoOfLayers(i, number_of_layers_);
  }
  specifier.setOversamplingStrategy(DSC_CONFIG_GET("msfem.oversampling_strategy", 1));
  MsFEM::MsFEMTraits::SubGridListType subgrid_list(specifier, DSC_CONFIG_GET("logging.subgrid_silent", false));

  typedef Dune::Fem::FiniteVolumeSpace<CommonTraits::FunctionSpaceType, CommonTraits::GridPartType, 0> FVSpace;
  typedef Dune::Fem::AdaptiveDiscreteFunction<FVSpace> FVFunc;
  FVSpace fv_space(gridPart);

  std::vector<std::unique_ptr<FVFunc>> oversampled_functions(subgrid_list.size());
  std::vector<std::unique_ptr<FVFunc>> functions(subgrid_list.size());


  auto oversampled_function_it = oversampled_functions.begin();
  auto function_it = functions.begin();
  // horrible, horrible complexity :)
  for(const auto& coarse_entity : discreteFunctionSpace_coarse)
  {
    const auto& subgrid = subgrid_list.getSubGrid(coarse_entity);
    const auto& id_set = discreteFunctionSpace_coarse.gridPart().grid().globalIdSet();
    const auto coarse_id = id_set.id(coarse_entity);
    auto& oversampled_function = (*oversampled_function_it++);
    oversampled_function = DSC::make_unique<FVFunc>(DSC::toString(coarse_id) + "_subgrid", fv_space);
    oversampled_function->clear();
    auto& function = (*function_it++);
    function = DSC::make_unique<FVFunc>(DSC::toString(coarse_id) + "_coarse_cell", fv_space);
    function->clear();

    for(const auto& fine_entity : discreteFunctionSpace)
    {
      if(subgrid.contains<0>(fine_entity)) {
        auto oversampled_local_function = oversampled_function->localFunction(fine_entity);
        for (const auto idx : DSC::valueRange(oversampled_local_function.size())) {
          oversampled_local_function[idx] = coarse_id+1;
        }
        if (coarse_id == subgrid_list.getEnclosingMacroCellId(CommonTraits::EntityPointerType(fine_entity)))
        {
          auto local_function = function->localFunction(fine_entity);
          for (const auto idx : DSC::valueRange(local_function.size())) {
            local_function[idx] = coarse_id+1;
          }
        }
      }
    }
  }

  output_all(oversampled_functions, gridPart, "subgrids");
  output_all(functions, gridPart, "coarse_cells");
}

void partition_vis(const std::string& macroGridName, int total_refinement_level_,
               int coarse_grid_level_) {
  partition_vis_single(macroGridName, "fine_grid", total_refinement_level_);
  partition_vis_single(macroGridName, "coarse_grid", coarse_grid_level_);
} // function algorithm


} // namespace Dune {
} // namespace Multiscale {

int main(int argc, char** argv) {
  try {
    using namespace Dune::Multiscale;
    using namespace Dune::Multiscale::MsFEM;
    init(argc, argv);

    const std::string datadir = DSC_CONFIG_GET("global.datadir", "data/");

    // generate directories for data output
    DSC::testCreateDirectory(datadir);
    DSC_LOG_INFO_0 << boost::format("Data will be saved under: %s\nLogs will be saved under: %s/%s/ms.log.log\n") %
                          datadir % datadir % DSC_CONFIG_GET("logging.dir", "log");
    int coarse_grid_level_ = DSC_CONFIG_GETV("msfem.coarse_grid_level", 4, DSC::ValidateLess<int>(-1));
    int number_of_layers_ = DSC_CONFIG_GET("msfem.oversampling_layers", 4);

    switch (DSC_CONFIG_GET("msfem.oversampling_strategy", 1)) {
      case 1:
        break;
      case 2:
        break;
      default:
        DUNE_THROW(Dune::InvalidStateException, "Oversampling Strategy must be 1 or 2.");
    }

    // data for the model problem; the information manager
    // (see 'problem_specification.hh' for details)
    auto info_ptr = Problem::getModelData();
    const auto& info = *info_ptr;

    // total_refinement_level denotes the (starting) grid refinement level for the global fine scale problem, i.e. it
    // describes 'h'
    int total_refinement_level_ =
        DSC_CONFIG_GETV("msfem.fine_grid_level", 4, DSC::ValidateLess<int>(coarse_grid_level_ - 1));

    // name of the grid file that describes the macro-grid:
    const std::string macroGridName = info.getMacroGridFile();
    partition_vis(macroGridName, total_refinement_level_, coarse_grid_level_);
    subgrid_vis(macroGridName, total_refinement_level_, coarse_grid_level_, number_of_layers_);
    return 0;
  }
  catch (Dune::Exception& e) {
    std::cerr << e.what() << std::endl;
  }
  catch (const std::exception& ex) {
    std::cerr << "Caught std::exception: " << ex.what() << "\n";
  }
  catch (...) {
    std::cerr << "Exception of non-known type thrown!\n";
  }
  return 1;
} // main
