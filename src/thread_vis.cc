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

void iterate(std::string macroGridName, std::string function_name, int level, int threadnum = 8)
{
  Dune::Fem::ThreadManager::setMaxNumberThreads(threadnum);
  typedef Dune::Fem::FiniteVolumeSpace<CommonTraits::FunctionSpaceType, CommonTraits::GridPartType, 0> FVSpace;
  typedef Dune::Fem::AdaptiveDiscreteFunction<FVSpace> FVFunc;

  CommonTraits::GridPointerType macro_grid_pointer(macroGridName);
  // refine the grid 'starting_refinement_level' times:
  macro_grid_pointer->globalRefine(level);
  CommonTraits::GridType& grid = *macro_grid_pointer;
  CommonTraits::GridPartType gridPart(grid);
  FVSpace fv_space(gridPart);
  CommonTraits::DiscreteFunctionSpaceType discreteFunctionSpace(gridPart);

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

//  OutputParameters outputparam;

//  outputparam.set_prefix(function_name);
//  Dune::Fem::DataWriter<CommonTraits::GridType, decltype(functions)> output(gridPart.grid(), functions, outputparam);
//  output.writeData(1.0 /*dummy*/, function_name);

  Dune::Fem::VTKIO<CommonTraits::GridPartType> vtkio(gridPart);

  for(auto& f : functions)
  {
    vtkio.addCellData(*f,f->name());
  }
  vtkio.write("all");
}

void algorithm(const std::string& macroGridName, int total_refinement_level_,
               int coarse_grid_level_, int number_of_layers_ ) {

  DSC_LOG_INFO_0 << "loading dgf: " << macroGridName << std::endl;

  iterate(macroGridName, "fine_grid", total_refinement_level_);
  iterate(macroGridName, "fine_grid", coarse_grid_level_);


//  const auto number_of_level_host_entities = grid_coarse.size(0 /*codim*/);

//  // number of layers per coarse grid entity T:  U(T) is created by enrichting T with n(T)-layers.
//  MsFEM::MsFEMTraits::MacroMicroGridSpecifierType specifier(discreteFunctionSpace_coarse, discreteFunctionSpace);
//  for (int i = 0; i < number_of_level_host_entities; ++i) {
//    specifier.setNoOfLayers(i, number_of_layers_);
//  }
//  specifier.setOversamplingStrategy(DSC_CONFIG_GET("msfem.oversampling_strategy", 1));

//  //! create subgrids:
//  { // this scopes subgridlist
//    MsFEM::MsFEMTraits::SubGridListType subgrid_list(specifier, DSC_CONFIG_GET("logging.subgrid_silent", false));

//  }



} // function algorithm


} // namespace Dune {
} // namespace Multiscale {

int main(int argc, char** argv) {
  try {
    using namespace Dune::Multiscale;
    using namespace Dune::Multiscale::MsFEM;
    init(argc, argv);

    //!TODO include base in config
    DSC_PROFILER.startTiming("msfem.all");

    const std::string datadir = DSC_CONFIG_GET("global.datadir", "data/");

    // generate directories for data output
    DSC::testCreateDirectory(datadir);
    DSC_LOG_INFO_0 << boost::format("Data will be saved under: %s\nLogs will be saved under: %s/%s/ms.log.log\n") %
                          datadir % datadir % DSC_CONFIG_GET("logging.dir", "log");

    // syntax: info_from_par_file / default  / validation of the value

    // coarse_grid_level denotes the (starting) grid refinement level for the global coarse scale problem, i.e. it
    // describes 'H'
    int coarse_grid_level_ = DSC_CONFIG_GETV("msfem.coarse_grid_level", 4, DSC::ValidateLess<int>(-1));

    // syntax: info_from_par_file / default
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
    algorithm(macroGridName, total_refinement_level_, coarse_grid_level_, number_of_layers_);
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
