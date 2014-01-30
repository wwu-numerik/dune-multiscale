#include <config.h>
// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include <dune/multiscale/common/main_init.hh>
#include <dune/multiscale/common/traits.hh>
#include <dune/multiscale/common/grid_creation.hh>
#include <dune/multiscale/msfem/msfem_traits.hh>
#include <dune/multiscale/problems/selector.hh>
#include <dune/multiscale/msfem/localproblems/subgrid-list.hh>
#include <dune/multiscale/tools/misc/outputparameter.hh>
#include <dune/multiscale/msfem/localproblems/localsolutionmanager.hh>

#include <dune/fem/space/finitevolume.hh>
#include <dune/fem/space/lagrange.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/misc/threads/threadmanager.hh>
#include <dune/fem/misc/threads/domainthreaditerator.hh>
#include <dune/fem/io/file/vtkio.hh>

#include <dune/stuff/grid/information.hh>
#include <dune/stuff/grid/structuredgridfactory.hh>
#include <dune/stuff/common/ranges.hh>
#include <dune/stuff/discretefunction/projection/heterogenous.hh>


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

void partition_vis_single(const std::shared_ptr<CommonTraits::GridType>& grid, std::string function_name, int threadnum = 8)
{
  Dune::Fem::ThreadManager::setMaxNumberThreads(threadnum);
  //Dune::Fem::Parameter::replace(std::string("fem.threads.partitioningmethod"), std::string("kway"));
  typedef Dune::Fem::FiniteVolumeSpace<CommonTraits::FunctionSpaceType, CommonTraits::GridPartType, 0> FVSpace;
  typedef Dune::Fem::AdaptiveDiscreteFunction<FVSpace> FVFunc;
  CommonTraits::GridPartType gridPart(*grid);
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

void subgrid_vis()
{
  auto grids = make_grids();
  CommonTraits::GridType& coarse_grid = *grids.first;
  CommonTraits::GridPartType coarse_gridPart(coarse_grid);
  CommonTraits::GridType& fine_grid = *grids.second;
  CommonTraits::GridPartType fine_gridPart(fine_grid);

  CommonTraits::DiscreteFunctionSpaceType coarse_space(coarse_gridPart);
  MsFEM::LocalGridList subgrid_list(coarse_space);
  typedef Dune::Fem::FiniteVolumeSpace<CommonTraits::FunctionSpaceType, CommonTraits::GridPartType, 0> FVSpace;
  typedef Dune::Fem::FiniteVolumeSpace<CommonTraits::FunctionSpaceType, DMM::MsFEMTraits::LocalGridPartType, 0> LocalFVSpace;
  typedef Dune::Fem::AdaptiveDiscreteFunction<FVSpace> FVFunc;
  typedef Dune::Fem::AdaptiveDiscreteFunction<LocalFVSpace> LocalFVFunc;
  FVSpace fv_space(fine_gridPart);

  std::vector<std::unique_ptr<FVFunc>> oversampled_functions(subgrid_list.size());
  std::vector<std::unique_ptr<FVFunc>> restricted_functions(subgrid_list.size());
  std::vector<std::unique_ptr<FVFunc>> weak_restricted_functions(subgrid_list.size());

  auto oversampled_function_it = oversampled_functions.begin();
  auto restricted_function_it = restricted_functions.begin();
  auto weak_restricted_function_it = weak_restricted_functions.begin();

  // horrible, horrible complexity :)
  for(const auto& coarse_entity : coarse_space)
  {
    DMM::LocalSolutionManager localSolManager(coarse_space, coarse_entity, subgrid_list);
    const auto& local_space = localSolManager.space();
    const auto& id_set = coarse_space.gridPart().grid().leafIndexSet();
    const auto coarse_id = id_set.index(coarse_entity);
    LocalFVSpace local_fv_space(local_space.gridPart());
    LocalFVFunc local_oversampled("", local_fv_space);
    LocalFVFunc local_restricted("", local_fv_space);
    LocalFVFunc local_weak_restricted("", local_fv_space);

    for(const auto& local_entity : local_space)
    {
      auto local_oversampled_lf = local_oversampled.localFunction(local_entity);
      for (const auto idx : DSC::valueRange(local_oversampled_lf.size())) {
          local_oversampled_lf[idx] = static_cast<unsigned long>(coarse_id+1);
      }
      if(subgrid_list.covers_strict(coarse_entity, local_entity)) {
        auto local_restricted_lf = local_restricted.localFunction(local_entity);
        for (const auto idx : DSC::valueRange(local_restricted_lf.size())) {
            local_restricted_lf[idx] = coarse_id+1;
        }
      }
      if(subgrid_list.covers(coarse_entity, local_entity)) {
        auto local_weak_restricted_lf = local_weak_restricted.localFunction(local_entity);
        for (const auto idx : DSC::valueRange(local_weak_restricted_lf.size())) {
            local_weak_restricted_lf[idx] = coarse_id+1;
        }
      }
    }

    auto& oversampled_function = (*oversampled_function_it++);
    oversampled_function = DSC::make_unique<FVFunc>(DSC::toString(coarse_id) + "_oversampled", fv_space);
    oversampled_function->clear();
    auto& restricted_function = (*restricted_function_it++);
    restricted_function = DSC::make_unique<FVFunc>(DSC::toString(coarse_id) + "_restricted", fv_space);
    restricted_function->clear();
    auto& weak_restricted_function = (*weak_restricted_function_it++);
    weak_restricted_function = DSC::make_unique<FVFunc>(DSC::toString(coarse_id) + "_weak_restricted", fv_space);
    weak_restricted_function->clear();
    DS::HeterogenousProjection<>::project(local_oversampled, *oversampled_function);
    DS::HeterogenousProjection<>::project(local_restricted, *restricted_function);
    DS::HeterogenousProjection<>::project(local_weak_restricted, *weak_restricted_function);
  }

  output_all(oversampled_functions, fine_gridPart, "oversampled_");
  output_all(restricted_functions, fine_gridPart, "restricted_");
  output_all(weak_restricted_functions, fine_gridPart, "weak_restricted_");
}

void partition_vis() {
  auto grids = make_grids();
  partition_vis_single(grids.second, "fine_grid");
  partition_vis_single(grids.first, "coarse_grid");
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

    partition_vis();
    subgrid_vis();
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
