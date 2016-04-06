#include <config.h>
// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include <dune/multiscale/common/main_init.hh>
#include <dune/multiscale/common/traits.hh>
#include <dune/multiscale/common/grid_creation.hh>
#include <dune/multiscale/msfem/msfem_traits.hh>
#include <dune/multiscale/problems/selector.hh>
#include <dune/multiscale/msfem/localproblems/localgridlist.hh>
#include <dune/multiscale/tools/misc/outputparameter.hh>

#include <dune/stuff/aliases.hh>
#include <dune/xt/common/ranges.hh>
#include <dune/xt/common/parallel/threadmanager.hh>
#include <dune/stuff/grid/output/entity_visualization.hh>

#include <dune/gdt/spaces/fv/default.hh>

namespace Dune {
namespace Multiscale {

template <class FunctionType>
void output_all(std::vector<std::unique_ptr<FunctionType>>& functions, CommonTraits::GridViewType& view,
                std::string name) {
  Dune::VTKWriter<CommonTraits::GridViewType> vtkio(view);

  for (auto& f : functions) {
    auto adapter = std::make_shared<DSFu::VisualizationAdapter<CommonTraits::GridViewType, 1, 1>>(*f);
    vtkio.addCellData(adapter);
  }
  const std::string datadir = DXTC_CONFIG_GET("global.datadir", "./data/");
  Dune::XT::Common::test_create_directory(datadir + "/piecefiles/" + name);
  vtkio.pwrite(name, datadir, "piecefiles");
}

void partition_vis_single(CommonTraits::GridType& grid, std::string function_name) {
  const auto threadnum = Dune::XT::CommonthreadManager().max_threads();
  // Dune::Fem::Parameter::replace(std::string("fem.threads.partitioningmethod"), std::string("kway"));
  auto view_ptr = grid.leafGridView();
  typedef GDT::Spaces::FV::Default<CommonTraits::GridViewType, double, 1, 1> FVSpace;
  typedef BackendChooser<FVSpace>::DiscreteFunctionType FVFunc;
  FVSpace fv_space(view_ptr);

  std::vector<std::unique_ptr<FVFunc>> functions(threadnum + 1);

  //  Dune::Fem::DomainDecomposedIteratorStorage< CommonTraits::GridPartType > iterators(gridPart);
  //  iterators.update();
  for (auto thread : Dune::XT::Common::value_range(threadnum)) {
    const auto fn = function_name + "_thread_" + Dune::XT::Common::to_string(thread);
    functions[thread] = Dune::XT::Common::make_unique<FVFunc>(fv_space, fn);
  }
  auto& combined = functions.back();
  combined = Dune::XT::Common::make_unique<FVFunc>(fv_space, function_name + "_combined");
  combined->vector() *= 0;

  //  #ifdef _OPENMP
  //  #pragma omp parallel
  //  #endif
  //  {
  //    const auto thread =  Dune::Fem::threadManager().thread();
  //    auto& function = functions[thread];
  //    function->clear();
  //    for(const auto& entity : iterators)
  //    {
  //      auto local_function = function->localFunction(entity);
  //      auto local_combined = combined->localFunction(entity);
  //      for (const auto idx : Dune::XT::Common::value_range(local_function.size())) {
  //        local_function[idx] = thread+1;
  //        local_combined[idx] = thread+1;
  //      }
  //    }
  //  }

  output_all(functions, view_ptr, function_name + "all");
}

template <class T>
unsigned long id_to_ulong(const T& id) {
  return static_cast<unsigned long>(id);
}

template <int k>
unsigned long id_to_ulong(const Dune::bigunsignedint<k>& id) {
  return id.touint();
}

void subgrid_vis(DMP::ProblemContainer& problem, CommonTraits::GridType& coarse_grid,
                 CommonTraits::GridType& fine_grid) {
  CommonTraits::GridProviderType coarse_grid_provider(coarse_grid);
  CommonTraits::GridProviderType fine_grid_provider(fine_grid);
  const CommonTraits::SpaceType coarseSpace(
      CommonTraits::SpaceChooserType::PartViewType::create(coarse_grid, CommonTraits::st_gdt_grid_level));
  const CommonTraits::SpaceType fineSpace(
      CommonTraits::SpaceChooserType::PartViewType::create(fine_grid, CommonTraits::st_gdt_grid_level));

  LocalGridList localgrid_list(problem, coarseSpace);

  auto fine_view_ptr = fine_grid.leafGridView();
  typedef GDT::Spaces::FV::Default<CommonTraits::GridViewType, double, 1, 1> FVSpace;
  typedef BackendChooser<FVSpace>::DiscreteFunctionType FVFunc;
  FVSpace fv_space(fine_view_ptr);

  std::vector<std::unique_ptr<FVFunc>> oversampled_functions(localgrid_list.size());
  std::vector<std::unique_ptr<FVFunc>> functions(localgrid_list.size());

  auto oversampled_function_it = oversampled_functions.begin();
  auto function_it = functions.begin();
  // horrible, horrible complexity :)
  for (const auto& coarse_entity : coarseSpace) {
    const auto& id_set = coarse_grid.globalIdSet();
    const auto coarse_id = id_set.id(coarse_entity);
    auto& oversampled_function = (*oversampled_function_it++);
    oversampled_function = Dune::XT::Common::make_unique<FVFunc>(fv_space, Dune::XT::Common::to_string(coarse_id) + "_subgrid");
    oversampled_function->vector() *= 0;
    auto& function = (*function_it++);
    function = Dune::XT::Common::make_unique<FVFunc>(fv_space, Dune::XT::Common::to_string(coarse_id) + "_coarse_cell");
    function->vector() *= 0;

    for (const auto& fine_entity : fineSpace) {
      if (localgrid_list.covers_strict(coarse_entity, fine_entity.geometry())) {
        auto oversampled_local_function = oversampled_function->local_discrete_function(fine_entity);
        for (const auto idx : Dune::XT::Common::value_range(oversampled_local_function->vector().size())) {
          oversampled_local_function->vector().set(idx, id_to_ulong(coarse_id + 1));
        }
        //        if (coarse_id == localgrid_list.getEnclosingMacroCellId(CommonTraits::EntityPointerType(fine_entity)))
        //        {
        //          auto local_function = function->local_discrete_function(fine_entity);
        //          for (const auto idx : Dune::XT::Common::value_range(local_function.size())) {
        //            local_function.vector().set(idx, coarse_id+1);
        //          }
        //        }
      }
    }
  }

  output_all(oversampled_functions, fine_view_ptr, "subgrids");
  output_all(functions, fine_view_ptr, "coarse_cells");
}

void partition_vis(CommonTraits::GridType& coarse_grid, CommonTraits::GridType& fine_grid) {
  partition_vis_single(coarse_grid, "coarse_grid");
  partition_vis_single(fine_grid, "fine_grid");
} // function algorithm

} // namespace Dune {
} // namespace Multiscale {

int main(int argc, char** argv) {
  using namespace Dune::Multiscale;
  try {
    init(argc, argv);

    assert(Dune::XT::CommonthreadManager().max_threads() == DXTC_CONFIG_GET("threading.max_count", 1u));
    const std::string datadir = DXTC_CONFIG_GET("global.datadir", "data/");

    // generate directories for data output
    Dune::XT::Common::test_create_directory(datadir);

    switch (DXTC_CONFIG_GET("msfem.oversampling_strategy", 1)) {
      case 1:
        break;
      case 2:
        break;
      default:
        DUNE_THROW(Dune::InvalidStateException, "Oversampling Strategy must be 1 or 2.");
    }

    const auto& comm = Dune::MPIHelper::getCommunicator();
    DSC_CONFIG.set("grids.dim", CommonTraits::world_dim);
    DMP::ProblemContainer problem(comm, comm, DSC_CONFIG);

    auto grids = Dune::Multiscale::make_grids(problem, false);
    const auto& coarse_grid = *grids.first;
    problem.getMutableModelData().prepare_new_evaluation(problem);

    if (DXTC_CONFIG_GET("global.vtk_output", false)) {
      problem.getDiffusion().visualize(
          coarse_grid.leafGridView(),
          OutputParameters(problem.config().get("global.datadir", "data")).fullpath("diffusion"));
    }
    //    partition_vis(*grids.first, *grids.second);
    //    subgrid_vis(*grids.first, *grids.second);
    DSG::ElementVisualization::all(*grids.second, datadir + "/fine_element_visualization");
    DSG::ElementVisualization::all(coarse_grid, datadir + "/coarse_element_visualization");
  } catch (Dune::Exception& e) {
    return handle_exception(e);
  } catch (std::exception& s) {
    return handle_exception(s);
  }
  return 0;
} // main
