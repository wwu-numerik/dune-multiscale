// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include <config.h>
#include <dune/multiscale/common/main_init.hh>
#include <dune/multiscale/fem/algorithm.hh>
#include <dune/stuff/common/parallel/helper.hh>
#include <dune/stuff/common/profiler.hh>
#include <thread>

#include <tbb/task_scheduler_init.h>

int main(int argc, char** argv) {
  using namespace Dune::Multiscale;
  try {
    init(argc, argv);

    DSC_PROFILER.startTiming("total_cpu");

    cgfem_algorithm();

    const auto cpu_time =
        DSC_PROFILER.stopTiming("total_cpu", DSC_CONFIG_GET("global.output_walltime", false)) / 1000.f;
    DSC_LOG_INFO_0 << "Total runtime of the program: " << cpu_time << "s" << std::endl;
    DSC_PROFILER.outputTimings("profiler");
  }
  catch (Dune::Exception& e) {
    return handle_exception(e);
  }
  catch (std::exception& s) {
    return handle_exception(s);
  }
  return 0;
} // main
