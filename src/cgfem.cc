// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include <config.h>
#include <dune/multiscale/common/main_init.hh>
#include <dune/multiscale/fem/algorithm.hh>
#include <dune/xt/common/parallel/helper.hh>
#include <dune/xt/common/timings.hh>
#include <dune/xt/common/logging.hh>
#include <dune/xt/common/configuration.hh>
#include <thread>

#include <tbb/task_scheduler_init.h>

int main(int argc, char** argv)
{
  using namespace Dune::Multiscale;
  try {
    init(argc, argv);
    const size_t max_threads = DXTC_CONFIG_GET("threading.max_count", 1);
    tbb::task_scheduler_init sched_init(max_threads);
    Dune::XT::Common::threadManager().set_max_threads(max_threads);

    DXTC_TIMINGS.start("total_cpu");

    cgfem_algorithm();

    const auto cpu_time = DXTC_TIMINGS.stop("total_cpu") / 1000.f;
    MS_LOG_INFO_0 << "Total runtime of the program: " << cpu_time << "s" << std::endl;
    DXTC_TIMINGS.output_per_rank("profiler");
    mem_usage();
    dump_environment();
  } catch (Dune::Exception& e) {
    return handle_exception(e);
  } catch (std::exception& s) {
    return handle_exception(s);
  }
  return 0;
} // main
