#include <config.h>
#include <config.h>
// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

// The following FEM code requires an access to the 'ModelProblemData' class,
// which provides us with information about f, A, \Omega, etc.

#include <dune/multiscale/common/main_init.hh>
#include <dune/multiscale/common/grid_creation.hh>
#include <dune/multiscale/fem/algorithm.hh>
#include <dune/multiscale/problems/selector.hh>
#include <dune/multiscale/fem/fem_traits.hh>

int main(int argc, char** argv) {
  try {
    using namespace Dune::Multiscale;
    using namespace Dune::Multiscale::FEM;
    init(argc, argv);

    DSC_PROFILER.startTiming("total_cpu");

    const std::string path = DSC_CONFIG_GET("global.datadir", "data/");

    // generate directories for data output
    DSC::testCreateDirectory(path);

    // name of the error file in which the data will be saved
    std::string filename_;
    const auto save_filename = std::string(path + "/logdata/ms.log.log");
    DSC_LOG_INFO << "LOG FILE\n" << "Data will be saved under: " << save_filename << std::endl;

    const auto grids = Dune::Multiscale::make_grids();
    algorithm(grids.second, filename_);

    const auto cpu_time = DSC_PROFILER.stopTiming("total_cpu",
                                                  DSC_CONFIG_GET("global.output_walltime", false)) / 1000.f;
    DSC_LOG_INFO << "Total runtime of the program: " << cpu_time << "s" << std::endl;
    return 0;
  }
  catch (Dune::Exception& e) {
    std::cerr << e.what() << std::endl;
  }
  return 1;
} // main
