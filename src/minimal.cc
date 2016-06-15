// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include <config.h>
#include <dune/multiscale/common/main_init.hh>
#include <dune/multiscale/fem/algorithm.hh>
#include <dune/xt/common/parallel/helper.hh>
#include <dune/xt/common/timings.hh>
#include <dune/xt/common/parallel/threadmanager.hh>
#include <thread>

#include <tbb/blocked_range.h>
#include <vector>
#include <tbb/task_scheduler_init.h>
using namespace std;
using namespace tbb;

class NumberPrinter {
private:
    static void develop(size_t i);

public:

    void operator()(const blocked_range<size_t>& r);
    NumberPrinter(NumberPrinter& x, split);
    NumberPrinter(int highbit);

    void join(const NumberPrinter& y);

};
void NumberPrinter::develop(size_t i)
{
  const auto mm = Dune::XT::Common::threadManager().thread();
  cout << "Thread " << mm << ": " << i << std::endl;

}


void NumberPrinter::operator()( const blocked_range<size_t>& r) {
    for(size_t i=r.begin(); i!=r.end(); ++i)
        develop(i);
}

NumberPrinter::NumberPrinter( NumberPrinter& x, split) {}
NumberPrinter::NumberPrinter(int highbit) {}

void NumberPrinter::join(const NumberPrinter& ) {
}

void printNumbers (int n) {
    NumberPrinter p(n);
    parallel_reduce(blocked_range<size_t>(0,n), p, auto_partitioner());
}

int main(int argc, char** argv) {
  using namespace Dune::Multiscale;
  try {
    init(argc, argv);
    const size_t max_threads = DXTC_CONFIG_GET("threading.max_count", 1);
    tbb::task_scheduler_init sched_init(max_threads);
    Dune::XT::Common::threadManager().set_max_threads(max_threads);

    //!TODO include base in config
    DXTC_TIMINGS.start("msfem.all");

    printNumbers(10000);

//    cgfem_algorithm();

//    const auto cpu_time =
//        DXTC_TIMINGS.stop("total_cpu") / 1000.f;
//    MS_LOG_INFO_0 << "Total runtime of the program: " << cpu_time << "s" << std::endl;
//    DXTC_TIMINGS.outputTimings("profiler");
//    mem_usage();
//    dump_environment();
  } catch (Dune::Exception& e) {
    return handle_exception(e);
  } catch (std::exception& s) {
    return handle_exception(s);
  }
  return 0;
} // main
