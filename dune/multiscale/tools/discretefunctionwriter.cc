#include <config.h>
#include "discretefunctionwriter.hh"

Dune::Multiscale::DiskBackend&
Dune::Multiscale::DiscreteFunctionIO::get_disk(const Stuff::Common::Configuration& config, std::string filename) {
  return *get(disk_, filename, config, filename);
}

Dune::Multiscale::MemoryBackend&
Dune::Multiscale::DiscreteFunctionIO::get_memory(std::string filename,
                                                 Dune::Multiscale::IOTraits::GridViewType& grid_view) {
  return *get(memory_, filename, grid_view, filename);
}

Dune::Multiscale::MemoryBackend&
Dune::Multiscale::DiscreteFunctionIO::memory(std::string filename,
                                             Dune::Multiscale::IOTraits::GridViewType& grid_view) {
  return instance().get_memory(filename, grid_view);
}

Dune::Multiscale::DiskBackend& Dune::Multiscale::DiscreteFunctionIO::disk(const DSC::Configuration& config,
                                                                          std::string filename) {
  return instance().get_disk(config, filename);
}

void Dune::Multiscale::DiscreteFunctionIO::clear() {
  auto& th = instance();
  DSC_LOG_DEBUG << (boost::format("cleared %d in-memory functions\ncleared %d "
                                  "on-disk   functions\nfor %s\n") %
                    th.memory_.size() % th.disk_.size() % DSC::getTypename(th)).str();
  th.memory_.clear();
  th.disk_.clear();
}
