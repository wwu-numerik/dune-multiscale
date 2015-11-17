#include <config.h>
#include "discretefunctionwriter.hh"
#include <dune/stuff/common/string.hh>

Dune::Multiscale::DiskBackend&
Dune::Multiscale::DiscreteFunctionIO::get_disk(const Stuff::Common::Configuration& config, std::string filename) {
  return *get(disk_, filename, config, filename);
}

Dune::Multiscale::MemoryBackend&
Dune::Multiscale::DiscreteFunctionIO::get_memory(std::string filename,
                                                 Dune::Multiscale::IOTraits::GridViewType& grid_view) {
  const auto tokens = DSC::tokenize(filename, "_");
  const size_t idx = DSC::fromString<size_t>(tokens.back());
  return *get(memory_, idx, grid_view, filename);
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
  MS_LOG_DEBUG << (boost::format("cleared %d in-memory functions\ncleared %d "
                                  "on-disk   functions\nfor %s\n") %
                    th.memory_.size() % th.disk_.size() % DSC::getTypename(th)).str();
  th.memory_.clear();
  th.disk_.clear();
}
