#include <config.h>
#include "df_io.hh"
#include <dune/xt/common/string.hh>

Dune::Multiscale::DiskBackend&
Dune::Multiscale::DiscreteFunctionIO::get_disk(const XT::Common::Configuration& config, std::string filename) {
  return *get(disk_, filename, config, filename);
}

Dune::Multiscale::MemoryBackend&
Dune::Multiscale::DiscreteFunctionIO::get_memory(std::string filename,
                                                 Dune::Multiscale::IOTraits::GridViewType& grid_view) {
  const auto tokens = Dune::XT::Common::tokenize(filename, "_");
  const size_t idx = Dune::XT::Common::from_string<size_t>(tokens.back());
  return *get(memory_, idx, grid_view, filename);
}

Dune::Multiscale::MemoryBackend&
Dune::Multiscale::DiscreteFunctionIO::memory(std::string filename,
                                             Dune::Multiscale::IOTraits::GridViewType& grid_view) {
  return instance().get_memory(filename, grid_view);
}

Dune::Multiscale::DiskBackend& Dune::Multiscale::DiscreteFunctionIO::disk(const Dune::XT::Common::Configuration& config,
                                                                          std::string filename) {
  return instance().get_disk(config, filename);
}

void Dune::Multiscale::DiscreteFunctionIO::clear() {
  auto& th = instance();
  MS_LOG_DEBUG << (boost::format("cleared %d in-memory functions\ncleared %d "
                                 "on-disk   functions\nfor %s\n") %
                   th.memory_.size() % th.disk_.size() % Dune::XT::Common::get_typename(th))
                      .str();
  th.memory_.clear();
  th.disk_.clear();
}
