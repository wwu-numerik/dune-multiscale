// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

/**
   *  \file discretefunctionwriter.hh
   *  \brief  write a bunch of discrete functions to one file and retrieve 'em
   **/

#ifndef DISCRETEFUNCTIONWRITER_HEADERGUARD
#define DISCRETEFUNCTIONWRITER_HEADERGUARD

#include <fstream>
#include <vector>
#include <cassert>
#include <memory>
#include <unordered_map>

#include <dune/multiscale/common/traits.hh>
#include <dune/common/deprecated.hh>
#include <dune/common/exceptions.hh>
#include <dune/stuff/common/parameter/configcontainer.hh>
#include <dune/stuff/common/filesystem.hh>
#include <dune/stuff/common/ranges.hh>
#include <dune/stuff/aliases.hh>
#include <dune/stuff/common/memory.hh>
#include <dune/fem/io/streams/xdrstreams.hh>

#include <boost/filesystem/path.hpp>

namespace Dune {
namespace Multiscale {

//! tiny struct to ensure i/o type don't diverge
struct IOTraits {
#if HAVE_SIONLIB&& HAVE_MPI
#define MULTISCALE_USE_SION
  typedef Dune::Fem::SIONlibOutStream OutstreamType;
  typedef Dune::Fem::SIONlibInStream InstreamType;
#endif
};

template <class DiscreteFunctionType>
class DiscreteFunctionIO : public boost::noncopyable {
  static_assert(std::is_base_of<Dune::Fem::IsDiscreteFunction, DiscreteFunctionType>::value, "");

  typedef DiscreteFunctionIO<DiscreteFunctionType> ThisType;
  typedef std::shared_ptr<DiscreteFunctionType> DiscreteFunction_ptr;
  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
  typedef std::vector<DiscreteFunction_ptr> Vector;
  typedef typename DiscreteFunctionType::GridPartType GridPartType;

  DiscreteFunctionIO() = default;

  class DiskBackend : public boost::noncopyable {

    void load_disk_functions() {
      DSC::testCreateDirectory(dir_.string());
      // if functions present, load em
    }

  public:
    /**
     * \brief DiscreteFunctionWriter
     * \param filename will open fstream at config["global.datadir"]/filename
     *  filename may include additional path components
     * \throws Dune::IOError if config["global.datadir"]/filename cannot be opened
     */
    DiskBackend(const std::string filename = "nonsense_default_for_map")
      : dir_(boost::filesystem::path(DSC_CONFIG_GET("global.datadir", "data")) / filename)
      , index_(0) {}

    void append(const DiscreteFunction_ptr& df) {
      const std::string fn = (dir_ / DSC::toString(index_++)).string();
      DSC::testCreateDirectory(fn);
#ifdef MULTISCALE_USE_SION
      IOTraits::OutstreamType stream(fn);
      df->write(stream);
#else
      Fem::XDRFileOutStream ss(fn);
      df->write(ss);
#endif
    }

    void read(const unsigned long index, DiscreteFunction_ptr& df) {
      const std::string fn = (dir_ / DSC::toString(index)).string();
#ifdef MULTISCALE_USE_SION
      IOTraits::InstreamType stream(fn);
      df->read(stream);
#else
      Fem::XDRFileInStream ss(fn);
      df->read(ss);
#endif
    }

  private:
    const boost::filesystem::path dir_;
    unsigned int index_;
  };

  /**
   * \brief simple discrete function to disk writer
   * this class isn't type safe in the sense that different appends may append
   * non-convertible discrete function implementations
   */
  class MemoryBackend : public boost::noncopyable {

  public:
    /**
     * \brief DiscreteFunctionWriter
     * \param filename will open fstream at config["global.datadir"]/filename
     *  filename may include additional path components
     * \throws Dune::IOError if config["global.datadir"]/filename cannot be opened
     */
    MemoryBackend(typename GridPartType::GridType& grid, const std::string filename = "nonsense_default_for_map")
      : dir_(boost::filesystem::path(DSC_CONFIG_GET("global.datadir", "data")) / filename)
      , grid_part_(grid)
      , space_(grid_part_) {}

    void append(const DiscreteFunction_ptr& df) { functions_.push_back(df); }

    void read(const unsigned long index, DiscreteFunction_ptr& df) {
      if (index < functions_.size()) {
        df = functions_.at(index);
      } else
        DUNE_THROW(InvalidStateException, "requesting function at oob index");
      assert(df != nullptr);
    }

    GridPartType& grid_part() { return grid_part_; }
    DiscreteFunctionSpaceType& space() { return space_; }

  private:
    const boost::filesystem::path dir_;
    GridPartType grid_part_;
    DiscreteFunctionSpaceType space_;
    Vector functions_;
  };

  static ThisType& instance() {
    static ThisType s_this;
    return s_this;
  }

  template <class IOMapType, class... Args>
  typename IOMapType::mapped_type& get(IOMapType& map, std::string filename, Args&&... ctor_args) {
    auto it = map.find(filename);
    if (it != map.end())
      return it->second;
    std::lock_guard<std::mutex> lock(mutex_);
    auto ptr = std::make_shared<typename IOMapType::mapped_type::element_type>(ctor_args...);
    auto ret = map.emplace(filename, std::move(ptr));
    assert(ret.second);
    return ret.first->second;
  }

  DiskBackend& get_disk(const std::string filename) { return *get(disk_, filename, filename); }

  MemoryBackend& get_memory(const std::string filename, typename GridPartType::GridType& grid) {
    return *get(memory_, filename, grid, filename);
  }

public:
  static MemoryBackend& memory(const std::string filename, typename GridPartType::GridType& grid) {
    return instance().get_memory(filename, grid);
  }

  static DiskBackend& disk(const std::string filename) { return instance().get_disk(filename); }

  //! this needs to be called before global de-init or else dune fem fails
  static void clear() {
    auto& th = instance();
    th.memory_.clear();
    th.disk_.clear();
  }

private:
  std::unordered_map<std::string, std::shared_ptr<MemoryBackend>> memory_;
  std::unordered_map<std::string, std::shared_ptr<DiskBackend>> disk_;
  std::mutex mutex_;

}; // class DiscreteFunctionIO

} // namespace Multiscale {
} // namespace Dune {

#endif // ifndef DISCRETEFUNCTIONWRITER_HEADERGUARD
