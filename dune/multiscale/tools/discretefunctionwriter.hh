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
#include <dune/multiscale/msfem/msfem_traits.hh>
#include <dune/common/deprecated.hh>
#include <dune/common/exceptions.hh>
#include <dune/stuff/common/configuration.hh>
#include <dune/stuff/common/filesystem.hh>
#include <dune/stuff/common/ranges.hh>
#include <dune/stuff/aliases.hh>
#include <dune/stuff/common/memory.hh>
#include <dune/stuff/common/type_utils.hh>

#include <boost/filesystem/path.hpp>

namespace Dune {
namespace Multiscale {

template <class DiscreteFunctionType>
class DiscreteFunctionIO : public boost::noncopyable {
//  static_assert(std::is_base_of<Dune::Fem::IsDiscreteFunction, DiscreteFunctionType>::value, "");

  typedef DiscreteFunctionIO<DiscreteFunctionType> ThisType;
  typedef std::shared_ptr<DiscreteFunctionType> DiscreteFunction_ptr;
  typedef typename DiscreteFunctionType::SpaceType DiscreteFunctionSpaceType;
  typedef std::vector<DiscreteFunction_ptr> Vector;
  typedef typename DiscreteFunctionSpaceType::GridViewType GridViewType;
  typedef std::shared_ptr<const GridViewType> GridViewPtrType;

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

      Fem::XDRFileOutStream ss(fn);
      df->write(ss);

    }

    void read(const unsigned long index, DiscreteFunction_ptr& df) {
      const std::string fn = (dir_ / DSC::toString(index)).string();

      Fem::XDRFileInStream ss(fn);
      df->read(ss);

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
    MemoryBackend(const GridViewPtrType& grid_view, const std::string filename = "nonsense_default_for_map")
      : dir_(boost::filesystem::path(DSC_CONFIG_GET("global.datadir", "data")) / filename)
      , local_grid_provider_(grid_view->grid())
      , space_(DMM::MsFEMTraits::SpaceProviderType::create(local_grid_provider_, CommonTraits::st_gdt_grid_level))
    {}

    void append(const DiscreteFunction_ptr& df) { functions_.push_back(df); }

    void read(const unsigned long index, DiscreteFunction_ptr& df) {
      if (index < functions_.size()) {
        df = functions_.at(index);
      } else
        DUNE_THROW(InvalidStateException, "requesting function at oob index");
      assert(df != nullptr);
    }

    DiscreteFunctionSpaceType& space() { return space_; }

  private:
    const boost::filesystem::path dir_;
    const DMM::MsFEMTraits::LocalGridProviderType local_grid_provider_;
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
#if HAVE_EMPLACE
    auto ret = map.emplace(filename, std::move(ptr));
#else 
    auto ret = map.insert(std::make_pair(filename, std::move(ptr)));
#endif
    assert(ret.second);
    return ret.first->second;
  }

  DiskBackend& get_disk(const std::string filename) { return *get(disk_, filename, filename); }

  MemoryBackend& get_memory(const std::string filename,  const GridViewPtrType& grid_view) {
    return *get(memory_, filename, grid_view, filename);
  }

public:
  static MemoryBackend& memory(const std::string filename, const GridViewPtrType& grid_view) {
    return instance().get_memory(filename, grid_view);
  }

  static DiskBackend& disk(const std::string filename) { return instance().get_disk(filename); }

  //! this needs to be called before global de-init or else dune fem fails
  static void clear() {
    auto& th = instance();
    DSC_LOG_DEBUG << (boost::format("cleared %d in-memory functions\ncleared %d on-disk   functions\nfor %s\n")
                      % th.memory_.size() % th.disk_.size() % DSC::getTypename(th)).str();
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
