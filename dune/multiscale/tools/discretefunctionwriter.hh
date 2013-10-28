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

#include <boost/filesystem/path.hpp>

namespace Dune {
namespace Multiscale {

//! tiny struct to ensure i/o type don't diverge
struct IOTraits {
#if HAVE_SIONLIB && HAVE_MPI
#define MULTISCALE_USE_SION
  typedef Dune::Fem::SIONlibOutStream OutstreamType;
  typedef Dune::Fem::SIONlibInStream InstreamType;
#endif
};

template <class DiscreteFunctionType>
class DiscreteFunctionIO {
  static_assert(std::is_base_of<Dune::Fem::IsDiscreteFunction, DiscreteFunctionType>::value, "");

  typedef DiscreteFunctionIO<DiscreteFunctionType> ThisType;
  typedef std::shared_ptr<DiscreteFunctionType> DiscreteFunction_ptr;
  typedef std::vector<DiscreteFunction_ptr> Vector;

  /**
   * \brief simple discrete function to disk writer
   * this class isn't type safe in the sense that different appends may append
   * non-convertible discrete function implementations
   */
  class DiscreteFunctionRW{

    void load_disk_functions()
    {
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
    DiscreteFunctionRW(const std::string filename = "nonsense_default_for_map")
      : dir_(boost::filesystem::path(DSC_CONFIG_GET("global.datadir", "data")) / filename) {
      load_disk_functions();
    }

    /**
     * \copydoc DiscreteFunctionReader()
     */
    DiscreteFunctionRW(const boost::filesystem::path& path)
      : dir_(boost::filesystem::path(DSC_CONFIG_GET("global.datadir", "data")) / path) {
      load_disk_functions();
    }

    ~DiscreteFunctionRW()
    {
      //dump remaining
      unsigned long id = 0;
      for(auto& df : functions_)
      {
        to_disk(id++, df);
      }
    }

    void append(const DiscreteFunction_ptr& df) {
      functions_.push_back(df);
    }

    void read(const unsigned long index, DiscreteFunction_ptr& df) {
      if(index<functions_.size()) {
        df = functions_.at(index);
        functions_.erase(functions_.begin()+index);
      } else {
        from_disk(index, df);
      }
      assert(df!=nullptr);
    }

    void to_disk(const unsigned long index, const DiscreteFunction_ptr& df) const
    {
      const std::string fn = (dir_ / DSC::toString(index)).string();
      DSC::testCreateDirectory(fn);
  #ifdef MULTISCALE_USE_SION
      IOTraits::OutstreamType stream(fn);
      df->write(stream);
  #else
      df->write_xdr(fn);
  #endif
    }

    void from_disk(const unsigned long index, const DiscreteFunction_ptr& df) const {
      const std::string fn = (dir_ / DSC::toString(index)).string();
  #ifdef MULTISCALE_USE_SION
      IOTraits::InstreamType stream(fn);
      df->read(stream);
  #else
      df->read_xdr(fn);
  #endif
    } // read

  private:
    const boost::filesystem::path dir_;
    Vector functions_;
  };

  static ThisType& instance() {
    static ThisType s_this;
    return s_this;
  }

  template <class IOMapType>
  typename IOMapType::mapped_type& get(IOMapType& map, std::string filename)
  {
    auto it = map.find(filename);
    if(it != map.end())
      return it->second;
    auto ret = map.emplace(filename,filename);
    assert(ret.second);
    return ret.first->second;
  }

  DiscreteFunctionRW& get_rw(const std::string filename) {
    return get(rws_, filename);
  }
public:
  static DiscreteFunctionRW& instance(const std::string filename) {
    return instance().get_rw(filename);
  }

private:
  std::unordered_map<std::string, DiscreteFunctionRW> rws_;

};//class DiscreteFunctionIO

} // namespace Multiscale {
} // namespace Dune {

#endif // ifndef DISCRETEFUNCTIONWRITER_HEADERGUARD
