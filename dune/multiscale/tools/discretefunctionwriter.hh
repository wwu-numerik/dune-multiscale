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

#include <dune/multiscale/common/traits.hh>
#include <dune/common/deprecated.hh>
#include <dune/common/exceptions.hh>
#include <dune/stuff/common/parameter/configcontainer.hh>
#include <dune/stuff/common/filesystem.hh>
#include <dune/stuff/common/ranges.hh>
#include <dune/stuff/aliases.hh>

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


class DiscreteFunctionIO {
  /**
   * \brief simple discrete function to disk writer
   * this class isn't type safe in the sense that different appends may append
   * non-convertible discrete function implementations
   */
  class DiscreteFunctionWriter {
  public:
    /**
     * \brief DiscreteFunctionWriter
     * \param filename will open fstream at config["global.datadir"]/filename
     *  filename may include additional path components
     * \throws Dune::IOError if config["global.datadir"]/filename cannot be opened
     */
    DiscreteFunctionWriter(const std::string filename)
      : dir_(boost::filesystem::path(DSC_CONFIG_GET("global.datadir", "data")) / filename)
      , size_(0) {
      DSC::testCreateDirectory(dir_.string());
    }

    /**
     * \copydoc DiscreteFunctionReader()
     */
    DiscreteFunctionWriter(const boost::filesystem::path& path)
      : dir_(boost::filesystem::path(DSC_CONFIG_GET("global.datadir", "data")) / path)
      , size_(0) {
      DSC::testCreateDirectory(dir_.string());
    }

    template <class DiscreteFunctionTraits>
    void append(const Dune::Fem::DiscreteFunctionInterface<DiscreteFunctionTraits>& df) {
      const std::string fn = (dir_ / DSC::toString(size_++)).string();
      DSC::testCreateDirectory(fn);
  #ifdef MULTISCALE_USE_SION
      IOTraits::OutstreamType stream(fn);
      df.write(stream);
  #else
      df.write_xdr(fn);
  #endif
    } // append

  private:
    const boost::filesystem::path dir_;
    unsigned int size_;
  };

  /**
   * \brief simple discrete function from disk reader
   * this class isn't type safe in the sense that different appends may append
   * non-convertible discrete function implementations
   * \todo base on discrete's functions write_xdr functionality
   */
  class DiscreteFunctionReader {

  public:
    DiscreteFunctionReader(const std::string filename)
      : dir_(boost::filesystem::path(DSC_CONFIG_GET("global.datadir", "data")) / filename) {}

    DiscreteFunctionReader(const boost::filesystem::path path)
      : dir_(boost::filesystem::path(DSC_CONFIG_GET("global.datadir", "data")) / path) {}



    template <class DiscreteFunctionTraits>
    void read(const unsigned long index, Dune::Fem::DiscreteFunctionInterface<DiscreteFunctionTraits>& df) {
      const std::string fn = (dir_ / DSC::toString(index)).string();
  #ifdef MULTISCALE_USE_SION
      IOTraits::InstreamType stream(fn);
      df.read(stream);
  #else
      df.read_xdr(fn);
  #endif
    } // read

  private:
    const boost::filesystem::path dir_;
  };


public:
  static DiscreteFunctionReader reader(const std::string filename) {
    return DiscreteFunctionReader(filename);
  }

  static DiscreteFunctionWriter writer(const std::string filename) {
    return DiscreteFunctionWriter(filename);
  }

};//class DiscreteFunctionIO

} // namespace Multiscale {
} // namespace Dune {

#endif // ifndef DISCRETEFUNCTIONWRITER_HEADERGUARD
