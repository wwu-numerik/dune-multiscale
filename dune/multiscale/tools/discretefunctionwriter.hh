// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

/**
   *  \file discretefunctionwriter.hh
   *  \brief  write a bunch of discrete functions to one file and retrieve 'em
   **/

#ifndef DISCRETEFUNCTIONWRITER_HEADERGUARD
#define DISCRETEFUNCTIONWRITER_HEADERGUARD

#ifdef HAVE_CMAKE_CONFIG
 #include "cmake_config.h"
#elif defined (HAVE_CONFIG_H)
 #include <config.h>
#endif // ifdef HAVE_CMAKE_CONFIG

#include <fstream>
#include <vector>
#include <cassert>
#include <memory>

#include <dune/common/deprecated.hh>
#include <dune/common/exceptions.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/stuff/common/parameter/configcontainer.hh>
#include <dune/stuff/common/filesystem.hh>
#include <dune/stuff/common/ranges.hh>
#include <dune/stuff/aliases.hh>

#include <boost/filesystem/path.hpp>

/**
 * \brief simple discrete function to disk writer
 * this class isn't type safe in the sense that different appends may append
 * non-convertible discrete function implementations
 */
class DiscreteFunctionWriter
{
public:
  /**
   * \brief DiscreteFunctionWriter
   * \param filename will open fstream at config["global.datadir"]/filename
   *  filename may include additional path components
   * \throws Dune::IOError if config["global.datadir"]/filename cannot be opened
   */
  DiscreteFunctionWriter(const std::string filename)
    : dir_(boost::filesystem::path(DSC_CONFIG_GET("global.datadir", "data")) / filename)
    , size_(0)
  {
    DSC::testCreateDirectory(dir_.string());
  }

  /**
   * \copydoc DiscreteFunctionReader()
   */
  DiscreteFunctionWriter(const boost::filesystem::path& path)
    : dir_(boost::filesystem::path(DSC_CONFIG_GET("global.datadir", "data")) / path)
    , size_(0)
  {
    DSC::testCreateDirectory(dir_.string());
  }

  template < class DiscreteFunctionTraits >
  void append(const Dune::DiscreteFunctionInterface< DiscreteFunctionTraits >& df) {
    const std::string fn = (dir_ / DSC::toString(size_++)).string();
    DSC::testCreateDirectory(fn);
    Dune::Fem::BinaryFileOutStream stream(fn);
    df.write(stream);
  } // append

  template < class DiscreteFunctionTraits >
  void append(const std::vector<const Dune::DiscreteFunctionInterface< DiscreteFunctionTraits >>& df_vec)  {
    for (const auto& df : df_vec)
      append(df);
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
class DiscreteFunctionReader
{

public:
  DiscreteFunctionReader(const std::string filename)
    : size_(0)
    , dir_(boost::filesystem::path(DSC_CONFIG_GET("global.datadir", "data")) / filename)
  {}

  DiscreteFunctionReader(const boost::filesystem::path path)
    : size_(0)
    , dir_(boost::filesystem::path(DSC_CONFIG_GET("global.datadir", "data")) / path)
  {}

  long size() const {
    return size_;
  }

  template < class DiscreteFunctionTraits >
  void read(const unsigned long index,
            Dune::DiscreteFunctionInterface< DiscreteFunctionTraits >& df) {
    const std::string fn = (dir_ / DSC::toString(index)).string();
    Dune::Fem::BinaryFileInStream stream(fn);
    df.read(stream);
  } // read

private:
  long size_;
  const boost::filesystem::path dir_;
};

#endif // ifndef DISCRETEFUNCTIONWRITER_HEADERGUARD
