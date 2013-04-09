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
#include <boost/filesystem/path.hpp>

/**
 * \brief simple discrete function to disk writer
 * this class isn't type safe in the sense that different appends may append
 * non-convertible discrete function implementations
 * \todo base on discrete's functions write_xdr functionality
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
    : filename_(filename)
    , dir_(DSC_CONFIG_GET("global.datadir", "data"))
    , file_(Dune::Stuff::Common
            ::make_ofstream(dir_ / filename_,
                            std::fstream::trunc | std::fstream::out | std::fstream::binary))
  {
    if(!file_->is_open())
      DUNE_THROW(Dune::IOError, boost::format("cannot open file %s in dir %s for writing") % filename_ % dir_ );
  }

  /**
   * \copydoc DiscreteFunctionReader()
   */
  DiscreteFunctionWriter(const boost::filesystem::path& path)
    : filename_(path.string())
    , dir_(DSC_CONFIG_GET("global.datadir", "data"))
    , file_(Dune::Stuff::Common
            ::make_ofstream(dir_ / filename_,
                            std::fstream::trunc | std::fstream::out | std::fstream::binary))
  {
    if(!file_->is_open())
      DUNE_THROW(Dune::IOError, boost::format("cannot open file %s in dir %s for writing") % filename_ % dir_ );
  }

  ~DiscreteFunctionWriter() {
    if (file_->is_open())
      file_->close();
  }

  bool is_open() const {
    return file_->is_open();
  }

  template < class DiscreteFunctionTraits >
  void append(const Dune::DiscreteFunctionInterface< DiscreteFunctionTraits >& df) {
    assert( file_->is_open() );
    typedef typename Dune::DiscreteFunctionInterface< DiscreteFunctionTraits >::DomainFieldType
      Field;
    for (const Field d : df)
    {
      file_->write( reinterpret_cast< const char* >(&d), sizeof(Field) );
    }
  } // append

  template < class DiscreteFunctionTraits >
  void append(const std::vector<const Dune::DiscreteFunctionInterface< DiscreteFunctionTraits >>& df_vec)  {
    for (const auto& df : df_vec)
      append(df);
  } // append

private:
  const std::string filename_;
  const boost::filesystem::path dir_;
  std::unique_ptr<boost::filesystem::ofstream> file_;
};

/**
 * \brief simple discrete function from disk reader
 * this class isn't type safe in the sense that different appends may append
 * non-convertible discrete function implementations
 * \todo base on discrete's functions write_xdr functionality
 */
class DiscreteFunctionReader
{
  void init() {
    if(file_->is_open())
    {
      // get size of file
      file_->seekg(0, std::fstream::end);
      size_ = file_->tellg();
      file_->seekg(0);
    }
  }

public:
  DiscreteFunctionReader(const std::string filename)
    : filename_(filename)
    , size_(-1)
    , dir_(DSC_CONFIG_GET("global.datadir", "data"))
    , file_(Dune::Stuff::Common
            ::make_ifstream(dir_ / filename_,
                            std::fstream::in | std::fstream::binary))
  {
    init();
  }

  DiscreteFunctionReader(const boost::filesystem::path path)
    : filename_(path.string())
    , size_(-1)
    , dir_(DSC_CONFIG_GET("global.datadir", "data"))
    , file_(Dune::Stuff::Common
            ::make_ifstream(dir_ / filename_,
                            std::fstream::in | std::fstream::binary))
  {
    init();
  }

  ~DiscreteFunctionReader() {
    if ( file_->is_open() )
      file_->close();
  }

  bool is_open() const {
    return file_->is_open();
  }

  long size() const {
    return size_;
  }

  template < class DiscreteFunctionTraits >
  void read(const unsigned long index,
            Dune::DiscreteFunctionInterface< DiscreteFunctionTraits >& df) {
    if(!is_open())
      DUNE_THROW(Dune::IOError, boost::format("cannot open file %s in dir %s for reading") % filename_ % dir_ );
    typedef typename Dune::DiscreteFunctionInterface< DiscreteFunctionTraits >::DomainFieldType
      Field;
    const unsigned long bytes = df.size() * sizeof(Field);

    assert( file_->is_open() );
    assert( size_ >= long( bytes * (index + 1) ) );
    file_->seekg(bytes * index);

    for (auto& dof : df)
    {
      file_->read( reinterpret_cast< char* >( &(dof) ), sizeof(Field) );
    }
  } // read

  /*template < class DFType >
     * void read( std::vector<DFType>& df_vec )
     * {
     *  for ( Iter it = df_vec.begin(); it != end; ++it ) {
     * append( *it );
     *  }
     * }*/

private:
  const std::string filename_;
  long size_;
  const boost::filesystem::path dir_;
  std::unique_ptr<boost::filesystem::ifstream> file_;
};

#endif // ifndef DISCRETEFUNCTIONWRITER_HEADERGUARD
