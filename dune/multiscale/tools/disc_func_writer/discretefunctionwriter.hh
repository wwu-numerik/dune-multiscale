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

#include <dune/common/deprecated.hh>
#include <dune/common/exceptions.hh>
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
    , file_(Dune::Stuff::Common::Filesystem
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
    , file_(Dune::Stuff::Common::Filesystem
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

  template< class DFType >
  void append(const DFType& df) {
    assert( file_->is_open() );
    typedef typename DFType::DomainFieldType
    Field;
    typedef typename DFType::ConstDofIteratorType
    Iter;
    Iter end = df.dend();
    for (Iter it = df.dbegin(); it != end; ++it)
    {
      double d = *it;
      file_->write( reinterpret_cast< char* >(&d), sizeof(Field) );
    }
  } // append

  template< class DFType >
  void append(const std::vector< DFType >& df_vec) {
    typedef typename std::vector< DFType >::const_iterator
    Iter;
    Iter end = df_vec.end();
    for (Iter it = df_vec.begin(); it != end; ++it)
    {
      append(*it);
    }
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
    else
      DUNE_THROW(Dune::IOError, boost::format("cannot open file %s in dir %s for reading") % filename_ % dir_ );
  }

public:
  DiscreteFunctionReader(const std::string filename)
    : filename_(filename)
    , size_(-1)
    , dir_(DSC_CONFIG_GET("global.datadir", "data"))
    , file_(Dune::Stuff::Common::Filesystem
            ::make_ifstream(dir_ / filename_,
                            std::fstream::in | std::fstream::binary))
  {
    init();
  }

  DiscreteFunctionReader(const boost::filesystem::path path)
    : filename_(path.string())
    , size_(-1)
    , dir_(DSC_CONFIG_GET("global.datadir", "data"))
    , file_(Dune::Stuff::Common::Filesystem
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

  template< class DFType >
  void read(const unsigned long index, DFType& df) {
    typedef typename DFType::DomainFieldType
    Field;
    const unsigned long bytes = df.size() * sizeof(Field);

    assert( file_->is_open() );
    assert( size_ >= long( bytes * (index + 1) ) );
    file_->seekg(bytes * index);

    typedef typename DFType::DofIteratorType
    Iter;
    Iter end = df.dend();
    for (Iter it = df.dbegin(); it != end; ++it)
    {
      file_->read( reinterpret_cast< char* >( &(*it) ), sizeof(Field) );
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
