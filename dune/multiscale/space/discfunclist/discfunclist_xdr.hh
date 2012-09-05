#ifndef DISCFUNCLIST_XDR_HH
#define DISCFUNCLIST_XDR_HH

#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>
#include <string>

#include "discfunclistinterface.hh"
#include <dune/fem/io/file/iointerface.hh>

namespace Dune {
namespace Multiscale {
using std::vector;
using std::pair;
using std::make_pair;

/** \class   DiscreteFunctionList_xdr
   *  \ingroup DiscFuncList
   *  \brief   Stores functions together with attributes in a list
   *
   *  The function is stored to disc using the functions memberfunction "write_xdr", that is:
   *  it is written with xdr encoding. You can also store an attribute together with each function
   *  (see interface documentation for details). This attribute is stored on disc aswell using the AttributeActions
   *  class. If you would like to use an AttributeType other than int, double, float, string or
   *  char* you will have to implement the AttributeActions class for that AttributeType.
   *  For convinience a header file is written which stores the following information about all
   *  functions handled by the current instance of this class:
   *
   *       index;Attributefilename;Datafilename
   *
   *  So one might change this file to change the list (remove, add, resort functions etc.)
   *
   *  \param DiscreteFunctionListTraits a traits class that defines the discrete function
   *                                    type
   */
template< class DiscreteFunctionListTraits >
class DiscreteFunctionList_xdr
  : public
    DiscreteFunctionListInterface< DiscreteFunctionListTraits >
{
  // typedefs
  typedef DiscreteFunctionListTraits                             Traits;
  typedef DiscreteFunctionList_xdr< DiscreteFunctionListTraits > ThisType;
  typedef DiscreteFunctionListInterface
  < DiscreteFunctionListTraits > BaseType;

public:
  typedef typename Traits::DiscreteFunctionType DiscreteFunctionType;
  typedef typename DiscreteFunctionType
    ::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
  /*the next line determines whether AttributeType is defined in Traits class, if it is not defined,
     * int is used as default.*/
  typedef typename SelectTypeAttributeTypeIfEnabled
  < DiscreteFunctionListTraits, int >::Type AttributeType;

private:
  void readOrCreateHeaderFile(const DiscreteFunctionSpaceType& discFuncSpace,
                              const std::string name, const std::string header,
                              const bool create) {
    IOInterface::createGlobalPath(discFuncSpace_.gridPart().grid().comm(), dataPath_);

    // save header file name for later use
    headerfileName_ = dataPath_ + header;
    // open inward file stream to read indices and attributes to this list...
    const char* headerfile = headerfileName_.c_str();
    std::ifstream in(headerfile);
    if (!in || create)
    {
      header_.open( headerfileName_.c_str() );
      if (!header_)
      {
        // ...and throw an exeption if the headerfile doesn't exist
        DUNE_THROW(IOError, "Could not open " << headerfile << " for reading!");
      }
      header_ << "# This is the header file for the discrete function list \n"
              << "# You may change the list by changing this file. Directory must contain\n"
              << "# a '.dat' and a '.att' file named as the entry in this list. You may also\n"
              << "# use the '[a:b:c]' syntax, meaning the list will contain each function from a.dat\n"
              << "# to c.dat, counting in b stepsize\n";
      header_.flush();

      return;
    }

    std::string row;
    std::string attributefile, datafile;
    AttributeType attr;
    std::string filename;
    unsigned int index;
    // read given header file
    while ( getline(in, row) )
    {
      // ignore blank lines
      if (strlen( row.c_str() ) == 0)
        continue;

      // ignore comments
      if (row[0] == '#')
        continue;

      // wrap matlab-like "[a:b:c]" terms
      if (row.find('[') != std::string::npos)
      {
        // we found an occurrence of a "[a:b:c]" term, delete all preceding text
        row.erase(0, row.find('[') + 1);
        int a, b, c;
        // \todo Only numbers of length one are covered here!!!
        a = atoi(&row[0]);
        b = atoi(&row[2]);
        c = atoi(&row[4]);

        for (int i = a; i <= c; i += b)
        {
          index = nextIndex();
          std::stringstream foo;
          foo << a;
          filename = foo.str();

          // read datafile name ...

          datafile = dataPath_;
          datafile += filename;
          datafile += ".dat";
          // and test if datafile exists in current path...
          std::ifstream test( datafile.c_str() );
          if (!test)
          {
            // ...throw an exeption if it can't be opened
            DUNE_THROW(IOError, "Could not open " << datafile << " for reading!");
          }
          test.close();

          // save filename
          filenames_.push_back( make_pair(index, filename) );

          // get attributefile name
          attributefile = dataPath_;
          attributefile += filename;
          attributefile += ".att";
          // and read attribute from attribute file, existence is checked in AttributeActions::read
          attrAct.read(attr, attributefile);

          attributes_.push_back( make_pair(index, attr) );
        }
        continue;
      }

      index = nextIndex();
      std::string::size_type pos = row.find('.');
      // if there is a line containing "." read only the name without extension
      if (pos != std::string::npos)
        filename = row.substr(0, pos);
      else
        filename = row;

      // construct datafile name ...
      datafile = dataPath_;
      datafile += filename;
      datafile += ".dat";
      // and test if datafile exists in current path...
      std::ifstream test( datafile.c_str() );
      if (!test)
      {
        // ...throw an exeption if it can't be opened
        DUNE_THROW(IOError, "Could not open " << datafile << " for reading!");
      }
      test.close();

      // save filename
      filenames_.push_back( make_pair(index, filename) );
      // construct attributefile name
      attributefile = dataPath_;
      attributefile += filename;
      attributefile += ".att";
      // and read attribute from attribute file, existence is checked in AttributeActions::read
      attrAct.read(attr, attributefile);

      attributes_.push_back( make_pair(index, attr) );
    }

    // if header file was empty, give a warning
    if (attributes_.size() == 0)
      DSC_LOG_ERROR << "Warning! Headerfile " << headerfile << " is empty!\n";

    // open header for this list
    header_.open(headerfile, std::ios::out | std::ios::app);

    DSC_LOG_INFO << "read discfunclist_xdr from headerfile, size = " << size() << std::endl;
  } // readOrCreateHeaderFile

public:
  /** \brief constructor
     *
     *  checks for existence of data path and tries to create it if it
     *  doesn't exist
     *  \param[in] discFuncSpace discrete function space where the functions are defined
     *  \param[in] name name for the list, this is also the pathname and the name for the header file
     *  \param[in] create  boolean flag indicating wether the discrete function list shall be created newly.
     */
  DiscreteFunctionList_xdr(const DiscreteFunctionSpaceType& discFuncSpace,
                           std::string name = "discFuncs",
                           const bool create = false)
    : dataPath_("./" + name + "/")
      , discFuncSpace_(discFuncSpace)
      , cleared_(false) {
    std::string headerfile(name + ".lst");

    readOrCreateHeaderFile(discFuncSpace, name, headerfile, create);
  }

private:
  /** \brief constructor
     *
     * checks for existence of data path and tries to create it if
     * it doesn't exist.
     * In addition, this constructor takes a number as argument. With this you can tell the list
     * the number of functions that will be stored. The list will still be able to store more functions,
     * but continuous reallocation of memory can be avoided.
     *  \param[in] discFuncSpace discrete function space where the functions are defined
     *  \param[in] length number of functions that will be stored approximatly in this list
     *  \param[in] name name for the list, this is also the pathname and the name for the header file
     *  \param[in] create  boolean flag indicating wether the discrete function list shall be created newly.
     */
  DiscreteFunctionList_xdr(const DiscreteFunctionSpaceType& discFuncSpace,
                           int length,
                           std::string name = "discFuncs",
                           const bool create = false)
    : dataPath_("./" + name + "/")
      , discFuncSpace_(discFuncSpace)
      , cleared_(false) {
    std::string headerfile(name + ".lst");

    readOrCreateHeaderFile(discFuncSpace, name, headerfile, create);
    this->attributes_.reserve(length);
  }

public:
  /** \brief constructor
     *
     * checks for existence of data path and tries to create it if
     * it doesn't exist.
     *  \param[in] discFuncSpace discrete function space where the functions are defined
     *  \param[in] path          the (absolute or relative) path where this list should be stored
     *  \param[in] name          name for the list, this is also the pathname and the name for the header file
     *  \param[in] create  boolean flag indicating wether the discrete function list shall be created newly.
     */
  DiscreteFunctionList_xdr(const std::string path,
                           const DiscreteFunctionSpaceType& discFuncSpace,
                           const std::string name = "discFuncs",
                           const bool create = false)
    : dataPath_(path + name + "/")
      , discFuncSpace_(discFuncSpace)
      , cleared_(false) {
    std::string headerfile(name + ".lst");

    readOrCreateHeaderFile(discFuncSpace, name, headerfile, create);
  }

  /** \brief constructor
     *
     *  reading a DiscreteFunctionList from a header file
     *  \param[in] discFuncSpace underlying discrete function space
     *  \param[in] name       name of the DiscreteFunctionList, this has to be also the name of the
     *                        folder containing the *.attr and *.dat files and the header file
     *  \param[in] headerfile name of the header file, just the name, no path! (see above)
     *  \param[in] create  boolean flag indicating wether the discrete function list shall be created newly.
     */
  DiscreteFunctionList_xdr(const DiscreteFunctionSpaceType& discFuncSpace,
                           const std::string name, const std::string headerfile,
                           const bool create = false)
    : dataPath_("./" + name + "/")
      , discFuncSpace_(discFuncSpace)
      , cleared_(false) {
    readOrCreateHeaderFile(discFuncSpace, name, headerfile, create);
  }

  /**destructor, close the header file*/
  ~DiscreteFunctionList_xdr() {
    if ( header_.is_open() )
      header_.close();
  }

  /**
     * \brief  store discrete function in list
     *
     *  writes given discrete function to disc and stores given attribute
     *  to the given discrete function.
     *  \param[in]  discFunc    The discrete function to be stored
     *  \param[in]  attr        The parameter to be stored with discrete function,
     *                          name for instance
     *
     *  \return returns the unique index of this function in the list
     */
  inline int push_back( const DiscreteFunctionType& discFunc, AttributeType attr = AttributeType() ) {
    // check if this function list has been cleared
    checkCleared();

    // we started a new run, get a unique index for the filename
    const int currentIndex = nextIndex();

    // transform currentIndex to string and append ".dat" for the datafile
    std::stringstream datafile;
    datafile << dataPath_ << currentIndex << ".dat";

    // append ".attr" to the current index for the attribute file
    std::stringstream attributefile;
    attributefile << dataPath_ << currentIndex << ".att";

    // write discFunc and attribute to disc
    discFunc.write_xdr( datafile.str() );
    attrAct.write( attr, attributefile.str().c_str() );

    // write information about this functions' location to header file
    header_ << currentIndex << ".dat" << std::endl;
    header_.flush();

    std::pair< int, AttributeType > paar = make_pair(currentIndex, attr);
    attributes_.push_back(paar);
    std::stringstream indexAsString;
    indexAsString << currentIndex;
    std::pair< int, std::string > filenamepaar = make_pair( currentIndex, indexAsString.str() );
    filenames_.push_back(filenamepaar);
    return currentIndex;
  } // push_back

  /**
     *  \brief returns the number of stored discrete functions
     *
     *  \return returns the number of stored discrete functions
     */
  inline unsigned int size() const {
    // check if this function list has been cleared
    checkCleared();
    return this->attributes_.size();
  }

private:
  /**
     *  \brief  returns discrete function with unique identifier index.
     *
     *  Assumes, that grid didn't change
     *  since discFunc was written to disk
     *  \param[in]  index index of the desired discrete function
     *  \param[out]  dest result is stored to dest
     *  \return returns true if discFunc was read without errors
     */
  virtual bool getFuncByIndex(const unsigned int index, DiscreteFunctionType& dest) const {
    // check if this function list has been cleared
    checkCleared();
    // typedef
    typedef typename vector< std::pair< int, std::string > >
      ::const_iterator IteratorType;

    std::string datafile = dataPath_;
    // look for matching index in filename vector
    for (IteratorType it = this->filenames_.begin(); it != this->filenames_.end(); ++it)
      if ( (unsigned int) (it->first) == index )
      {
        datafile += it->second;
      }


    datafile += ".dat";

    dest.clear();
    return dest.read_xdr(datafile);
  } // getFuncByIndex

public:
  /**
     *  \brief  returns discrete function number i.
     *
     *  Get the i'th function,
     *  numbering may change during runtime due to deletion of functions.
     *  Assumes, that grid didn't change since discFunc was written to disk.
     *
     *  \param[in]  i number of the desired discrete function
     *  \param[out]  dest result is stored to dest
     *  \return returns true if discFunc was read without errors
     */
  inline bool getFunc(const unsigned int i, DiscreteFunctionType& dest) const {
    // check if this function list has been cleared
    checkCleared();
    if ( i >= filenames_.size() )
      DSC_LOG_ERROR << "Error: i=" << i << " !!!\n";
    assert( i < filenames_.size() );
    std::string filename = filenames_[i].second;

    dest.clear();
    std::string arg = dataPath_;
    arg += filename;
    arg += ".dat";
    return dest.read_xdr(arg);
  } // getFunc

  inline void print(std::ostream& out) const {
    // check if this function list has been cleared
    checkCleared();
    int size = this->size();
    DiscreteFunctionType buffer("buffer", discFuncSpace_);
    for (int i = 0; i != size; ++i)
    {
      this->getFunc(i, buffer);
      buffer.print(out);
    }
  } // print

  /**
     *  \brief sets the i'th function.
     *
     *  Set the i'th function,
     *  numbering may change during runtime due to deletion of functions.
     *  Assumes, that grid didn't change since discFunc was written to disk.
     *
     *  \param[in]  i number of the discrete function that is to be altered
     *  \param[in] arg the function that should replace the i'th function
     *  \return returns true if discFunc was written without errors
     */
  inline bool setFunc(const unsigned int i, DiscreteFunctionType& arg) {
    // check if this function list has been cleared
    checkCleared();
    assert( i >= 0 && i < filenames_.size() );

    std::string filename = filenames_[i].second;

    std::string dest = dataPath_;
    dest += filename;
    dest += ".dat";
    if ( arg.write_xdr(dest) )
    {
      return true;
    }
    return false;
  } // setFunc

private:
  /**
     *  \brief  returns discrete function number i.
     *
     *  Get the i'th function,
     *  numbering may change during runtime due to deletion of functions.
     *  Assumes, that grid didn't change
     *  since discFunc was written to disk.
     *
     *  \param[in]  i number of the desired discrete function
     *  \return returns a pointer to the i'th discrete function
     */
  inline DiscreteFunctionType* getFunc(int i) const {
    // check if this function list has been cleared
    checkCleared();
    DiscreteFunctionType* dest = new DiscreteFunctionType("", discFuncSpace_);
    std::string filename = filenames_[i].second;

    dest->clear();
    std::string arg;
    arg += dataPath_;
    arg += filename;
    arg += ".dat";
    if ( dest->read_xdr(arg) )
    {
      return dest;
    }
    return 0;
  } // getFunc

public:
  /** \brief return the i'th attribute
     *
     *  The i'th attribute is returned,
     *  numbering may change during runtime due to deletion of functions and attributes.
     *
     *  \param[in]  i number of the desired discrete attribute
     *  \param[out]  dest result is stored to dest
     *  \return returns true if attribute was read without errors
     */
  inline bool getAttribute(int i, AttributeType& dest) const {
    // check if this function list has been cleared
    checkCleared();
    std::string filename = dataPath_;
    filename += filenames_[i].second;
    filename += ".att";
    if ( attrAct.read(dest, filename) )
    {
      return true;
    }
    return false;
  } // getAttribute

  /**
     *  \brief  get discrete function by attribute
     *
     *  returns the first discrete function in functions_ with matching parameter, i.e. if attribute is not
     *  unique, you may not get the expected discfunc back! Assumes, that grid didn't change since discFunc
     *  was written. Also, if you didn't push back the functions togther with an attribute, this doesn't make
     *  sense!
     *
     *  \param[in]  attr the attribute of the desired discrete function
     *  \param[out]  dest result is stored to dest
     *  \return returns true if discFunc was read without errors
     */
  inline bool getFuncByAttribute(const AttributeType& attr, DiscreteFunctionType& dest) const {
    // check if this function list has been cleared
    checkCleared();
    // typedef
    typedef typename vector< std::pair< int, AttributeType > >
      ::const_iterator IteratorType;

    // look for matching attribute in attribute vector
    for (IteratorType it = this->attributes_.begin(); it != this->attributes_.end(); ++it)
      if ( (*it).second == attr )
      {
        // it->first==index for correspondig discrete function
        if ( this->getFuncByIndex(it->first, dest) )
          return true;

        return false;
      }


    return false;
  } // getFuncByAttribute

private:
  static int nextIndex(int i = -1) {
    static int id = 0;

    if (i != -1)
      id = i;
    return id++;
  } // nextIndex

public:
  //!return underlying space
  const DiscreteFunctionSpaceType& space() const {
    // check if this function list has been cleared
    checkCleared();
    return discFuncSpace_;
  }

  /** @brief delete all functions in this list
     *
     *  That means, that all functions, managed by this list
     *  are deleted from the hard disk! If a file belonging to
     *  to this list can not be deleted, an error message is
     *  given to std:cerr.
     *
     *  @return returns 0 if list was deleted without errors, -1 if
     *          there where errors.
     */
  int clear() {
    // check if this function list has been cleared
    checkCleared();
    int size = filenames_.size();
    int returnarg = 0;
    // delete a attribute and xdr files
    for (int i = 0; i != size; ++i)
    {
      std::string xdrfile = dataPath_;
      xdrfile += filenames_[i].second;
      std::string attributefile = xdrfile;
      attributefile += ".att";
      xdrfile += ".dat";
      if (remove( attributefile.c_str() ) != 0)
      {
        std::string errormsg;
        errormsg += "Error clearing function list (in file";
        errormsg += filenames_[i].second;
        errormsg += ".att)";
        perror( errormsg.c_str() );
        returnarg = -1;
      }
      if (remove( xdrfile.c_str() ) != 0)
      {
        std::string errormsg;
        errormsg += "Error clearing function list (in file";
        errormsg += filenames_[i].second;
        errormsg += ".dat)";
        perror( errormsg.c_str() );
        returnarg = -1;
      }
    }
    // close header file
    header_.close();
    // now, delete header file
    if (remove( headerfileName_.c_str() ) != 0)
    {
      std::string errormsg;
      errormsg += "Error clearing function list (in file";
      errormsg += headerfileName_;
      errormsg += ")";
      perror( errormsg.c_str() );
      returnarg = -1;
    }
    // finally, delete the containing directory
    // rmdir only deletes the directory if it is empty
    if (rmdir( dataPath_.c_str() ) != 0)
    {
      std::string errormsg;
      errormsg += "Error deleting folder containing function list (";
      errormsg += dataPath_;
      errormsg += ")";
      perror( errormsg.c_str() );
      returnarg = -1;
    }
    cleared_ = true;
    return returnarg;
  } // clear

private:
  void checkCleared() const {
    // check if function list has been cleared
    if (cleared_)
      DUNE_THROW(InvalidStateException,
                 "This discrete function list has been cleared already, further operations on it are prohibited!");
  }

private:
  std::string headerfileName_;
  std::vector< std::pair< int, AttributeType > > attributes_;
  std::vector< std::pair< int, std::string > > filenames_;
  const std::string dataPath_;
  const DiscreteFunctionSpaceType& discFuncSpace_;
  AttributeActions< AttributeType > attrAct;
  std::ofstream header_;
  bool cleared_;
};
} // end of namespace Dune::RB
} // end of namespace Dune

#endif // ifndef DISCFUNCLIST_XDR_HH
