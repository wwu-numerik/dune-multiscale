#ifndef DISCFUNCLIST_INTERFACE
#define DISCFUNCLIST_INTERFACE

#include <dune/common/bartonnackmanifcheck.hh>
#include "utility.hh"
#include <iostream>
#include <fstream>
#include <vector>

namespace Dune {
namespace Multiscale {
// the following macro adds structures
// HasTypeAttributeType and SelectTypeAttributeTypeIfEnabled
CHECK_FOR_ATTRIBUTE(AttributeType)

/** \class   AttributeActions
   *  \ingroup DiscFuncList
   *  \brief   Write attribute to disk and read it from disk
   *
   *  The AttributeActions template struct provides functionality to write an attribute to disk
   *  and read it from disk. This is the standart implementation for all AttributeTypes that
   *  can be written to disk using std::ofstream and the write-memberfunction of std::ofstream.
   *  If you like to use another type you have to specialize this class for your attribute type
   *
   *  \param AttributeType the type of the attribute to be handled by this
   *                       class
   */
template< typename AttributeType >
struct AttributeActions
{
  /**
     *  \brief write attribute to disk
     *
     *  \param[in] arg   attribute to be written to disk
     *  \param[in] file  filename where attribute should be written
     *
     *  \return returns true if attribute was read without errors
     */
  bool write(AttributeType& arg, const std::string& file) {
    std::ofstream foo( file.c_str() );

    foo.precision(10);
    // convert AttributeType to c string
    // std::ostringstream bar;
    // bar << arg;
    // const char * bar_c_str = bar.str().c_str();
    // //write arg to disk
    // foo.write(bar_c_str, strlen(bar_c_str));
    foo << arg;
    return true;
  } // write

  /**
     *  \brief read attribute from disk
     *
     *  \param[out] dest destination where attribute should be
     *                   written after beeing read from disk
     *  \param[in]  file filename to be used for reading from disk
     *
     *  \return          returns true if attribute was read without
     *                   errors
     *
     */
  bool read(AttributeType& dest, const std::string& file) {
    std::ifstream foo( file.c_str() );

    if (!foo)
    {
      DUNE_THROW(IOError, "Could not open " << file << " for reading!");
    }
    foo.precision(10);
    // char * bar;
    // //TODO: use istringstream
    // foo.read(bar, sizeof (T));
    foo >> dest;
    return true;
  } // read
};

/** \class   DiscreteFunctionListInterface
   *  \ingroup DiscFuncList
   *  \brief   Stores functions together with attributes in a list
   *
   *  You have to provide a class DiscreteFunctionListTraits that provides
   *  a type definition
   *  DiscreteFunctionType
   *  for the type of the discrete functions to be stored by this list.
   *  You can also store an attribute together with the function. For this you have to
   *  define AttributeType in DiscreteFunctionListTraits, e.g.
   *    typedef double AttributeType;
   *  if you want to store time information together with the function.
   *
   *  \param DiscreteFunctionListTraits a traits class that defines the
   *                                    discrete function type
   *  \param DiscreteFunctionListImp    derived class (crtp)
   *
   */
template< class DiscreteFunctionListTraits >
class DiscreteFunctionListInterface
{
public:
  // typedefs
  typedef DiscreteFunctionListTraits            Traits;
  typedef typename Traits::DiscreteFunctionType DiscreteFunctionType;
  typedef typename DiscreteFunctionType
    ::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
  typedef DiscreteFunctionListInterface
  < DiscreteFunctionListTraits > ThisType;
  /*the next line determines whether AttributeType is defined in Traits class, if it is not defined,
     * int is used as default*/
  typedef typename SelectTypeAttributeTypeIfEnabled
  < DiscreteFunctionListTraits, int >::Type AttributeType;

public:
  virtual ~DiscreteFunctionListInterface() {}

  /**
     * \brief  store discrete function in list
     *
     *  writes given discrete function to disc and stores given attribute
     *  to the given discrete function.
     *  \param[in]  discFunc    The discrete function to be stored
     *  \param[in]  attr        The parameter to be stored with discrete function,
     *                          name for instance
     *  \return returns unique id connected to stored discrete function
     */
  virtual int push_back( const DiscreteFunctionType& discFunc, AttributeType attr = AttributeType() ) = 0;

  /**
     * \brief append a whole DiscreteFunctionList
     *
     * inserts the given DiscreteFunctionList at the end
     * \param[in] other the other DiscreteFunctionList
     */
  void push_back_list(const ThisType& other) {
    const int otherSize = other.size();
    DiscreteFunctionType temp( "temp", this->space() );

    for (int i = 0; i != otherSize; ++i)
    {
      temp.clear();
      other.getFunc(i, temp);
      this->push_back(temp);
    }
  } // push_back_list

  /**
     *  \brief returns the number of stored discrete functions
     *
     *  \return returns the number of stored discrete functions
     */
  virtual unsigned int size() const = 0;

  virtual int clear() {
    return 0;
  }

private:
  /**
     * \brief  returns discrete function with identifier id.
     *
     * Assumes, that grid didn't change
     * since discFunc was written to disk
     * \param[in]   index  index of the desired discrete function
     * \param[out]  dest   result is stored to dest
     * \return returns true if discFunc was read without errors
     */
  virtual bool getFuncByIndex(const unsigned int index, DiscreteFunctionType& dest) const = 0;

public:
  /**
     * \brief  get discrete function by attribute
     *
     * returns the first discrete function in functions_ with matching parameter, i.e. if attribute is not
     * unique, you may not get the expected discfunc back! Assumes, that grid didn't change since discFunc
     * was written.
     * \param[in]  attr the attribute of the desired discrete function
     * \param[out]  dest result is stored to dest
     * \return returns true if discFunc was read without errors
     */
  virtual bool getFuncByAttribute(const AttributeType& attr, DiscreteFunctionType& dest) const = 0;

private:
  virtual bool getFunc(const unsigned int i, DiscreteFunctionType& dest) const = 0;

public:
  //!return underlying space
  virtual const DiscreteFunctionSpaceType& space() const = 0;
};
} // end of namespace Dune::RB
} // end of namespace Dune

#endif // ifndef DISCFUNCLIST_INTERFACE
/*vim: set et sw=2:*/
