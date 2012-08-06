#ifndef DISCFUNCLIST_MEM_HH
#define DISCFUNCLIST_MEM_HH

#include <vector>
#include <sstream>
#include "discfunclistinterface.hh"

namespace Dune {
namespace Multiscale {
using std::vector;
using std::pair;
using std::make_pair;

/** \class   DiscreteFunctionList_mem
   *  \ingroup DiscFuncList
   *  \brief   Stores functions together with attributes in a list
   *
   *  The function is stored in memory together with a given
   *  attribute using a std::vector
   *
   *  \param DiscreteFunctionListTraits a traits class that defines the
   *                                    the type of the discrete function
   */
template< class DiscreteFunctionListTraits >
class DiscreteFunctionList_mem
  : public DiscreteFunctionListInterface< DiscreteFunctionListTraits >
{
public:
  // typedefs
  typedef DiscreteFunctionListTraits            Traits;
  typedef typename Traits::DiscreteFunctionType DiscreteFunctionType;

  typedef typename DiscreteFunctionType
    ::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
  // see interface header for description
  typedef typename SelectTypeAttributeTypeIfEnabled
  < DiscreteFunctionListTraits, int >::Type AttributeType;

private:
  typedef DiscreteFunctionList_mem< DiscreteFunctionListTraits > ThisType;
  typedef DiscreteFunctionListInterface
  < DiscreteFunctionListTraits > BaseType;

public:
  //! constructor taking discFuncSpace as argument
  DiscreteFunctionList_mem(const DiscreteFunctionSpaceType& discFuncSpace)
    : current_index_(-1)
      , discFuncSpace_(discFuncSpace) {}

  /**
     *  \brief constructor
     *
     *  by giving the constructor an approximate number of functions that will
     *  be stored, you can avoid repeated reallocating of memory, which will make
     *  storing functions quicker
     *  \param[in] discFuncSpace Discrete function space for the functions stored by
     *             this list
     *  \param[in] length number of functions to be stored
     */
  DiscreteFunctionList_mem(const DiscreteFunctionSpaceType& discFuncSpace, int length)
    : discFuncSpace_(discFuncSpace) {
    this->functions_.reserve(length);
    current_index_ = -1;
  }

  //!destructor, free all memory allocated by the list
  ~DiscreteFunctionList_mem() {
    typedef typename vector< std::pair< DiscreteFunctionType*, AttributeType > >::iterator IteratorType;
    for (IteratorType it = functions_.begin(); it != functions_.end(); ++it)
      delete it->first;
  }

  /**
     * \brief  store discrete function in list
     *
     *  stores given discrete function together with attribute in new memory
     *  \param[in]  discFunc    The discrete function to be stored
     *  \param[in]  attr        The parameter to be stored with discrete function,
     *                          name for instance
     *  \return returns unique index connected to stored discrete function
     */
  inline int push_back( const DiscreteFunctionType& discFunc, AttributeType attr = AttributeType() ) {
    // we started a new run
    ++current_index_;
    // copy discrete function to new memory and get a pointer to it
    DiscreteFunctionType* discFuncPtr = new DiscreteFunctionType(discFunc);
    // pair discrete function and attribute
    std::pair< DiscreteFunctionType*, AttributeType > paar = make_pair(discFuncPtr, attr);
    functions_.push_back(paar);
    // store index together with position of function pointer to be able to find function pointer in the sequel
    map_index_.push_back( make_pair( ( functions_.size() ) - 1, current_index_ ) );
    return current_index_;
  } // push_back

  /**
     *  \brief returns the number of stored discrete functions
     */
  inline unsigned int size() const {
    return this->functions_.size();
  }

private:
  /**
     * \brief  returns discrete function with identifier index.
     *
     * Assumes, that grid didn't change
     * since discFunc was written to disk
     * \param[in]  index index of the desired discrete function
     * \param[out]  dest result is stored to dest
     * \return returns true if discFunc was read without errors
     */
  inline bool getFuncByIndex(const unsigned int index, DiscreteFunctionType& dest) const {
    typedef vector< pair< int, int > >::const_iterator IteratorType;
    for (IteratorType it = this->map_index_.begin(); it != this->map_index_.end(); ++it)
      if ( (unsigned int) (it->second) == index )
      {
        dest.assign(*this->functions_[it->first].first);
        return true;
      }


    return false;
  } // getFuncByIndex

public:
  /**
     * \brief  get i'th discrete function.
     *
     * Assumes, that grid didn't change
     * since discFunc was written to disk
     * \param[in]  i desired discrete function
     * \param[out]  dest result is stored to dest
     */
  inline bool getFunc(const unsigned int i, DiscreteFunctionType& dest) const {
    assert( i < this->size() );
    dest.assign(*functions_[i].first);
    return true;
  }

private:
  /**
     * \brief  returns discrete function number i.
     *
     * Get the i'th function,
     * numbering may change during runtime due to deletion of functions.
     * Assumes, that grid didn't change
     * since discFunc was written to disk.
     * \param[in]  i number of the desired discrete function
     * \return returns a pointer to the i'th discrete function
     */
  inline const DiscreteFunctionType* getFunc(const unsigned int i) const {
    typedef vector< pair< int, int > >::const_iterator IteratorType;
    for (IteratorType it = this->map_index_.begin(); it != this->map_index_.end(); ++it)
      if (it->second == i)
      {
        return this->functions_[it->first].first;
      }


    return 0;
  } // getFunc

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
  inline bool getFuncByAttribute(const AttributeType& attr, DiscreteFunctionType& dest) const {
    // get iterator type for functions_-vector
    typedef typename vector< std::pair< DiscreteFunctionType*,
                                        AttributeType > >
      ::const_iterator IteratorType;

    for (IteratorType it = this->functions_.begin(); it != this->functions_.end(); ++it)
      if (it->second == attr)
      {
        dest.assign(*it->first);
        return true;
      }


    return false;
  } // getFuncByAttribute

  /**
     *  \brief sets the i'th function.
     *
     *  Set the i'th function,
     *  numbering may change during runtime due to deletion of functions.
     *  Assumes, that grid didn't change since discFunc was written to disk.
     *
     *  \param[in] i    number of the discrete function that is to be altered
     *  \param[in] arg  the function that should replace the i'th function
     *  \param[in] attr an attribute that should be stored together with the function
     *  \return         returns true if discFunc was written without errors
     */
  inline bool setFunc( const unsigned int i, const DiscreteFunctionType& arg, const AttributeType attr
                         = AttributeType() ) {
    {
      // copy discrete function to new memory and get a pointer to it
      DiscreteFunctionType* discFuncPtr = new DiscreteFunctionType(arg);
      // pair discrete function and attribute
      std::pair< DiscreteFunctionType*, AttributeType > paar = make_pair(discFuncPtr, attr);
      if ( i < functions_.size() )
      {
        delete functions_[i].first;
      }
      functions_[i] = paar;
      return true;
    }
  } // setFunc

  /**
     * \brief deletes the discfunc with given index
     * \param[in] i index of the function to be deleted
     * \return      returns true if function was deleted successfully,
     *              false otherwise
     */
  inline bool deleteFunc(const unsigned int i) {
// typedef vector<pair<int,int> >::iterator IteratorType;
// for (IteratorType it=this->map_index_.begin(); it!=this->map_index_.end(); ++it)
// if (it->second==index)
// {

    if (this->size() > i)
    {
      delete this->functions_[i].first;
      this->functions_[i].erase();
      return true;
    }
    return false;
  } // deleteFunc

  //!return underlying space
  const DiscreteFunctionSpaceType& space() const {
    return discFuncSpace_;
  }

protected:
  vector< pair< DiscreteFunctionType*, AttributeType > > functions_;
  // first int is number in vec functions_, second is index
  vector< pair< int, int > > map_index_;

private:
  int current_index_;
  const DiscreteFunctionSpaceType& discFuncSpace_;
};
} // end of namespace Dune::RB
} // end of namespace Dune

#endif // ifndef DISCFUNCLIST_MEM_HH
