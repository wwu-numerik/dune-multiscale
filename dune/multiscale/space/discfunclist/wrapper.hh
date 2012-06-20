#ifndef __DISCFUNCLISTWRAPPER_HH__
#define __DISCFUNCLISTWRAPPER_HH__

// includes
#include "discfunclist_xdr.hh"
#include "discfunclist_mem.hh"

namespace Dune {
namespace Multiscale {
/** @class DiscreteFunctionList_Wrapper
   *  @ingroup DiscFuncList
   *  @brief Wrap a given discrete function list
   *
   *  This is a wrapper for a discrete function list. It gives the
   *  possibility to keep a given number of functions in memory
   *  while the rest is stored on hard disk.
   */
template< class DiscreteFunctionListTraits >
class DiscreteFunctionList_Wrapper
  : public DiscreteFunctionList_xdr< DiscreteFunctionListTraits >
{
  typedef typename DiscreteFunctionListTraits::DiscreteFunctionType DiscreteFunctionType;
  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType  DiscreteFunctionSpaceType;
  typedef DiscreteFunctionList_xdr< DiscreteFunctionListTraits >    BaseType;

public:
  /** @brief constructor
     *
     *
     */
  DiscreteFunctionList_Wrapper(const int& blocksize,
                               const DiscreteFunctionSpaceType& discFuncSpace,
                               std::string name = "discFuncs")
    : DiscreteFunctionList_xdr< DiscreteFunctionListTraits >(discFuncSpace, name)
      , blocksize_(blocksize)
      , range_( std::make_pair< int, int >(0, 0) ) {
// buffer_.reserve(blocksize_);
// for (int i=0; i!=blocksize_; ++i)
// buffer_[i] = std::make_pair<int, DiscreteFunctionType>(0, DiscreteFunctionType("noName", BaseType::discFuncSpace_));
  }

  /** @brief constructor
     *
     *
     */
  DiscreteFunctionList_Wrapper(const int& blocksize,
                               const DiscreteFunctionSpaceType& discFuncSpace,
                               int length,
                               std::string name = "discFuncs")
    : DiscreteFunctionList_xdr< DiscreteFunctionListTraits >(discFuncSpace, length, name)
      , blocksize_(blocksize)
      , range_( std::make_pair< int, int >(0, 0) ) {
// buffer_.reserve(blocksize_);
// for (int i=0; i!=blocksize_; ++i)
// buffer_[i] = std::pair<int, DiscreteFunctionType>(0, DiscreteFunctionType("noName", BaseType::discFuncSpace_));
  }

  /** @brief constructor
     *
     *
     */
  DiscreteFunctionList_Wrapper(const int& blocksize,
                               const DiscreteFunctionSpaceType& discFuncSpace,
                               const std::string name,
                               const std::string headerfile)
    : DiscreteFunctionList_xdr< DiscreteFunctionListTraits >(discFuncSpace, name, headerfile)
      , blocksize_(blocksize)
      , range_( std::make_pair< int, int >(0, 0) ) {
// buffer_.reserve(blocksize_);
// for (int i=0; i!=blocksize_; ++i)
// buffer_[i] = std::pair<int, DiscreteFunctionType>(0, DiscreteFunctionType("noName", BaseType::discFuncSpace_));
  }

  // !destructor, free all memory allocated by the list
  ~DiscreteFunctionList_Wrapper() {
    typedef typename std::map< int, DiscreteFunctionType* >::iterator IteratorType;
    for (IteratorType it = buffer_.begin(); it != buffer_.end(); ++it)
      delete it->second;
  }

  bool getFuncSingle(const int& i, DiscreteFunctionType& discFunc) {
    return BaseType::getFunc(i, discFunc);
  }

  /** @brief get the i'th discrete function
     *
     *
     */
  inline DiscreteFunctionType& getFunc(unsigned int i) {
    // check if i is in the range of functions currently stored in buffer_
    if ( (i >= range_.first) && (i < range_.second) )
      return *buffer_[i];
    else {
      // get the next blocksize_-functions
      int begin = std::floor(i / blocksize_);
      begin *= blocksize_;
      unsigned int end = begin + blocksize_;
      end = std::min( end, BaseType::size() );
      range_.first = begin;
      range_.second = end;
      typedef typename std::map< int, DiscreteFunctionType* >::iterator IteratorType;
      for (IteratorType it = buffer_.begin(); it != buffer_.end(); ++it)
        if (it->second != NULL)
          delete it->second;
      buffer_.clear();
      for (int func = begin; func != end; ++func)
      {
        buffer_[func] = new DiscreteFunctionType("noName", BaseType::discFuncSpace_);
        BaseType::getFunc(func, *buffer_[func]);
      }
    }
    return *buffer_[i];
  } // getFunc

private:
  const int blocksize_;
  std::map< int, DiscreteFunctionType* > buffer_;
  std::map< int, int > indexMapping_;
  std::pair< int, int > range_;
};

template< class DiscreteFunctionListImp >
class FuncListBlock
{
  typedef DiscreteFunctionListImp                             DiscreteFunctionList;
  typedef typename DiscreteFunctionList::DiscreteFunctionType DiscreteFunctionType;
  typedef typename std::vector< DiscreteFunctionType* >
    ::iterator VecIterator;
  typedef typename DiscreteFunctionType
    ::DiscreteFunctionSpaceType DiscreteFunctionSpace;

public:
  FuncListBlock(DiscreteFunctionList& list, unsigned int size)
    : list_(list)
      , size_(size)
      , discFuncSpace_( list_.space() ) {
    funcs_.resize(size_);
    VecIterator end = funcs_.end();
    for (VecIterator it = funcs_.begin(); it != end; ++it)
      *it = new DiscreteFunctionType("noName", discFuncSpace_);
  }

  ~FuncListBlock() {
    VecIterator end = funcs_.end();

    for (VecIterator it = funcs_.begin(); it != end; ++it)
      delete *it;
  }

private:
  FuncListBlock(FuncListBlock& other) {}

public:
  unsigned int first() const {
    assert( map_.begin() != map_.end() );
    return map_.begin()->second;
  }

  unsigned int last() const {
    assert( map_.begin() != map_.end() );
    return map_.rbegin()->second;
  }

  /** @brief load functions into this block
     *
     * Loads the functions given by the range [begin,end] from hard
     * disk.
     *
     * @param[in] begin index in discfunclist for the first function
     * @param[in] end   index in discfunclist for one past the last function
     */
  void load(const unsigned int begin, const unsigned int end) {
    VecIterator it = funcs_.begin();

    map_.clear();
    unsigned int vecPos = 0;
    unsigned int realEnd = std::min( end, list_.size() );
    for (unsigned int i = begin; i < realEnd; ++i, ++it, ++vecPos)
    {
      list_.getFunc(i, **it);
      map_[vecPos] = i;
    }
  } // load

  /** @brief move to next block
     *
     * @return true if the block was move, false if it was not moved
     *         because the block already show the end of the list.
     */
  bool nextBlock() {
    if (map_.rbegin()->second != list_.size() - 1)
    {
      int currentEnd = map_.rbegin()->second;
      load( currentEnd + 1, std::min(list_.size(), currentEnd + size_ + 1) );
      return true;
    } else
      return false;
  } // nextBlock

  /** @brief map local function index to global one, that is, the
     *         the functions position in the function list.
     */
  unsigned int mapToGlobal(const unsigned int i) const {
    std::map< unsigned int, unsigned int >::const_iterator ret = map_.find(i);
    assert( ret != map_.end() );

    return ret->second;
  }

  /** @brief return how many functions are really stored in this block.
     */
  unsigned int usedSize() const {
    return map_.size();
  }

  /** @brief get func
     *
     * @param[in] i the LOCAL index of the function
     */
  DiscreteFunctionType& getFunc(const unsigned int i) {
    assert( i < funcs_.size() );
    assert(funcs_[i] != NULL);
    return *funcs_[i];
  }

private:
  void clearFuncs() {
    VecIterator end = funcs_.end();

    for (VecIterator it = funcs_.begin();
         it != end; ++it)
    {
      delete *it;
      *it = NULL;
    }
  } // clearFuncs

private:
  // the function list wrapped by this class
  DiscreteFunctionList& list_;
  // size of this block of functions
  const unsigned int size_;
  // Pointers for the stored function
  std::vector< DiscreteFunctionType* > funcs_;
  // map local indices to global ones
  std::map< unsigned int, unsigned int > map_;
  const DiscreteFunctionSpace& discFuncSpace_;
};
} // end of namespace Dune::RB
} // end of namespace Dune

#endif // ifndef __DISCFUNCLISTWRAPPER_HH__
