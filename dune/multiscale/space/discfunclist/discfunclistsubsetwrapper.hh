#ifndef __DISCFUNCLISTSUBSETWRAPPER__
#define __DISCFUNCLISTSUBSETWRAPPER__

#include <iostream>
#include <cmath>

/** \class DiscFuncListSubSetWrapper
   *  \ingroup RB
   *  \brief Wraps a discrete function list, such that only each k-th
   *         function is accessed.
   */
template< class DiscFuncListImp >
class DiscFuncListSubSetWrapper
{
  typedef DiscFuncListImp DiscFuncListType;

public:
  typedef typename DiscFuncListType::DiscreteFunctionType      DiscreteFunctionType;
  typedef typename DiscFuncListType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
  typedef typename DiscFuncListType::AttributeType             AttributeType;

  DiscFuncListSubSetWrapper(DiscFuncListType& list,
                            const unsigned int stepSize)
    : list_(list)
      , stepSize_(stepSize)
  {}

  ~DiscFuncListSubSetWrapper() {}

private:
  template< class DiscreteFunctionType, class AttributeType >
  int push_back( const DiscreteFunctionType& discFunc,
                 AttributeType attr = AttributeType() ) {
    return list_.push_back(discFunc, attr);
  }

public:
  unsigned int size() const {
    return static_cast< unsigned int >
           ( std::floor( double( list_.size() ) / double(stepSize_) ) );
  }

  template< class DiscreteFunctionType >
  bool getFunc(unsigned int i, DiscreteFunctionType& dest) const {
    return list_.getFunc(stepSize_ * i, dest);
  }

  void print(std::ostream& out) const { list_.print(out); }

  template< class DiscreteFunctionType >
  bool setFunc(unsigned int i, DiscreteFunctionType& arg) const {
    return list_.setFunc(stepSize_ * i, arg);
  }

  template< class AttributeType >
  bool getAttribute(int i, AttributeType& dest) const {
    return list_.getAttribute(stepSize_ * i, dest);
  }

  template< class AttributeType, class DiscreteFunctionType >
  bool getFuncByAttribute(AttributeType attr, DiscreteFunctionType& dest) {
    return list_.getFuncByAttribute(attr, dest);
  }

  const DiscreteFunctionSpaceType& space() { return list_.space(); }

  int clear() { return list_.clear(); }

private:
  DiscFuncListType& list_;
  const unsigned int stepSize_;
};

#endif // #ifdef __DISCFUNCLISTSUBSETWRAPPER__
