#ifndef  __UTILITY_HH__
#define  __UTILITY_HH__

namespace Dune {
namespace Multiscale {
#define CHECK_FOR_ATTRIBUTE(__AttributeName__) \
  template< typename T > \
  struct HasType ## __AttributeName__ \
  { \
    template< typename TType > \
    static TType makeT();  \
    template< typename Attribute > \
    struct YesType  \
    { char small; }; \
    struct NoType \
    { char big[2]; }; \
    static NoType test(...); \
 \
    template< typename TType > \
    static YesType< typename TType::__AttributeName__ > test(TType); \
 \
    static const bool Yes = ( sizeof( test( makeT< T >() ) ) == sizeof(YesType< T >) ); \
  }; \
\
  template< typename Class, typename DefaultType > \
  struct SelectType ## __AttributeName__ ## IfEnabled \
  { \
    template< bool, typename TType > \
    struct Select \
    { \
      typedef DefaultType Type; \
    }; \
  \
    template< typename TType > \
    struct Select< true, TType > \
    { \
      typedef typename TType::__AttributeName__ Type; \
    }; \
  \
    typedef typename Select< HasType ## __AttributeName__< Class > \
                               ::Yes, Class >::Type Type; \
  };

// ! helper function for template meta programs.
// !
// ! During compile time, the compiler see the return value as an object of type
// ! T. However, there is never code generated for this object, because there is
// ! no implementation for this function.
// !
// ! \tparam T type of object to be created for compiler.
// ! \return object of type T, that is never initialized.
template< typename T >
static T makeT();
template< typename T >
static T& makeRT();
} // end namespace RB
} // end namespace Dune

#endif  /*__UTILITY_HH__*/
