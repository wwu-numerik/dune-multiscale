#ifndef DUNE_MODEL_PROBLEM_ALL_HH
#define DUNE_MODEL_PROBLEM_ALL_HH

#include <dune/stuff/common/parameter/configcontainer.hh>

namespace Problem {
struct Constants
{
  const double epsilon;
  const double epsilon_est;
  const double delta;
  Constants(double def_epsilon, double def_epsilon_est, double def_delta)
    : epsilon( Dune::Stuff::Common::Parameter::Config().get("problem.epsilon", def_epsilon) )
      , epsilon_est( Dune::Stuff::Common::Parameter::Config().get("hmm.epsilon_guess", def_epsilon_est) )
      , delta( Dune::Stuff::Common::Parameter::Config().get("hmm.delta", def_delta) )
  {}

  template< typename T, class Validator = Dune::Stuff::Common::Parameter::ValidateAny< T > >
  T get( const std::string& key, const T& def,
         Dune::Stuff::Common::Parameter::ValidatorInterface< T, Validator > validator = Validator() ) const {
    return Dune::Stuff::Common::Parameter::Config().get(std::string("problem.") + key, def, validator);
  }
};

#define CONSTANTSFUNCTION(d, e, f) \
  static const Constants &constants() \
  { \
    static Constants c(d, e, f); \
    return c; \
  }
} // namespace Problem

#endif // DUNE_MODEL_PROBLEM_ALL_HH
