#ifndef DUNE_MODEL_PROBLEM_ALL_HH
#define DUNE_MODEL_PROBLEM_ALL_HH

#include <dune/stuff/configcontainer.hh>

namespace Problem {
    struct Constants {
        const double epsilon;
        const double epsilon_est;
        const double delta;
        Constants( double def_epsilon, double def_epsilon_est, double def_delta )
            : epsilon(Stuff::Config().get("problem.epsilon", def_epsilon))
            , epsilon_est(Stuff::Config().get("problem.epsilon_est", def_epsilon_est))
            , delta(Stuff::Config().get("problem.delta", def_delta))
        {}

        template < typename T, class Validator = Stuff::ValidateAny<T> >
        T get( const std::string& key, const T& def ,
               Stuff::ValidatorInterface< T, Validator > validator = Validator() ) const
        {
            return Stuff::Config().get( std::string("problem.") + key, def, validator);
        }
    };

#define CONSTANTSFUNCTION(d,e,f) \
    static const Constants& constants() {\
    static Constants c(d,e,f); \
    return c; }
}


#endif // DUNE_MODEL_PROBLEM_ALL_HH
