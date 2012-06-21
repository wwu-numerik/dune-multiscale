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
    };

#define CONSTANTSFUNCTION(d,f,e) \
    static const Constants& constants() {\
    static Constants c(d,f,e); \
    return c; }
}


#endif // DUNE_MODEL_PROBLEM_ALL_HH
