#ifndef ERROR_CALC_HH
#define ERROR_CALC_HH

#include <dune/multiscale/common/traits.hh>

namespace Dune {
namespace Multiscale {

class ErrorCalculator {

public:
    ErrorCalculator(const CommonTraits::DiscreteFunctionType* msfem_solution,
                    const CommonTraits::DiscreteFunctionType* fem_solution);

    void print(std::ostream& out);
private:
    const CommonTraits::DiscreteFunctionType* msfem_solution_;
    const CommonTraits::DiscreteFunctionType* fem_solution_;
};

} // namespace Multiscale
} // namespace Dune

#endif // ERROR_CALC_HH