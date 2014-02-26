#ifndef ERROR_CALC_HH
#define ERROR_CALC_HH

#include <dune/multiscale/common/traits.hh>
#include <iosfwd>

namespace Dune {
namespace Multiscale {

class ErrorCalculator {

public:
  ErrorCalculator(const CommonTraits::DiscreteFunctionType* const msfem_solution,
                  const CommonTraits::PdelabVectorType* const fem_solution);

  void print(std::ostream& out);

private:
  const CommonTraits::DiscreteFunctionType* const msfem_solution_;
  const CommonTraits::PdelabVectorType* const fem_solution_;
};

} // namespace Multiscale
} // namespace Dune

#endif // ERROR_CALC_HH
