#ifndef ERROR_CALC_HH
#define ERROR_CALC_HH

#include <dune/multiscale/common/traits.hh>
#include <dune/multiscale/msfem/localsolution_proxy.hh>
#include <iosfwd>
#include <memory>

namespace Dune {
namespace Multiscale {

namespace MsFEM {
class LocalsolutionProxy;
}

class ErrorCalculator {

public:
  ErrorCalculator(const std::unique_ptr<MsFEM::LocalsolutionProxy>& msfem_solution,
                  const CommonTraits::ConstDiscreteFunctionType* const fem_solution);

  std::map<std::string, double> print(std::ostream& out);

private:
  const std::unique_ptr<MsFEM::LocalsolutionProxy>& msfem_solution_;
  const CommonTraits::ConstDiscreteFunctionType* const fem_solution_;
};

} // namespace Multiscale
} // namespace Dune

#endif // ERROR_CALC_HH
