#ifndef PRINT_INFO_HH
#define PRINT_INFO_HH

#include <dune/multiscale/common/traits.hh>
#include <dune/multiscale/problems/base.hh>
#include <dune/multiscale/problems/selector.hh>
#include <dune/stuff/common/ranges.hh>


namespace Dune {
namespace Multiscale {
namespace FEM {

void print_info(const CommonTraits::ModelProblemDataType& info, std::ostream& out);
void write_discrete_function(typename CommonTraits::DiscreteFunctionType& discrete_solution, const std::string prefix);

template <class DiscreteFunctionType>
void boundaryTreatment(DiscreteFunctionType& rhs) {
  using namespace Dune::Stuff;
  const auto& discreteFunctionSpace = rhs.space();
  static constexpr unsigned int faceCodim = 1;
  for (const auto& entity : discreteFunctionSpace) {
    for (const auto& intersection : DSC::intersectionRange(discreteFunctionSpace.gridPart(), entity)) {
      if (!intersection.boundary())
        continue;
      if (intersection.boundary() && (intersection.boundaryId() != 1))
        continue;

      auto rhsLocal = rhs.localFunction(entity);
      const auto face = intersection.indexInInside();
      for (auto loc_point : DSC::lagrangePointSetRange<faceCodim>(rhs.space(), entity, face))
        rhsLocal[loc_point] = 0;
    }
  }
} // boundaryTreatment

//! set the dirichlet points to the Dirichlet BC
template <class DirichletBC, class DiscreteFunctionType>
void setDirichletValues(DirichletBC& dirichlet_func, DiscreteFunctionType& func) {
  using namespace Dune::Stuff;
  const auto& discreteFunctionSpace = func.space();
  static constexpr unsigned int faceCodim = 1;
  for (const auto& entity : discreteFunctionSpace) {
    for (const auto& intersection : DSC::intersectionRange(discreteFunctionSpace.gridPart(), entity)) {
      if (Dune::Multiscale::Problem::isDirichletBoundary(intersection)) {
        auto funcLocal = func.localFunction(entity);
        const auto face = intersection.indexInInside();
        for (auto loc_point : DSC::lagrangePointSetRange<faceCodim>(func.space(), entity, face)) {
          const auto& global_point =
              entity.geometry().global(discreteFunctionSpace.lagrangePointSet(entity).point(loc_point));
          CommonTraits::RangeType dirichlet_value(0.0);
          dirichlet_func.evaluate(global_point, dirichlet_value);
          funcLocal[loc_point] = dirichlet_value;
        }
      }
    }
  }
} // setDirichletValues

} // namespace FEM {
} // namespace Multiscale {
} // namespace Dune {

#endif // PRINT_INFO_HH
