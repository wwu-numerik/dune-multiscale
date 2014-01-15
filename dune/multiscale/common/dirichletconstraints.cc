#include <config.h>
#include <dune/fem/function/common/gridfunctionadapter.hh>
#include <dune/stuff/discretefunction/projection/heterogenous.hh>
#include <memory>

#include "dirichletconstraints.hh"
#include "dune/multiscale/common/traits.hh"
#include "dune/multiscale/problems/base.hh"
#include "dune/multiscale/problems/selector.hh"

namespace Dune {
namespace Multiscale {

DirichletConstraints<CommonTraits::DiscreteFunctionSpaceType>&
getConstraintsCoarse(const CommonTraits::DiscreteFunctionSpaceType& space) {
  // set dirichlet dofs to zero
  static const auto boundary = Problem::getModelData()->boundaryInfo();
  static DirichletConstraints<CommonTraits::DiscreteFunctionSpaceType> constraints(*boundary, space);
  return constraints;
}

DirichletConstraints<CommonTraits::DiscreteFunctionSpaceType>&
getConstraintsFine(const CommonTraits::DiscreteFunctionSpaceType& space) {
  // set dirichlet dofs to zero
  static const auto boundary = Problem::getModelData()->boundaryInfo();
  static DirichletConstraints<CommonTraits::DiscreteFunctionSpaceType> constraints(*boundary, space);
  return constraints;
}




} // namespace Multiscale
} // namespace Dune
