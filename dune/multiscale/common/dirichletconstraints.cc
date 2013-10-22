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

/** Copy the dirichlet values to a given discrete function on the fine mesh.
 *
 *  This method projects the dirichlet values to a function on the coarse mesh.
 *  In a second step, that function is projected to a given function on the fine mesh.
 *
 * @attention The given function will be cleared!
 *
 * @param[in] coarseSpace the discrete function space on the coarse grid, needed
 *                        to perfom the correct projection.
 * @param[in, out] function The function to which the values will be projected.
 */
void copyDirichletValues(const CommonTraits::DiscreteFunctionSpaceType &coarseSpace,
                         CommonTraits::DiscreteFunctionType &function) {
  static bool initialized = false;
  static CommonTraits::DiscreteFunctionType dirichletExtension(
      "Dirichlet Extension", function.space());
  if (!initialized) {
    initialized = true;

    auto dirichletDataPtr = Problem::getDirichletData();
    const auto &dirichletData = *dirichletDataPtr;

    Dune::Fem::GridFunctionAdapter<
        CommonTraits::DirichletDataType,
        CommonTraits::DiscreteFunctionSpaceType::GridPartType> gf(
        "Dirichlet Data", dirichletData, coarseSpace.gridPart());
    CommonTraits::DiscreteFunctionType dirichletExtensionCoarse(
        "Dirichlet Extension Coarse", coarseSpace);
    dirichletExtensionCoarse.clear();
    getConstraintsCoarse(coarseSpace)(gf, dirichletExtensionCoarse);

    static Dune::Stuff::HeterogenousProjection<> projection;
    projection.project(dirichletExtensionCoarse, dirichletExtension);
  }
  function.assign(dirichletExtension);
}


} // namespace Multiscale
} // namespace Dune
