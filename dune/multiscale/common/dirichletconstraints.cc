#include "dirichletconstraints.hh"

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

/** Project the dirichlet values to a given discrete function on the  mesh
*
* @param[in, out] function The function to which the values will be projected.
*/
void projectDirichletValues(const CommonTraits::DiscreteFunctionSpaceType& coarseSpace,
                            CommonTraits::DiscreteFunctionType& function) {
  static bool initialized = false;
  static CommonTraits::DiscreteFunctionType dirichletExtension("Dirichlet Extension", function.space());
  if (!initialized) {
    initialized = true;

    auto dirichletDataPtr = Problem::getDirichletData();
    const auto& dirichletData = *dirichletDataPtr;

    Dune::Fem::GridFunctionAdapter<CommonTraits::DirichletDataType,
                                   CommonTraits::DiscreteFunctionSpaceType::GridPartType> gf("Dirichlet Data",
                                                                                             dirichletData,
                                                                                             coarseSpace.gridPart());
    static CommonTraits::DiscreteFunctionType dirichletExtensionCoarse("Dirichlet Extension Coarse", coarseSpace);
    dirichletExtensionCoarse.clear();
    getConstraintsCoarse(coarseSpace)(gf, dirichletExtensionCoarse);

    static Dune::Stuff::HeterogenousProjection<> projection;
    projection.project(dirichletExtensionCoarse, dirichletExtension);
  }
  getConstraintsFine(function.space())(dirichletExtension, function);
}

} // namespace Multiscale
} // namespace Dune
