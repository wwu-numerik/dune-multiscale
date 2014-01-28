#include <config.h>
// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "localoperator.hh"

#include <assert.h>
#include <boost/assert.hpp>
#include <dune/common/exceptions.hh>
#include <dune/multiscale/problems/selector.hh>
#include <dune/stuff/common/parameter/configcontainer.hh>
#include <dune/stuff/fem/matrix_object.hh>
#include <dune/stuff/fem/localmatrix_proxy.hh>
#include <dune/stuff/discretefunction/projection/heterogenous.hh>
#include <dune/stuff/fem/functions/integrals.hh>
#include <dune/multiscale/tools/misc.hh>
#include <dune/multiscale/msfem/localproblems/localproblemsolver.hh>


namespace Dune {
namespace Multiscale {
namespace MsFEM {

LocalProblemOperator::LocalProblemOperator(const CoarseSpaceType& coarse_space, const LocalGridDiscreteFunctionSpaceType& subDiscreteFunctionSpace,
                                           const DiffusionOperatorType& diffusion_op)
  : subDiscreteFunctionSpace_(subDiscreteFunctionSpace)
  , diffusion_operator_(diffusion_op)
  , coarse_space_(coarse_space)
{}

// is a given 'point' in the convex hull of corner 0, corner 1 and corner 2 (which forms a codim 0 entity)
bool LocalProblemOperator::point_is_in_element(const DomainType& corner_0, const DomainType& corner_1,
                                               const DomainType& corner_2, const DomainType& point) const {
  DomainType v = corner_0 - corner_2;
  DomainType w = corner_1 - corner_2;
  DomainType p = point - corner_2;

  double lambda_0, lambda_1;

  if (v[0] != 0.0) {
    double val0 = p[1] - ((v[1] / v[0]) * p[0]);
    double val1 = w[1] - ((v[1] / v[0]) * w[0]);

    if (val1 != 0.0) {
      lambda_1 = val0 / val1;
      lambda_0 = (p[0] - (lambda_1 * w[0])) / v[0];
      if ((0.0 <= lambda_0) && (1.0 >= lambda_0) && (0.0 <= lambda_1) && (1.0 >= lambda_1) &&
          (lambda_0 + lambda_1 <= 1.0))
        return true;
      else
        return false;
    } else {
      DUNE_THROW(Dune::InvalidStateException,
                 "... in method 'point_is_in_element': Given corners do not span a codim 0 entity in 2D.");
    }
  } else {
    if ((w[0] != 0.0) && (v[1] != 0.0)) {
      lambda_1 = p[0] / w[0];
      lambda_0 = (p[1] - (lambda_1 * w[1])) / v[1];
      if ((0.0 <= lambda_0) && (1.0 >= lambda_0) && (0.0 <= lambda_1) && (1.0 >= lambda_1) &&
          (lambda_0 + lambda_1 <= 1.0))
        return true;
      else
        return false;
    } else {
      DUNE_THROW(Dune::InvalidStateException,
                 "... in method 'point_is_in_element': Given corners do not span a codim 0 entity in 2D.");
    }
  }
}

//! stiffness matrix for a linear elliptic diffusion operator
// for oversampling strategy 1 (no constraints)
void LocalProblemOperator::assemble_matrix(LocalProblemSolver::LocProbLinearOperatorTypeType& global_matrix) const
    // x_T is the barycenter of the macro grid element T
{
  global_matrix.reserve(DSFe::diagonalAndNeighborStencil(global_matrix));
  global_matrix.clear();

  // local grid basis functions:
  std::vector<RangeType> phi(subDiscreteFunctionSpace_.mapper().maxNumDofs());

  // gradient of micro scale base function:
  std::vector<typename BasisFunctionSetType::JacobianRangeType> gradient_phi(
      subDiscreteFunctionSpace_.mapper().maxNumDofs());
  typename BasisFunctionSetType::JacobianRangeType diffusion_in_gradient_phi;

  for (const auto& sub_grid_entity : subDiscreteFunctionSpace_) {
    const auto& sub_grid_geometry = sub_grid_entity.geometry();

    DSFe::LocalMatrixProxy<LocalProblemSolver::LocProbLinearOperatorTypeType> local_matrix(
        global_matrix, sub_grid_entity, sub_grid_entity);

    const auto& baseSet = local_matrix.domainBasisFunctionSet();
    const auto numBaseFunctions = baseSet.size();

    // for constant diffusion "2*discreteFunctionSpace_.order()" is sufficient, for the general case, it is better to
    // use a higher order quadrature:
    const auto quadrature = DSFe::make_quadrature(sub_grid_entity, subDiscreteFunctionSpace_);

    const auto numQuadraturePoints = quadrature.nop();
    for (size_t quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint) {
      // local (barycentric) coordinates (with respect to local grid entity)
      const auto& local_point = quadrature.point(quadraturePoint);
      const auto global_point = sub_grid_geometry.global(local_point);

      const double weight = quadrature.weight(quadraturePoint) * sub_grid_geometry.integrationElement(local_point);

      baseSet.jacobianAll(quadrature[quadraturePoint], gradient_phi);
      baseSet.evaluateAll(quadrature[quadraturePoint], phi);

      for (unsigned int i = 0; i < numBaseFunctions; ++i) {
        // A( x, \nabla \phi(x) )
        diffusion_operator_.diffusiveFlux(global_point, gradient_phi[i], diffusion_in_gradient_phi);
        for (unsigned int j = 0; j < numBaseFunctions; ++j) {
          // stiffness contribution
          local_matrix.add(j, i, weight * (diffusion_in_gradient_phi[0] * gradient_phi[j][0]));
          // mass contribution (just for stabilization!)
          // local_matrix.add( j, i, 0.00000001 * weight * (phi[ i ][ 0 ] * phi[ j ][ 0 ]) );
        }
      }
    }
  }
  global_matrix.communicate();
} // assemble_matrix

//! stiffness matrix for a linear elliptic diffusion operator
void LocalProblemOperator::assemble_matrix(LocalProblemSolver::LocProbLinearOperatorTypeType& global_matrix,
                                           const LocalGridList::CoarseNodeVectorType& coarse_node_vector) const
    // x_T is the barycenter of the macro grid element T
{
  global_matrix.reserve(DSFe::diagonalAndNeighborStencil(global_matrix));
  global_matrix.clear();

  // local grid basis functions:
  std::vector<RangeType> phi(subDiscreteFunctionSpace_.mapper().maxNumDofs());

  // gradient of micro scale base function:
  std::vector<typename BasisFunctionSetType::JacobianRangeType> gradient_phi(
      subDiscreteFunctionSpace_.mapper().maxNumDofs());

  for (const auto& sub_grid_entity : subDiscreteFunctionSpace_) {
    const auto& sub_grid_geometry = sub_grid_entity.geometry();
    assert(sub_grid_entity.partitionType() == InteriorEntity);

    std::vector<int> sub_grid_entity_corner_is_relevant;
    for (int c = 0; c < sub_grid_geometry.corners(); ++c) {
      for (size_t coarse_node_local_id = 0; coarse_node_local_id < coarse_node_vector.size(); ++coarse_node_local_id) {
        // if the subgrid corner 'c' is in the 'relevant coarse node vector' and if 'c' was not yet added to the
        // vector 'sub_grid_entity_corner_is_relevant' then add it to the vector
        if ((coarse_node_vector[coarse_node_local_id] == sub_grid_geometry.corner(c)) &&
            (std::find(sub_grid_entity_corner_is_relevant.begin(), sub_grid_entity_corner_is_relevant.end(), c) ==
             sub_grid_entity_corner_is_relevant.end())) {
          sub_grid_entity_corner_is_relevant.push_back(c);
        }
      }
    }

    DSFe::LocalMatrixProxy<LocalProblemSolver::LocProbLinearOperatorTypeType> local_matrix(
        global_matrix, sub_grid_entity, sub_grid_entity);

    const auto& baseSet = local_matrix.domainBasisFunctionSet();
    const auto numBaseFunctions = baseSet.size();

    std::vector<RangeType> value_phi(numBaseFunctions);
    // for constant diffusion "2*discreteFunctionSpace_.order()" is sufficient, for the general case, it is better to
    // use a higher order quadrature:
    const auto quadrature = DSFe::make_quadrature(sub_grid_entity, subDiscreteFunctionSpace_);
    const auto numQuadraturePoints = quadrature.nop();
    for (size_t quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint) {
      // local (barycentric) coordinates (with respect to local grid entity)
      const auto& local_point = quadrature.point(quadraturePoint);
      const auto global_point = sub_grid_geometry.global(local_point);

      const double weight = quadrature.weight(quadraturePoint) * sub_grid_geometry.integrationElement(local_point);

      baseSet.jacobianAll(quadrature[quadraturePoint], gradient_phi);
      baseSet.evaluateAll(quadrature[quadraturePoint], phi);

      for (size_t sgec = 0; sgec < sub_grid_entity_corner_is_relevant.size(); ++sgec) {
        baseSet.evaluateAll(sub_grid_geometry.local(sub_grid_geometry.corner(sub_grid_entity_corner_is_relevant[sgec])),
                            value_phi);
        for (unsigned int i = 0; i < numBaseFunctions; ++i) {
          if (value_phi[i] == 1.0) {
            assert(dimension == 2);
            phi[i][0] = 0.0;
            gradient_phi[i][0][0] = 0.0;
            gradient_phi[i][0][1] = 0.0;
          }
        }
      }

      for (unsigned int i = 0; i < numBaseFunctions; ++i) {
        // A( x, \nabla \phi(x) )
        typename BasisFunctionSetType::JacobianRangeType diffusion_in_gradient_phi;
        diffusion_operator_.diffusiveFlux(global_point, gradient_phi[i], diffusion_in_gradient_phi);
        for (unsigned int j = 0; j < numBaseFunctions; ++j) {
          // stiffness contribution
          local_matrix.add(j, i, weight * (diffusion_in_gradient_phi[0] * gradient_phi[j][0]));
          // mass contribution (just for stabilization!)
          // local_matrix.add( j, i, 0.00000001 * weight * (phi[ i ][ 0 ] * phi[ j ][ 0 ]) );
        }
      }
    }
  }
  global_matrix.communicate();
} // assemble_matrix

void LocalProblemOperator::set_zero_boundary_condition_RHS(const LocalGridDiscreteFunctionSpaceType& local_space,
                                                           LocalGridDiscreteFunctionType& rhs) const {
  const auto& discreteFunctionSpace = rhs.space();
  const auto& localGridPart = local_space.gridPart();

  // set Dirichlet Boundary to zero
  for (const auto& localgrid_entity : discreteFunctionSpace) {
    const auto& host_entity = localgrid_entity;

    auto iit = localGridPart.ibegin(host_entity);
    const auto endiit = localGridPart.iend(host_entity);
    for (; iit != endiit; ++iit) {
      if (iit->neighbor() && iit->boundaryId() != 1) {
        continue;
      }

      auto rhs_local = rhs.localFunction(localgrid_entity);
      const auto& lagrangePointSet = discreteFunctionSpace.lagrangePointSet(localgrid_entity);

      const int face = (*iit).indexInInside();
      for (const auto& lp : DSC::lagrangePointSetRange(lagrangePointSet, face)) {
        rhs_local[lp] = 0;
      }
    }
  }
} // end method

double LocalProblemOperator::normRHS(const LocalProblemOperator::LocalGridDiscreteFunctionType& rhs) const {
  const auto& discreteFunctionSpace = rhs.space();

  double norm = 0.0;
  for (const auto& entity : discreteFunctionSpace) {
    // create quadrature for given geometry type
    const auto quadrature = DSFe::make_quadrature(entity, discreteFunctionSpace);
    const auto& geo = entity.geometry();
    const auto localRHS = rhs.localFunction(entity);

    // integrate
    for (auto quadraturePoint : DSC::valueRange(quadrature.nop())) {
      const double weight =
          quadrature.weight(quadraturePoint) * geo.integrationElement(quadrature.point(quadraturePoint));
      RangeType value(0.0);
      localRHS.evaluate(quadrature[quadraturePoint], value);

      norm += weight * value * value;
    }
  }
  return norm;
} // end method

void LocalProblemOperator::assemble_local_RHS(const JacobianRangeType& e,
                                              LocalProblemOperator::LocalGridDiscreteFunctionType& local_problem_RHS) const {
  const auto& discreteFunctionSpace = local_problem_RHS.space();
  local_problem_RHS.clear();
  // gradient of micro scale base function:
  std::vector<JacobianRangeType> gradient_phi(discreteFunctionSpace.mapper().maxNumDofs());

  for (const auto& local_grid_entity : discreteFunctionSpace) {
    const auto& geometry = local_grid_entity.geometry();
    auto elementOfRHS = local_problem_RHS.localFunction(local_grid_entity);
    const auto& baseSet = elementOfRHS.basisFunctionSet();
    const auto numBaseFunctions = baseSet.size();

    const auto quadrature = DSFe::make_quadrature(local_grid_entity, discreteFunctionSpace);
    const auto numQuadraturePoints = quadrature.nop();
    for (size_t quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint) {
      const auto& local_point = quadrature.point(quadraturePoint);
      // remember, we are concerned with: - \int_{U(T)} (A^eps)(x) e · ∇ \phi(x)
      // global point in the subgrid
      const auto global_point = geometry.global(local_point);
      const double weight = quadrature.weight(quadraturePoint) * geometry.integrationElement(local_point);

      // A^eps(x) e
      // diffusion operator evaluated in 'x' multiplied with e
      JacobianRangeType diffusion_in_e;
      diffusion_operator_.diffusiveFlux(global_point, e, diffusion_in_e);

      baseSet.jacobianAll(quadrature[quadraturePoint], gradient_phi);
      for (unsigned int i = 0; i < numBaseFunctions; ++i) {
        elementOfRHS[i] -= weight * (diffusion_in_e[0] * gradient_phi[i][0]);
      }
    }
  }
} // assemble_local_RHS


void LocalProblemOperator::assembleAllLocalRHS(const CoarseEntityType& coarseEntity,
                                               MsFEMTraits::LocalSolutionVectorType& allLocalRHS) const {
  BOOST_ASSERT_MSG(allLocalRHS.size() > 0, "You need to preallocate the necessary space outside this function!");

  //! @todo correct the error message below (+1 for simplecial, +2 for arbitrary), as there's no finespace any longer
//  BOOST_ASSERT_MSG(
//      (DSG::is_simplex_grid(coarse_space_) && allLocalRHS.size() == GridType::dimension + 1) ||
//          (!(DSG::is_simplex_grid(coarse_space_)) &&
//           static_cast<long long>(allLocalRHS.size()) ==
//               static_cast<long long>(specifier.fineSpace().mapper().maxNumDofs() + 2)),
//      "You need to allocate storage space for the correctors for all unit vector/all coarse basis functions"
//      " and the dirichlet- and neuman corrector");

  // build unit vectors (needed for cases where rhs is assembled for unit vectors instead of coarse
  // base functions)
  JacobianRangeType unitVectors[dimension];
  for (int i = 0; i < dimension; ++i) {
    for (int j = 0; j < dimension; ++j) {
      if (i == j) {
        unitVectors[i][0][j] = 1.0;
      } else {
        unitVectors[i][0][j] = 0.0;
      }
    }
  }

  // get dirichlet and neumann data
  auto neumannDataPtr = Dune::Multiscale::Problem::getNeumannData();
  const auto& neumannData = *neumannDataPtr;
  const auto& discreteFunctionSpace = allLocalRHS[0]->space();

  //! @todo we should use the dirichlet constraints here somehow
  LocalGridDiscreteFunctionType dirichletExtension("dirichletExtension", discreteFunctionSpace);
  dirichletExtension.clear();
  CommonTraits::DiscreteFunctionType dirichletExtensionCoarse("Dirichlet Extension Coarse", coarse_space_);
  dirichletExtensionCoarse.clear();
  //! @todo is this needed or could it be replaced by a method from dirichletconstraints.hh?
  projectDirichletValues(dirichletExtensionCoarse);
  Dune::Stuff::HeterogenousProjection<> projection;
  projection.project(dirichletExtensionCoarse, dirichletExtension);

  // set entries to zero:
  for (auto& rhs : allLocalRHS)
    rhs->clear();

  // get the base function set of the coarse space for the given coarse entity
  const auto& coarseBaseSet = coarse_space_.basisFunctionSet(coarseEntity);
  std::vector<CoarseBaseFunctionSetType::JacobianRangeType> coarseBaseFuncJacs(coarseBaseSet.size());

  // gradient of micro scale base function:
  std::vector<JacobianRangeType> gradient_phi(discreteFunctionSpace.mapper().maxNumDofs());
  std::vector<RangeType> phi(discreteFunctionSpace.mapper().maxNumDofs());

  const auto numBoundaryCorrectors = DSG::is_simplex_grid(coarse_space_) ? 1u : 2u;
  const auto numInnerCorrectors = allLocalRHS.size() - numBoundaryCorrectors;

  for (auto& localGridCell : discreteFunctionSpace) {
    const auto& geometry = localGridCell.geometry();
    const bool hasBoundaryIntersection = localGridCell.hasBoundaryIntersections();
    auto dirichletLF = dirichletExtension.localFunction(localGridCell);
    JacobianRangeType dirichletJac(0.0);

    for (std::size_t coarseBaseFunc = 0; coarseBaseFunc < allLocalRHS.size(); ++coarseBaseFunc) {
      auto rhsLocalFunction = allLocalRHS[coarseBaseFunc]->localFunction(localGridCell);

      const auto& baseSet = rhsLocalFunction.basisFunctionSet();
      const auto numBaseFunctions = baseSet.size();

      // correctors with index < numInnerCorrectors are for the basis functions, corrector at
      // position numInnerCorrectors is for the neumann values, corrector at position numInnerCorrectors+1
      // for the dirichlet values.
      if (coarseBaseFunc < numInnerCorrectors || coarseBaseFunc == numInnerCorrectors + 1) {
        const auto quadrature = DSFe::make_quadrature(localGridCell, discreteFunctionSpace);
        const auto numQuadraturePoints = quadrature.nop();
        for (size_t quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint) {
          const auto& local_point = quadrature.point(quadraturePoint);
          // global point in the subgrid
          const auto global_point = geometry.global(local_point);

          const double weight = quadrature.weight(quadraturePoint) * geometry.integrationElement(local_point);

          JacobianRangeType diffusion(0.0);
          if (coarseBaseFunc < numInnerCorrectors) {
            if (DSG::is_simplex_grid(coarse_space_))
              diffusion_operator_.diffusiveFlux(global_point, unitVectors[coarseBaseFunc], diffusion);
            else {
              const DomainType quadInCoarseLocal = coarseEntity.geometry().local(global_point);
              coarseBaseSet.jacobianAll(quadInCoarseLocal, coarseBaseFuncJacs);
              diffusion_operator_.diffusiveFlux(global_point, coarseBaseFuncJacs[coarseBaseFunc], diffusion);
            }
          } else {
            dirichletLF.jacobian(local_point, dirichletJac);
            diffusion_operator_.diffusiveFlux(global_point, dirichletJac, diffusion);
          }
          baseSet.jacobianAll(quadrature[quadraturePoint], gradient_phi);
          for (unsigned int i = 0; i < numBaseFunctions; ++i) {
            rhsLocalFunction[i] -= weight * (diffusion[0] * gradient_phi[i][0]);
          }
        }
      }

      // boundary integrals
      if (coarseBaseFunc == numInnerCorrectors && hasBoundaryIntersection) {
        const auto intEnd = discreteFunctionSpace.gridPart().iend(localGridCell);
        for (auto iIt = discreteFunctionSpace.gridPart().ibegin(localGridCell); iIt != intEnd; ++iIt) {
          const auto& intersection = *iIt;
          if (DMP::is_neumann(intersection)) {
            const auto orderOfIntegrand = (CommonTraits::polynomial_order - 1) + 2 * (CommonTraits::polynomial_order + 1);
            const auto quadOrder = std::ceil((orderOfIntegrand + 1) / 2);
            const auto faceQuad = DSFe::make_quadrature(intersection, discreteFunctionSpace, quadOrder);
            RangeType neumannValue(0.0);
            const auto numQuadPoints = faceQuad.nop();
            // loop over all quadrature points
            for (unsigned int iqP = 0; iqP < numQuadPoints; ++iqP) {
              // get local coordinate of quadrature point
              const auto& xLocal = faceQuad.localPoint(iqP);
              const auto& faceGeometry = intersection.geometry();

              // the following does not work because subgrid does not implement geometryInInside()
              // const auto& insideGeometry    = intersection.geometryInInside();
              // const typename FaceQuadratureType::CoordinateType& xInInside = insideGeometry.global(xLocal);
              // therefore, we have to do stupid things:
              const auto& xGlobal = faceGeometry.global(xLocal);
              auto insidePtr = intersection.inside();
              const auto& insideEntity = *insidePtr;
              const auto& xInInside = insideEntity.geometry().local(xGlobal);
              const double factor = faceGeometry.integrationElement(xLocal) * faceQuad.weight(iqP);

              neumannData.evaluate(xGlobal, neumannValue);
              baseSet.evaluateAll(xInInside, phi);
              for (unsigned int i = 0; i < numBaseFunctions; ++i) {
                rhsLocalFunction[i] -= factor * (neumannValue * phi[i]);
              }
            }
          }
        }
      }
    }
  }

  return;
}

void LocalProblemOperator::assemble_local_RHS(
    const JacobianRangeType& e,
    const LocalGridList::CoarseNodeVectorType& coarse_node_vector,
    LocalGridDiscreteFunctionType& local_problem_RHS) const {

  const auto& discreteFunctionSpace = local_problem_RHS.space();
  local_problem_RHS.clear();

  // gradient of micro scale base function:
  std::vector<JacobianRangeType> gradient_phi(discreteFunctionSpace.mapper().maxNumDofs());

  for (const auto& local_grid_entity : discreteFunctionSpace) {
    const auto& geometry = local_grid_entity.geometry();
    assert(local_grid_entity.partitionType() == InteriorEntity);

    // if the we use oversampling strategy 2 or 3/Lagrange, we need to sort out some coarse grid nodes:
    std::vector<int> sub_grid_entity_corner_is_relevant;
    for (int c = 0; c < geometry.corners(); ++c) {
      for (size_t coarse_node_local_id = 0; coarse_node_local_id < coarse_node_vector.size();
           ++coarse_node_local_id) {
        // if the subgrid corner 'c' is in the 'relevant coarse node vector' and if 'c' was not yet added to the
        // vector 'sub_grid_entity_corner_is_relevant' then add it to the vector
        if ((coarse_node_vector[coarse_node_local_id] == geometry.corner(c)) &&
            (std::find(sub_grid_entity_corner_is_relevant.begin(), sub_grid_entity_corner_is_relevant.end(), c) ==
             sub_grid_entity_corner_is_relevant.end())) {
          sub_grid_entity_corner_is_relevant.push_back(c);
          // std :: cout << std ::endl << "geometry.corner(" << c << ") = " << geometry.corner(c) << " is relevant."
          // << std ::endl;
        }
      }
    }

    auto elementOfRHS = local_problem_RHS.localFunction(local_grid_entity);

    const auto& baseSet = elementOfRHS.basisFunctionSet();
    const auto numBaseFunctions = baseSet.size();
    const auto quadrature = DSFe::make_quadrature(local_grid_entity, discreteFunctionSpace);
    const auto numQuadraturePoints = quadrature.nop();
    for (size_t quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint) {
      const auto& local_point = quadrature.point(quadraturePoint);

      // remember, we are concerned with: - \int_{U(T)} (A^eps)(x) e · ∇ \phi(x)
      const auto global_point = geometry.global(local_point);
      const double weight = quadrature.weight(quadraturePoint) * geometry.integrationElement(local_point);

      // A^eps(x) e
      // diffusion operator evaluated in 'x' multiplied with e
      JacobianRangeType diffusion_in_e;
      diffusion_operator_.diffusiveFlux(global_point, e, diffusion_in_e);
      baseSet.jacobianAll(quadrature[quadraturePoint], gradient_phi);
      std::vector<std::vector<RangeType>> phi_values(sub_grid_entity_corner_is_relevant.size());
      for (auto j : DSC::valueRange(sub_grid_entity_corner_is_relevant.size())) {
        baseSet.evaluateAll(geometry.local(geometry.corner(sub_grid_entity_corner_is_relevant[j])), phi_values[j]);
      }

      for (unsigned int i = 0; i < numBaseFunctions; ++i) {
        bool zero_entry = false;
        for (size_t sgec = 0; sgec < sub_grid_entity_corner_is_relevant.size(); ++sgec) {
          const auto& value_phi_i = phi_values[sgec][i];
          if (value_phi_i == 1.0) {
            zero_entry = true;
          }
        }
        if (!zero_entry)
          elementOfRHS[i] -= weight * (diffusion_in_e[0] * gradient_phi[i][0]);
        else
          elementOfRHS[i] = 0.0;
      }
    }
  }
} // assemble_local_RHS

/** Set the dirichlet values to a given discrete function on the sub mesh
*
* @param[in, out] function The function in which the values will be set.
*/
void LocalProblemOperator::projectDirichletValues(CommonTraits::DiscreteFunctionType &function) const {
  /*  // make sure, we are on a hexahedral element
    BOOST_ASSERT_MSG(function.space().gridPart().grid().leafIndexSet().geomTypes(0).size()==1 &&
           function.space().gridPart().grid().leafIndexSet().geomTypes(0)[0].isCube(),
           "This method only works for hexahedral elements at the moment!");*/

  const auto& gridPart = function.space().gridPart();
  auto dirichletDataPtr = Multiscale::Problem::getDirichletData();
  const auto& dirichletData = *dirichletDataPtr;
  //  Fem::GridFunctionAdapter<Multiscale::Problem::DirichletDataBase,
  //  typename LocalGridDiscreteFunctionType::LocalGridDiscreteFunctionSpaceType::GridPartType> gf("dirichlet", dirichletData ,
  // gridPart);
  RangeType dirichletVal(0.0);
  for (const auto& localCell : function.space()) {
    if (localCell.hasBoundaryIntersections())
      for (const auto& intersection : DSC::intersectionRange(gridPart, localCell)) {
        if (DMP::is_dirichlet(intersection)) {
          auto funcLocal = function.localFunction(localCell);
          const auto& lagrangePointSet = function.space().lagrangePointSet(localCell);
          const auto faceNumber = intersection.indexInInside();
          for (auto lp : DSC::lagrangePointSetRange<1>(function.space(), localCell, faceNumber)) {
            auto lagrangePoint = lagrangePointSet.point(lp);
            auto lagrangePointGlobal = localCell.geometry().global(lagrangePoint);
            dirichletData.evaluate(lagrangePointGlobal, dirichletVal);
            funcLocal[lp] = dirichletVal;
          }
        }
      }
  }
  return;
}

} // namespace MsFEM {
} // namespace Multiscale {
} // namespace Dune {
