#include <config.h>
// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include <assert.h>
#include <dune/common/exceptions.hh>
#include <dune/stuff/common/parameter/configcontainer.hh>
#include <dune/stuff/fem/localmatrix_proxy.hh>
#include <dune/stuff/fem/functions/integrals.hh>
#include <iostream>
#include <sstream>
#include <vector>

#include "discrete_cell_operator.hh"
#include "dune/multiscale/hmm/cell_problem_solver.hh"

// artificical mass coefficient to guarantee uniqueness and existence of the cell problem solution
// (should be as small as possible)
#define CELL_MASS_WEIGHT 0.0000001

namespace Dune {
namespace Multiscale {
namespace HMM {

void DiscreteCellProblemOperator::operator()(const DiscreteFunction& /* u*/, DiscreteFunction& /*w*/) const {
  std::cout << "the ()-operator of the DiscreteCellProblemOperator class is not yet implemented and still a dummy."
            << std::endl;
  DUNE_THROW(Dune::Exception, "not implemented");
}

//!// x_T is the barycenter of the macro grid element T
void DiscreteCellProblemOperator::assemble_matrix(const DomainType& x_T,
                                                  CellProblemSolver::CellLinearOperatorType& global_matrix) const {
  const double delta = DSC_CONFIG_GET("hmm.delta", 1.0f);

  global_matrix.reserve();
  global_matrix.clear();

  // micro scale base function:
  std::vector<RangeType> phi(periodicDiscreteFunctionSpace_.blockMapper().maxNumDofs());

  // gradient of micro scale base function:
  std::vector<typename BaseFunctionSet::JacobianRangeType> gradient_phi(
      periodicDiscreteFunctionSpace_.blockMapper().maxNumDofs());

  for (const auto& cell_grid_entity : periodicDiscreteFunctionSpace_) {
    const auto& cell_grid_geometry = cell_grid_entity.geometry();
    assert(cell_grid_entity.partitionType() == InteriorEntity);

    DSFe::LocalMatrixProxy<CellProblemSolver::CellLinearOperatorType> local_matrix(global_matrix, cell_grid_entity,
                                                                                   cell_grid_entity);

    const auto& baseSet = local_matrix.domainBasisFunctionSet();
    const auto numBaseFunctions = baseSet.size();

    // for constant diffusion "2*discreteFunctionSpace_.order()" is sufficient, for the general case, it is better to
    // use a higher order quadrature:
    const auto quadrature = DSFe::make_quadrature(cell_grid_entity, periodicDiscreteFunctionSpace_);
    const auto numQuadraturePoints = quadrature.nop();
    for (size_t quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint) {
      // local (barycentric) coordinates (with respect to cell grid entity)
      const auto& local_point = quadrature.point(quadraturePoint);

      // global point in the unit cell Y
      auto global_point = cell_grid_geometry.global(local_point);

      // x_T + (delta * global_point)
      DomainType x_T_delta_global_point;
      for (int k = 0; k < dimension; ++k) {
        x_T_delta_global_point[k] = x_T[k] + (delta * global_point[k]);
      }

      const double weight = quadrature.weight(quadraturePoint) * cell_grid_geometry.integrationElement(local_point);

      baseSet.jacobianAll(quadrature[quadraturePoint], gradient_phi);
      baseSet.evaluateAll(quadrature[quadraturePoint], phi);

      for (unsigned int i = 0; i < numBaseFunctions; ++i) {
        // A( x_T + \delta y, \nabla \phi )
        // diffusion operator evaluated in (x_T + \delta y , \nabla \phi)
        JacobianRangeType diffusion_in_gradient_phi;
        diffusion_operator_.diffusiveFlux(x_T_delta_global_point, gradient_phi[i], diffusion_in_gradient_phi);
        for (unsigned int j = 0; j < numBaseFunctions; ++j) {
          // stiffness contribution
          local_matrix.add(j, i, weight * (diffusion_in_gradient_phi[0] * gradient_phi[j][0]));

          // mass contribution
          local_matrix.add(j, i, CELL_MASS_WEIGHT * weight * (phi[i][0] * phi[j][0]));
        }
      }
    }
  }
} // assemble_matrix

void DiscreteCellProblemOperator::assemble_jacobian_matrix(
    const DomainType& x_T, const JacobianRangeType& grad_coarse_function, const DiscreteFunction& old_fine_function,
    CellProblemSolver::CellLinearOperatorType& global_matrix) const {
  const double delta = DSC_CONFIG_GET("hmm.delta", 1.0f);

  global_matrix.reserve();
  global_matrix.clear();

  // micro scale base function:
  std::vector<RangeType> phi(periodicDiscreteFunctionSpace_.blockMapper().maxNumDofs());

  // gradient of micro scale base function:
  std::vector<typename BaseFunctionSet::JacobianRangeType> gradient_phi(
      periodicDiscreteFunctionSpace_.blockMapper().maxNumDofs());

  for (const auto& cell_grid_entity : periodicDiscreteFunctionSpace_) {
    const auto& cell_grid_geometry = cell_grid_entity.geometry();
    assert(cell_grid_entity.partitionType() == InteriorEntity);

    DSFe::LocalMatrixProxy<CellProblemSolver::CellLinearOperatorType> local_matrix(global_matrix, cell_grid_entity,
                                                                                   cell_grid_entity);
    auto local_fine_function = old_fine_function.local_function(cell_grid_entity);

    const auto& baseSet = local_matrix.domainBasisFunctionSet();
    const auto numBaseFunctions = baseSet.size();

    // for constant diffusion "2*periodicDiscreteFunctionSpace_.order()" is sufficient, for the general case, it is
    // better to use a higher order quadrature:
    const auto quadrature = DSFe::make_quadrature(cell_grid_entity, periodicDiscreteFunctionSpace_);
    const auto numQuadraturePoints = quadrature.nop();
    for (size_t quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint) {
      // local (barycentric) coordinates (with respect to entity)
      const auto& local_point = quadrature.point(quadraturePoint);

      // global point in the unit cell Y
      const auto global_point = cell_grid_geometry.global(local_point);

      // x_T + (delta * global_point)
      DomainType x_T_delta_global_point;
      for (int k = 0; k < dimension; ++k) {
        x_T_delta_global_point[k] = x_T[k] + (delta * global_point[k]);
      }

      const double weight = quadrature.weight(quadraturePoint) * cell_grid_geometry.integrationElement(local_point);

      baseSet.jacobianAll(quadrature[quadraturePoint], gradient_phi);
      baseSet.evaluateAll(quadrature[quadraturePoint], phi);

      for (unsigned int i = 0; i < numBaseFunctions; ++i) {
        // JA( x_T + \delta y, \nabla_x PHI_H(x_T) + \nabla_y old_fine_function ) \nabla \phi
        // jacobian matrix of the diffusion operator evaluated in (x_T + \delta y , \nabla_x PHI_H(x_T) + \nabla_y
        // old_fine_function ) in direction \nabla \phi
        typename BaseFunctionSet::JacobianRangeType grad_local_fine_function;
        local_fine_function.jacobian(quadrature[quadraturePoint], grad_local_fine_function);
        // here: no multiplication with jacobian inverse transposed required!

        typename BaseFunctionSet::JacobianRangeType position_vector;
        for (int k = 0; k < dimension; ++k) {
          position_vector[0][k] = grad_coarse_function[0][k] + grad_local_fine_function[0][k];
        }

        // jacobian of diffusion operator evaluated in (x,grad coarse + grad fine) in direction of the gradient of the
        // current base function
        JacobianRangeType jac_diffusion_flux;
        diffusion_operator_.jacobianDiffusiveFlux(x_T_delta_global_point, position_vector, gradient_phi[i],
                                                  jac_diffusion_flux);

        for (unsigned int j = 0; j < numBaseFunctions; ++j) {
          // stiffness contribution
          local_matrix.add(j, i, weight * (jac_diffusion_flux[0] * gradient_phi[j][0]));

          // mass contribution
          local_matrix.add(j, i, CELL_MASS_WEIGHT * weight * (phi[i][0] * phi[j][0]));
        }
      }
    }
  }
} // assemble_jacobian_matrix

void DiscreteCellProblemOperator::assembleCellRHS_linear(const DomainType& x_T, const JacobianRangeType& gradient_PHI_H,
                                                         DiscreteFunction& cell_problem_RHS) const {
  const auto& discreteFunctionSpace = cell_problem_RHS.space();

  // set entries to zero:
  cell_problem_RHS.clear();

  // get edge length of cell:
  const double delta = DSC_CONFIG_GET("hmm.delta", 1.0f);

  // gradient of micro scale base function:
  std::vector<JacobianRangeType> gradient_phi(discreteFunctionSpace.blockMapper().maxNumDofs());

  for (const auto& cell_grid_entity : periodicDiscreteFunctionSpace_) {
    const auto& geometry = cell_grid_entity.geometry();
    assert(cell_grid_entity.partitionType() == InteriorEntity);

    auto elementOfRHS = cell_problem_RHS.local_function(cell_grid_entity);

    const auto& baseSet = elementOfRHS.basisFunctionSet();
    const auto numBaseFunctions = baseSet.size();

    const auto quadrature = DSFe::make_quadrature(cell_grid_entity, discreteFunctionSpace);
    const auto numQuadraturePoints = quadrature.nop();
    for (size_t quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint) {
      const auto& local_point = quadrature.point(quadraturePoint);

      // global point in the unit cell Y:
      const auto global_point = geometry.global(local_point);

      // x_T + (delta * global_point)
      DomainType x_T_delta_global_point;
      for (int k = 0; k < dimension; ++k) {
        x_T_delta_global_point[k] = x_T[k] + (delta * global_point[k]);
      }

      const double weight = quadrature.weight(quadraturePoint) * geometry.integrationElement(local_point);

      // A^{\eps}( x_T + \delta y) \nabla_x PHI_H(x_T)
      // diffusion operator evaluated in (x_T + \delta y) multiplied with \nabla_x PHI_H(x_T)
      JacobianRangeType diffusion_in_gradient_PHI_H;
      diffusion_operator_.diffusiveFlux(x_T_delta_global_point, gradient_PHI_H, diffusion_in_gradient_PHI_H);

      baseSet.jacobianAll(quadrature[quadraturePoint], gradient_phi);

      for (unsigned int i = 0; i < numBaseFunctions; ++i) {
        elementOfRHS[i] -= weight * (diffusion_in_gradient_PHI_H[0] * gradient_phi[i][0]);
      }
    }
  }
} // assembleCellRHS_linear

void DiscreteCellProblemOperator::assembleCellRHS_nonlinear(const DomainType& x_T,
                                                            const JacobianRangeType& grad_coarse_function,
                                                            const DiscreteFunction& old_fine_function,
                                                            DiscreteFunction& cell_problem_RHS) const {
  const auto& discreteFunctionSpace = cell_problem_RHS.space();

  // set entries to zero:
  cell_problem_RHS.clear();
  // get edge length of cell:
  const double delta = DSC_CONFIG_GET("hmm.delta", 1.0f);

  // gradient of micro scale base function:
  std::vector<JacobianRangeType> gradient_phi(discreteFunctionSpace.blockMapper().maxNumDofs());

  for (const auto& cell_grid_entity : discreteFunctionSpace) {
    const auto& geometry = cell_grid_entity.geometry();
    assert(cell_grid_entity.partitionType() == InteriorEntity);

    auto local_old_fine_function = old_fine_function.local_function(cell_grid_entity);
    auto elementOfRHS = cell_problem_RHS.local_function(cell_grid_entity);

    const auto& baseSet = elementOfRHS.basisFunctionSet();
    const auto numBaseFunctions = baseSet.size();

    const auto quadrature = DSFe::make_quadrature(cell_grid_entity, discreteFunctionSpace);
    const auto numQuadraturePoints = quadrature.nop();
    for (size_t quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint) {
      const auto& local_point = quadrature.point(quadraturePoint);

      // global point in the unit cell Y:
      const auto global_point = geometry.global(local_point);

      // x_T + (delta * global_point)
      DomainType x_T_delta_global_point;
      for (int k = 0; k < dimension; ++k) {
        x_T_delta_global_point[k] = x_T[k] + (delta * global_point[k]);
      }

      JacobianRangeType grad_old_fine_function;
      local_old_fine_function.jacobian(quadrature[quadraturePoint], grad_old_fine_function);
      // here: no multiplication with jacobian inverse transposed required!

      JacobianRangeType position_vector;
      for (int k = 0; k < dimension; ++k) {
        position_vector[0][k] = grad_coarse_function[0][k] + grad_old_fine_function[0][k];
      }

      // A^{\eps}( x_T + \delta y, \nabla_x grad_coarse_function(x_T) + \nabla_y old_fine_function )
      JacobianRangeType diffusive_flux;
      diffusion_operator_.diffusiveFlux(x_T_delta_global_point, position_vector, diffusive_flux);

      const double weight = quadrature.weight(quadraturePoint) * geometry.integrationElement(local_point);
      baseSet.jacobianAll(quadrature[quadraturePoint], gradient_phi);
      for (unsigned int i = 0; i < numBaseFunctions; ++i) {
        elementOfRHS[i] -= weight * (diffusive_flux[0] * gradient_phi[i][0]);
      }
    }
  }
} // assembleCellRHS_nonlinear

void DiscreteCellProblemOperator::assemble_jacobian_corrector_cell_prob_RHS( // the global quadrature point in the macro
                                                                             // grid element T
    const DomainType& x_T,                                                   // barycenter of macro entity T
    // gradient of the old coarse function (old means last iteration step)
    const JacobianRangeType& grad_old_coarse_function, // \nabla_x u_H^{(n-1)}(x_T)
    // gradient of the corrector of the old coarse function
    const DiscreteFunction& corrector_of_old_coarse_function, /*Q_h(u_H^{(n-1)})*/
    // gradient of the current macroscopic base function
    const JacobianRangeType& grad_coarse_base_function, // \nabla_x \Phi_H(x_T)
    // rhs cell problem:
    DiscreteFunction& jac_corrector_cell_problem_RHS) const {
  const auto& discreteFunctionSpace = corrector_of_old_coarse_function.space();

  // set entries of right hand side to zero:
  jac_corrector_cell_problem_RHS.clear();
  const double delta = DSC_CONFIG_GET("hmm.delta", 1.0f);

  // gradient of micro scale base function:
  std::vector<JacobianRangeType> gradient_phi(discreteFunctionSpace.blockMapper().maxNumDofs());

  for (const auto& cell_grid_entity : discreteFunctionSpace) {
    const auto& geometry = cell_grid_entity.geometry();
    assert(cell_grid_entity.partitionType() == InteriorEntity);

    // local Q_h(u_H^{(n-1)}):
    const auto local_Q_old_u_H = corrector_of_old_coarse_function.local_function(cell_grid_entity);
    auto elementOfRHS = jac_corrector_cell_problem_RHS.local_function(cell_grid_entity);

    const auto& baseSet = elementOfRHS.basisFunctionSet();
    const auto numBaseFunctions = baseSet.size();

    const auto quadrature = DSFe::make_quadrature(cell_grid_entity, discreteFunctionSpace);
    const auto numQuadraturePoints = quadrature.nop();
    for (size_t quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint) {
      const auto& local_point = quadrature.point(quadraturePoint);

      // global point in the unit cell Y:
      const auto global_point_in_Y = geometry.global(local_point);

      // x_T + (delta * global_point_in_Y)
      DomainType x_T_plus_delta_y;
      for (int k = 0; k < dimension; ++k) {
        x_T_plus_delta_y[k] = x_T[k] + (delta * global_point_in_Y[k]);
      }

      // grad_y Q_h( u_H^{(n-1)})
      JacobianRangeType grad_Q_old_u_H;
      local_Q_old_u_H.jacobian(quadrature[quadraturePoint], grad_Q_old_u_H);
      // here: no multiplication with jacobian inverse transposed required!

      JacobianRangeType position_vector;
      for (int k = 0; k < dimension; ++k) {
        position_vector[0][k] = grad_old_coarse_function[0][k] + grad_Q_old_u_H[0][k];
      }

      JacobianRangeType direction_vector;
      for (int k = 0; k < dimension; ++k) {
        direction_vector[0][k] = grad_coarse_base_function[0][k];
      }

      // DA^{\eps}( x_T + \delta y, \nabla_x u_H^{(n-1)})(x_T) + \nabla_y Q_h( u_H^{(n-1)})(y) )( \nabla_x \Phi_H(x_T) )
      JacobianRangeType jacobian_diffusive_flux;
      diffusion_operator_.jacobianDiffusiveFlux(x_T_plus_delta_y, position_vector, direction_vector,
                                                jacobian_diffusive_flux);

      const double weight = quadrature.weight(quadraturePoint) * geometry.integrationElement(local_point);
      baseSet.jacobianAll(quadrature[quadraturePoint], gradient_phi);

      for (unsigned int i = 0; i < numBaseFunctions; ++i) {
        elementOfRHS[i] -= weight * (jacobian_diffusive_flux[0] * gradient_phi[i][0]);
      }
    }
  }
} // assemble_jacobian_corrector_cell_prob_RHS

} // namespace HMM {
} // namespace Multiscale {
} // namespace Dune {
