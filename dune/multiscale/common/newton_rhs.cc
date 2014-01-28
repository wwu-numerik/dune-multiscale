#include <config.h>

#include <dune/stuff/functions/global.hh>
#include <dune/stuff/common/ranges.hh>
#include <dune/stuff/fem/functions/integrals.hh>
#include <dune/multiscale/problems/base.hh>
#include "dune/multiscale/common/traits.hh"
#include "newton_rhs.hh"

void Dune::Multiscale::NewtonRightHandSide::assemble_for_Newton_method(const Dune::Multiscale::CommonTraits::FirstSourceType &f,
                                                                       const Dune::Multiscale::CommonTraits::DiffusionType &A,
                                                                       const Dune::Multiscale::CommonTraits::LowerOrderTermType &F,
                                                                       const Dune::Multiscale::CommonTraits::DiscreteFunctionType &old_u_H,
                                                                       Dune::Multiscale::CommonTraits::DiscreteFunctionType &rhsVector) {
  typedef typename CommonTraits::DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionSpaceType::RangeType RangeType;
  typedef typename DiscreteFunctionSpaceType::JacobianRangeType JacobianRangeType;

  rhsVector.clear();

  for (const auto& entity : rhsVector.space()) {
    const auto& geometry = entity.geometry();
    auto elementOfRHS = rhsVector.localFunction(entity);
    const auto baseSet = rhsVector.space().basisFunctionSet(entity);

    const auto old_u_H_loc = old_u_H.localFunction(entity);
    const auto quadrature = DSFe::make_quadrature(entity, rhsVector.space());

    const auto numDofs = elementOfRHS.numDofs();
    // the return values:
    RangeType f_x;
    std::vector<RangeType> phi_x(numDofs);
    // gradient of base function and gradient of old_u_H
    std::vector<JacobianRangeType> grad_phi_x(numDofs);
    JacobianRangeType grad_old_u_H;
    RangeType value_old_u_H;
    // Let A denote the diffusion operator, then we save
    // A( \gradient old_u_H )
    JacobianRangeType diffusive_flux_in_grad_old_u_H;
    for (auto quadraturePoint : DSC::valueRange(quadrature.nop())) {
      // local (barycentric) coordinates (with respect to entity)
      const auto& local_point = quadrature.point(quadraturePoint);
      const auto global_point = geometry.global(local_point);
      const double det = geometry.integrationElement(local_point);
      // evaluate the Right Hand Side Function f at the current quadrature point and save its value in 'f_y':
      f.evaluate(global_point, f_x);
      // evaluate the current base function at the current quadrature point and save its value in 'z':
      baseSet.evaluateAll(quadrature[quadraturePoint], phi_x); // i = i'te Basisfunktion;
      // evaluate the gradient of the current base function at the current quadrature point and save its value in
      // 'returnGradient':
      baseSet.jacobianAll(quadrature[quadraturePoint], grad_phi_x);
      // get value of old u_H:
      old_u_H_loc.evaluate(quadrature[quadraturePoint], value_old_u_H);
      // get gradient of old u_H:
      old_u_H_loc.jacobian(quadrature[quadraturePoint], grad_old_u_H);
      // evaluate diffusion operator in x(=global_point) and grad_old_u_H
      A.diffusiveFlux(global_point, grad_old_u_H, diffusive_flux_in_grad_old_u_H);

      RangeType F_x;
      F.evaluate(global_point, value_old_u_H, grad_old_u_H, F_x);

      for (int i = 0; i < numDofs; ++i) {
        elementOfRHS[i] += det * quadrature.weight(quadraturePoint) * (f_x * phi_x[i]);
        elementOfRHS[i] -=
            det * quadrature.weight(quadraturePoint) * (diffusive_flux_in_grad_old_u_H[0] * grad_phi_x[i][0]);
        elementOfRHS[i] -= det * quadrature.weight(quadraturePoint) * (F_x * phi_x[i]);
      }
    }
  }
}
