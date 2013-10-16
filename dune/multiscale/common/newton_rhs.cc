#include "newton_rhs.hh"

void Dune::Multiscale::NewtonRightHandSide::assemble_for_Newton_method(const Dune::Multiscale::CommonTraits::FirstSourceType &f, const Dune::Multiscale::CommonTraits::DiffusionType &A, const Dune::Multiscale::NewtonRightHandSide::DiscreteFunctionType &old_u_H, Dune::Multiscale::NewtonRightHandSide::DiscreteFunctionType &rhsVector) {
  rhsVector.clear();

  for (const auto& entity : rhsVector.space()) {
    const auto& geometry = entity.geometry();
    auto elementOfRHS = rhsVector.localFunction(entity);
    const auto baseSet = rhsVector.space().basisFunctionSet(entity);

    const auto old_u_H_loc = old_u_H.localFunction(entity);
    const auto quadrature = make_quadrature(entity, rhsVector.space(), quadratureOrder);

    const auto numDofs = elementOfRHS.numDofs();
    const auto numQuadraturePoints = quadrature.nop();
    // the return values:
    RangeType f_x;
    std::vector<RangeType> phi_x(numDofs);
    // gradient of base function and gradient of old_u_H
    std::vector<JacobianRangeType> grad_phi_x(numDofs);
    JacobianRangeType grad_old_u_H;
    // Let A denote the diffusion operator, then we save
    // A( \gradient old_u_H )
    JacobianRangeType diffusive_flux_in_grad_old_u_H;
    for (auto  quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint) {
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
      // get gradient of old u_H:
      old_u_H_loc.jacobian(quadrature[quadraturePoint], grad_old_u_H);
      // evaluate diffusion operator in x(=global_point) and grad_old_u_H
      A.diffusiveFlux(global_point, grad_old_u_H, diffusive_flux_in_grad_old_u_H);
      for (int i = 0; i < numDofs; ++i) {
        elementOfRHS[i] += det * quadrature.weight(quadraturePoint) * (f_x * phi_x[i]);
        elementOfRHS[i] -=
            det * quadrature.weight(quadraturePoint) * (diffusive_flux_in_grad_old_u_H[0] * grad_phi_x[i][0]);
      }
    }
  }
}



void Dune::Multiscale::NewtonRightHandSide::assemble_for_Newton_method(const Dune::Multiscale::CommonTraits::FirstSourceType &f, const Dune::Multiscale::CommonTraits::DiffusionType &A, const Dune::Multiscale::CommonTraits::LowerOrderTermType &F, const Dune::Multiscale::NewtonRightHandSide::DiscreteFunctionType &old_u_H, Dune::Multiscale::NewtonRightHandSide::DiscreteFunctionType &rhsVector) {
  rhsVector.clear();

  for (const auto& entity : rhsVector.space()) {
    const auto& geometry = entity.geometry();
    auto elementOfRHS = rhsVector.localFunction(entity);
    const auto baseSet = rhsVector.space().basisFunctionSet(entity);

    const auto old_u_H_loc = old_u_H.localFunction(entity);
    const auto quadrature = make_quadrature(entity, rhsVector.space(), quadratureOrder);

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


void Dune::Multiscale::NewtonRightHandSide::assemble_for_Newton_method(const Dune::Multiscale::CommonTraits::FirstSourceType &f, const Dune::Multiscale::CommonTraits::DiffusionType &A, const Dune::Multiscale::CommonTraits::LowerOrderTermType &F, const Dune::Multiscale::NewtonRightHandSide::DiscreteFunctionType &old_u_H, const Dune::Multiscale::NewtonRightHandSide::DiscreteFunctionType &dirichlet_extension, const Dune::Multiscale::CommonTraits::NeumannBCType &neumann_bc, Dune::Multiscale::NewtonRightHandSide::DiscreteFunctionType &rhsVector) {
  rhsVector.clear();

  for (const auto& entity : rhsVector.space()) {
    const auto& geometry = entity.geometry();
    auto elementOfRHS = rhsVector.localFunction(entity);
    const auto baseSet = rhsVector.space().basisFunctionSet(entity);

    const int numDofs = elementOfRHS.numDofs();

    std::vector<RangeType> phi_x(numDofs);
    // gradient of base function and gradient of old_u_H
    std::vector<JacobianRangeType> grad_phi_x(numDofs);

    const auto old_u_H_loc = old_u_H.localFunction(entity);
    const auto loc_dirichlet_extension = dirichlet_extension.localFunction(entity);
    const auto quadrature = make_quadrature(entity, rhsVector.space(), quadratureOrder);

    const auto& lagrangePointSet = rhsVector.space().lagrangePointSet(entity);

    for (const auto& intersection : DSC::intersectionRange(rhsVector.space().gridPart(), entity)) {
      if (!intersection.boundary())
        continue;
      // boundaryId 1 = Dirichlet face; boundaryId 2 = Neumann face;
      if (intersection.boundary() && (intersection.boundaryId() != 2))
        continue;

      const auto face = intersection.indexInInside();
      const auto faceQuadrature = make_quadrature(intersection, rhsVector.space(), quadratureOrder);
      for (auto faceQuadraturePoint : DSC::valueRange(faceQuadrature.nop())) {
        baseSet.evaluateAll(faceQuadrature[faceQuadraturePoint], phi_x);
        baseSet.jacobianAll(faceQuadrature[faceQuadraturePoint], grad_phi_x);

        const auto local_point_entity = faceQuadrature.point(faceQuadraturePoint);
        const auto global_point = geometry.global(local_point_entity);
        const auto local_point_face = intersection.geometry().local(global_point);

        RangeType neumann_value(0.0);
        neumann_bc.evaluate(global_point, neumann_value);

        const double face_weight =
            intersection.geometry().integrationElement(local_point_face) * faceQuadrature.weight(faceQuadraturePoint);

        for (const auto& lp : DSC::lagrangePointSetRange(lagrangePointSet, face)) {
          elementOfRHS[lp] += neumann_value * face_weight * phi_x[lp];
        }
      }
    }

    // the return values:
    RangeType f_x;

    JacobianRangeType gradient_dirichlet_extension;
    RangeType value_dirichlet_extension;
    JacobianRangeType grad_old_u_H;
    RangeType value_old_u_H;
    // Let A denote the diffusion operator, then we save
    // A( \gradient old_u_H )
    JacobianRangeType diffusive_flux;

    JacobianRangeType direction;

    for (auto quadraturePoint : DSC::valueRange(quadrature.nop())) {
      // local (barycentric) coordinates (with respect to entity)
      const auto& local_point = quadrature.point(quadraturePoint);
      const auto global_point = geometry.global(local_point);
      const double weight = geometry.integrationElement(local_point) * quadrature.weight(quadraturePoint);
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
      // get value of dirichlet extension:
      loc_dirichlet_extension.evaluate(quadrature[quadraturePoint], value_dirichlet_extension);
      // get gradient of dirichlet extension:
      loc_dirichlet_extension.jacobian(quadrature[quadraturePoint], gradient_dirichlet_extension);
      direction[0] = grad_old_u_H[0] + gradient_dirichlet_extension[0];
      // evaluate diffusion operator in x(=global_point) and grad_old_u_H
      A.diffusiveFlux(global_point, direction, diffusive_flux);

      RangeType F_x;
      F.evaluate(global_point, value_old_u_H + value_dirichlet_extension, direction, F_x);

      for (int i = 0; i < numDofs; ++i) {
        elementOfRHS[i] += weight * (f_x * phi_x[i]);
        elementOfRHS[i] -= weight * (diffusive_flux[0] * grad_phi_x[i][0]);
        elementOfRHS[i] -= weight * (F_x * phi_x[i]);
      }
    }
  }
}

void Dune::Multiscale::NewtonRightHandSide::assemble_for_HMM_Newton_method(const Dune::Multiscale::CommonTraits::FirstSourceType &f, const Dune::Multiscale::CommonTraits::DiffusionType &A, const Dune::Multiscale::NewtonRightHandSide::DiscreteFunctionType &old_u_H, const Dune::Multiscale::HMM::CellProblemNumberingManager &cp_num_manager, const Dune::Multiscale::HMM::HMMTraits::PeriodicDiscreteFunctionType &dummy_func, Dune::Multiscale::NewtonRightHandSide::DiscreteFunctionType &rhsVector) {
  const std::string cell_solution_location_baseSet = "/cell_problems/_cellSolutions_baseSet";
  const std::string cell_solution_location_discFunc = "/cell_problems/_cellSolutions_discFunc";

  // reader for the cell problem data file:
  DiscreteFunctionReader discrete_function_reader_baseSet(cell_solution_location_baseSet);
  // reader for the cell problem data file:
  DiscreteFunctionReader discrete_function_reader_discFunc(cell_solution_location_discFunc);

  const double delta = DSC_CONFIG_GET("hmm.delta", 1.0f);
  const double epsilon_estimated = DSC_CONFIG_GET("hmm.epsilon_guess", 1.0f);

  const DiscreteFunctionSpaceType& discreteFunctionSpace = rhsVector.space();
  const auto& periodicDiscreteFunctionSpace = dummy_func.space();

  // set rhsVector to zero:
  rhsVector.clear();

  int number_of_entity = 0;

  const auto macro_grid_endit = discreteFunctionSpace.end();
  for (auto macro_grid_it = discreteFunctionSpace.begin(); macro_grid_it != macro_grid_endit; ++macro_grid_it) {
    // it* Pointer auf ein Element der Entity
    const auto& macro_grid_geometry = (*macro_grid_it).geometry(); // Referenz auf Geometrie
    auto elementOfRHS = rhsVector.localFunction(*macro_grid_it);

    const auto macro_grid_baseSet = discreteFunctionSpace.basisFunctionSet(*macro_grid_it);
    const auto old_u_H_loc = old_u_H.localFunction(*macro_grid_it);
    // for \int_{\Omega} f \Phi
    const auto macro_quadrature = make_quadrature(*macro_grid_it, discreteFunctionSpace, quadratureOrder);
    // for - \int_{\Omega} \in_Y A^{\epsilon}( gradient reconstruction ) \nabla \Phi
    // the fine scale reconstructions are only available for the barycenter of the macro grid entity
    const auto macro_entity_barycenter = macro_grid_geometry.center();
    const auto barycenter_local = macro_grid_geometry.local(macro_entity_barycenter);
    const double macro_entity_volume = macro_grid_geometry.volume();
    const auto numDofs = elementOfRHS.numDofs();
    // gradient of base function and gradient of old_u_H
    std::vector<JacobianRangeType> grad_Phi_x_vec(numDofs);
    std::vector<RangeType> phi_x(numDofs);
    JacobianRangeType grad_old_u_H_x;

    // get gradient of old u_H:
    old_u_H_loc.jacobian(barycenter_local, grad_old_u_H_x);

    // Q_h(u_H^{(n-1}))(x_T,y):
    HMM::HMMTraits::PeriodicDiscreteFunctionType corrector_old_u_H("Corrector of u_H^(n-1)", periodicDiscreteFunctionSpace);
    corrector_old_u_H.clear();

    HMM::HMMTraits::PeriodicDiscreteFunctionType corrector_Phi_i("Corrector of Phi_i", periodicDiscreteFunctionSpace);
    discrete_function_reader_discFunc.read(number_of_entity, corrector_old_u_H);
    macro_grid_baseSet.jacobianAll(barycenter_local, grad_Phi_x_vec);

    for (int i = 0; i < numDofs; ++i) {
      const auto& grad_Phi_x = grad_Phi_x_vec[i];
      // --------------- the source contribution ( \int_{\Omega} f \Phi ) -------------------------------

      // the return values:
      RangeType f_x;

      const auto numMacroQuadraturePoints = macro_quadrature.nop();
      for (size_t quadraturePoint = 0; quadraturePoint < numMacroQuadraturePoints; ++quadraturePoint) {
        // local (barycentric) coordinates (with respect to entity)
        const auto& local_point = macro_quadrature.point(quadraturePoint);
        const auto global_point = macro_grid_geometry.global(local_point);
        const double quad_weight =
            macro_grid_geometry.integrationElement(local_point) * macro_quadrature.weight(quadraturePoint);
        // evaluate the Right Hand Side Function f at the current quadrature point and save its value in 'f_y':
        f.evaluate(global_point, f_x);
        //!TODO order of loops sucks
        macro_grid_baseSet.evaluateAll(macro_quadrature[quadraturePoint], phi_x);
        elementOfRHS[i] += quad_weight * (f_x * phi_x[i]);
      }

      // --------------- end of source contribution -----------------------------------------------------

      // --------------- the contribution of the jacobian of the diffusion operator, evaluated in the old
      // reconstructed macro solution -------------------------------

      corrector_Phi_i.clear();
      if (!DSC_CONFIG_GET("hmm.petrov_galerkin", true)) {
        typename EntityType::EntityPointer entity_ptr(macro_grid_it);
        discrete_function_reader_baseSet.read(cp_num_manager.get_number_of_cell_problem(entity_ptr, i),
                                              corrector_Phi_i);
      }

      RangeType fine_scale_contribution = 0.0;
      for (const auto& micro_grid_entity : periodicDiscreteFunctionSpace) {
        const auto& micro_grid_geometry = micro_grid_entity.geometry();
        assert(micro_grid_entity.partitionType() == InteriorEntity);

        auto loc_corrector_old_u_H = corrector_old_u_H.localFunction(micro_grid_entity);
        auto loc_corrector_Phi_i = corrector_Phi_i.localFunction(micro_grid_entity);

        // higher order quadrature, since A^{\epsilon} is highly variable
        const auto micro_grid_quadrature = make_quadrature(micro_grid_entity, periodicDiscreteFunctionSpace, quadratureOrder);
        const auto numQuadraturePoints = micro_grid_quadrature.nop();

        for (size_t microQuadraturePoint = 0; microQuadraturePoint < numQuadraturePoints; ++microQuadraturePoint) {
          // local (barycentric) coordinates (with respect to entity)
          const auto& local_micro_point = micro_grid_quadrature.point(microQuadraturePoint);

          const auto global_point_in_Y = micro_grid_geometry.global(local_micro_point);

          const double weight_micro_quadrature = micro_grid_quadrature.weight(microQuadraturePoint) *
                                                 micro_grid_geometry.integrationElement(local_micro_point);

          JacobianRangeType grad_corrector_old_u_H;
          loc_corrector_old_u_H.jacobian(micro_grid_quadrature[microQuadraturePoint], grad_corrector_old_u_H);

          JacobianRangeType grad_corrector_Phi_i;
          if (!DSC_CONFIG_GET("hmm.petrov_galerkin", true))
            loc_corrector_Phi_i.jacobian(micro_grid_quadrature[microQuadraturePoint], grad_corrector_Phi_i);

          // x_T + (delta * y)
          DomainType current_point_in_macro_grid;
          for (int k = 0; k < dimension; ++k)
            current_point_in_macro_grid[k] = macro_entity_barycenter[k] + (delta * global_point_in_Y[k]);

          // evaluate jacobian matrix of diffusion operator in 'position_vector' in direction 'direction_vector':

          JacobianRangeType direction_vector;
          for (int k = 0; k < dimension; ++k)
            direction_vector[0][k] = grad_old_u_H_x[0][k] + grad_corrector_old_u_H[0][k];

          JacobianRangeType diffusive_flux;
          A.diffusiveFlux(current_point_in_macro_grid, direction_vector, diffusive_flux);

          double cutting_function = 1.0;
          for (int k = 0; k < dimension; ++k) {
            // is the current quadrature point in the relevant cell?
            if (fabs(global_point_in_Y[k]) > (0.5 * (epsilon_estimated / delta))) {
              cutting_function *= 0.0;
            }
          }

          // if test function reconstruction = non-Petrov-Galerkin
          if (!DSC_CONFIG_GET("hmm.petrov_galerkin", true)) {
            JacobianRangeType grad_reconstruction_Phi_i;
            for (int k = 0; k < dimension; ++k)
              grad_reconstruction_Phi_i[0][k] = grad_Phi_x[0][k] + grad_corrector_Phi_i[0][k];

            fine_scale_contribution +=
                cutting_function * weight_micro_quadrature * (diffusive_flux[0] * grad_reconstruction_Phi_i[0]);
          } else {
            fine_scale_contribution +=
                cutting_function * weight_micro_quadrature * (diffusive_flux[0] * grad_Phi_x[0]);
          }
        }
      }

      elementOfRHS[i] -= pow(delta / epsilon_estimated, dimension) * macro_entity_volume * fine_scale_contribution;

      // --------------- end of diffusion contribution -----------------------------------------------------
    }

    number_of_entity += 1;
  }
}
