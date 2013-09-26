#include "elliptic_hmm_matrix_assembler.hh"

#include <dune/multiscale/hmm/cell_problem_numbering.hh>
#include <dune/multiscale/tools/discretefunctionwriter.hh>
#include <dune/stuff/fem/localmatrix_proxy.hh>
#include <dune/stuff/common/memory.hh>
#include <dune/multiscale/problems/base.hh>

namespace Dune {
namespace Multiscale {
namespace HMM {

void DiscreteEllipticHMMOperator::boundary_treatment(CommonTraits::FEMMatrix& global_matrix) const {
  const auto& gridPart = discreteFunctionSpace_.gridPart();
  for (const auto& entity : discreteFunctionSpace_) {
    if (!entity.hasBoundaryIntersections())
      continue;

    auto local_matrix = global_matrix.localMatrix(entity, entity);

    const auto& lagrangePointSet = discreteFunctionSpace_.lagrangePointSet(entity);

    const IntersectionIterator iend = gridPart.iend(entity);
    for (IntersectionIterator iit = gridPart.ibegin(entity); iit != iend; ++iit) {
      const Intersection& intersection = *iit;
      if (!intersection.boundary())
        continue;

      // boundaryId 1 = Dirichlet face; boundaryId 2 = Neumann face;
      if (intersection.boundary() && (intersection.boundaryId() == 2))
        continue;

      const int face = intersection.indexInInside();
      const auto fdend = lagrangePointSet.endSubEntity<1>(face);
      for (auto fdit = lagrangePointSet.beginSubEntity<1>(face); fdit != fdend; ++fdit)
        local_matrix.unitRow(*fdit);
    }
  }
}

void DiscreteEllipticHMMOperator::assemble_matrix(CommonTraits::FEMMatrix& global_matrix) const {
  // if test function reconstruction
  if (!DSC_CONFIG_GET("hmm.petrov_galerkin", true))
    DSC_LOG_INFO << "Assembling classical (non-Petrov-Galerkin) HMM Matrix." << std::endl;
  else
    DSC_LOG_INFO << "Assembling (Petrov-Galerkin) HMM Matrix." << std::endl;

  // place, where we saved the solutions of the cell problems
  const std::string cell_solution_location = filename_ + "/cell_problems/_cellSolutions_baseSet";
  const double delta = DSC_CONFIG_GET("hmm.delta", 1.0f);
  const double epsilon_estimated = DSC_CONFIG_GET("hmm.epsilon_guess", 1.0f);

  // reader for the cell problem data file:
  DiscreteFunctionReader discrete_function_reader(cell_solution_location);

  global_matrix.reserve();
  global_matrix.clear();

  std::vector<typename BaseFunctionSet::JacobianRangeType> gradient_Phi(discreteFunctionSpace_.mapper().maxNumDofs());

  for (const auto& macro_grid_entity : discreteFunctionSpace_) {
    const auto& macro_grid_geometry = macro_grid_entity.geometry();
    assert(macro_grid_entity.partitionType() == InteriorEntity);

    DSFe::LocalMatrixProxy<CommonTraits::FEMMatrix> local_matrix(global_matrix, macro_grid_entity, macro_grid_entity);

    const auto& macro_grid_baseSet = local_matrix.domainBasisFunctionSet();
    const auto numMacroBaseFunctions = macro_grid_baseSet.size();
    const auto macro_entity_barycenter = macro_grid_geometry.center();
    const auto center_local = macro_grid_geometry.local(macro_entity_barycenter);
    const double macro_entity_volume = macro_grid_geometry.volume();

    std::vector<int> cell_problem_id(numMacroBaseFunctions, 0);

    typedef std::unique_ptr<PeriodicDiscreteFunction> PeriodicDiscreteFunctionPointer;
    std::vector<PeriodicDiscreteFunctionPointer> corrector_Phi(discreteFunctionSpace_.mapper().maxNumDofs());

    macro_grid_baseSet.jacobianAll(center_local, gradient_Phi);

    //!TODO generator functions
    for (unsigned int i = 0; i < numMacroBaseFunctions; ++i) {
      // get number of cell problem from entity and number of base function
      typename Entity::EntityPointer macro_entity_pointer(macro_grid_entity);
      cell_problem_id[i] = cp_num_manager_.get_number_of_cell_problem(macro_entity_pointer, i);

      corrector_Phi[i] =
          DSC::make_unique<PeriodicDiscreteFunction>("Corrector Function of Phi", periodicDiscreteFunctionSpace_);
      corrector_Phi[i]->clear();
      discrete_function_reader.read(cell_problem_id[i], *(corrector_Phi[i]));
    }

    for (unsigned int i = 0; i < numMacroBaseFunctions; ++i) {
      for (unsigned int j = 0; j < numMacroBaseFunctions; ++j) {
        RangeType fine_scale_average = 0.0;

        // nur checken ob der momentane Quadraturpunkt in der Zelle delta/epsilon_estimated*Y ist (0 bzw. 1 =>
        // Abschneidefunktion!)

        for (const auto& micro_grid_entity : periodicDiscreteFunctionSpace_) {
          const auto& micro_grid_geometry = micro_grid_entity.geometry();
          assert(micro_grid_entity.partitionType() == InteriorEntity);

          const auto localized_corrector_i = corrector_Phi[i]->localFunction(micro_grid_entity);
          const auto localized_corrector_j = corrector_Phi[j]->localFunction(micro_grid_entity);

          // higher order quadrature, since A^{\epsilon} is highly variable
          const auto micro_grid_quadrature = make_quadrature(micro_grid_entity, periodicDiscreteFunctionSpace_);
          const auto numQuadraturePoints = micro_grid_quadrature.nop();

          for (size_t microQuadraturePoint = 0; microQuadraturePoint < numQuadraturePoints; ++microQuadraturePoint) {
            // local (barycentric) coordinates (with respect to entity)
            const auto& local_micro_point = micro_grid_quadrature.point(microQuadraturePoint);
            const auto global_point_in_Y = micro_grid_geometry.global(local_micro_point);
            const double weight_micro_quadrature = micro_grid_quadrature.weight(microQuadraturePoint) *
                                                   micro_grid_geometry.integrationElement(local_micro_point);

            JacobianRangeType grad_corrector_i, grad_corrector_j;
            localized_corrector_i.jacobian(micro_grid_quadrature[microQuadraturePoint], grad_corrector_i);
            localized_corrector_j.jacobian(micro_grid_quadrature[microQuadraturePoint], grad_corrector_j);

            // x_T + (delta * y)
            DomainType current_point_in_macro_grid;
            for (int k = 0; k < dimension; ++k)
              current_point_in_macro_grid[k] = macro_entity_barycenter[k] + (delta * global_point_in_Y[k]);

            JacobianRangeType direction_of_diffusion;
            for (int k = 0; k < dimension; ++k)
              direction_of_diffusion[0][k] = gradient_Phi[i][0][k] + grad_corrector_i[0][k];

            JacobianRangeType diffusion_in_gradient_Phi_reconstructed;
            diffusion_operator_.diffusiveFlux(current_point_in_macro_grid, direction_of_diffusion,
                                              diffusion_in_gradient_Phi_reconstructed);

            double cutting_function = 1.0;
            for (int k = 0; k < dimension; ++k) {
              // is the current quadrature point in the relevant cell?
              if (fabs(global_point_in_Y[k]) > (0.5 * (epsilon_estimated / delta))) {
                cutting_function *= 0.0;
              }
            }

            // if test function reconstruction = non-Petrov-Galerkin HMM
            if (!DSC_CONFIG_GET("hmm.petrov_galerkin", true)) {
              JacobianRangeType grad_reconstruction_Phi_j;
              for (int k = 0; k < dimension; ++k)
                grad_reconstruction_Phi_j[0][k] = gradient_Phi[j][0][k] + grad_corrector_j[0][k];

              fine_scale_average += cutting_function * weight_micro_quadrature *
                                    (diffusion_in_gradient_Phi_reconstructed[0] * grad_reconstruction_Phi_j[0]);
            } else { // if Petrov-Galerkin-HMM
              fine_scale_average += cutting_function * weight_micro_quadrature *
                                    (diffusion_in_gradient_Phi_reconstructed[0] * gradient_Phi[j][0]);
            }
          }
        }

        // add |T| * (delta/epsilon)^N \int_Y ...
        local_matrix.add(j, i, pow(delta / epsilon_estimated, dimension) * macro_entity_volume * fine_scale_average);
      }
    }
  }
  boundary_treatment(global_matrix);
} // assemble_matrix

//! assemble stiffness matrix for HMM with Newton Method
void DiscreteEllipticHMMOperator::assemble_jacobian_matrix(DiscreteFunction& old_u_H /*u_H^(n-1)*/,
                                                           CommonTraits::FEMMatrix& global_matrix) const {

  // if test function reconstruction
  if (!DSC_CONFIG_GET("hmm.petrov_galerkin", true))
    DSC_LOG_INFO << "Assembling classical (non-Petrov-Galerkin) HMM Matrix for Newton Iteration." << std::endl;
  else
    DSC_LOG_INFO << "Assembling (Petrov-Galerkin) HMM Matrix for Newton Iteration." << std::endl;

  // place, where we saved the solutions of the cell problems
  const std::string cell_solution_location_baseSet = filename_ + "/cell_problems/_cellSolutions_baseSet";
  const std::string cell_solution_location_discFunc = filename_ + "/cell_problems/_cellSolutions_discFunc";
  const std::string jac_cor_cell_solution_location_baseSet_discFunc =
      filename_ + "/cell_problems/_JacCorCellSolutions_baseSet_discFunc";
  const double delta = DSC_CONFIG_GET("hmm.delta", 1.0f);
  const double epsilon_estimated = DSC_CONFIG_GET("hmm.epsilon_guess", 1.0f);

  // reader for the cell problem data file:
  DiscreteFunctionReader discrete_function_reader_baseSet(cell_solution_location_baseSet);
  DiscreteFunctionReader discrete_function_reader_discFunc(cell_solution_location_discFunc);
  DiscreteFunctionReader discrete_function_reader_jac_cor(jac_cor_cell_solution_location_baseSet_discFunc);

  global_matrix.reserve();
  global_matrix.clear();

  std::vector<typename BaseFunctionSet::JacobianRangeType> gradient_Phi(discreteFunctionSpace_.mapper().maxNumDofs());
  std::vector<typename BaseFunctionSet::JacobianRangeType> gradient_Phi_new(
      discreteFunctionSpace_.mapper().maxNumDofs());

  int number_of_macro_entity = 0;

  for (const auto& macro_grid_entity : discreteFunctionSpace_) {
    const auto& macro_grid_geometry = macro_grid_entity.geometry();
    assert(macro_grid_entity.partitionType() == InteriorEntity);

    DSFe::LocalMatrixProxy<CommonTraits::FEMMatrix> local_matrix(global_matrix, macro_grid_entity, macro_grid_entity);
    auto local_old_u_H = old_u_H.localFunction(macro_grid_entity);

    const auto& macro_grid_baseSet = local_matrix.domainBasisFunctionSet();
    const auto numMacroBaseFunctions = macro_grid_baseSet.size();
    const auto macro_entity_barycenter = macro_grid_geometry.center();
    const auto center_local = macro_grid_geometry.local(macro_entity_barycenter);
    const double macro_entity_volume = macro_grid_geometry.volume();

    std::vector<int> cell_problem_id(numMacroBaseFunctions, -1);

    // \nabla_x u_H^{(n-1})(x_T)
    typename BaseFunctionSet::JacobianRangeType grad_old_u_H;
    local_old_u_H.jacobian(center_local, grad_old_u_H);
    // here: no multiplication with jacobian inverse transposed required!

    // Q_h(u_H^{(n-1}))(x_T,y):
    PeriodicDiscreteFunction corrector_old_u_H("Corrector of u_H^(n-1)", periodicDiscreteFunctionSpace_);
    corrector_old_u_H.clear();

    discrete_function_reader_discFunc.read(number_of_macro_entity, corrector_old_u_H);

    std::vector<std::unique_ptr<PeriodicDiscreteFunction>> corrector_Phi(discreteFunctionSpace_.mapper().maxNumDofs());
    macro_grid_baseSet.jacobianAll(center_local, gradient_Phi);

    // gradients of macrocopic base functions:
    // TODO generator
    for (unsigned int i = 0; i < numMacroBaseFunctions; ++i) {
      // get number of cell problem from entity and number of base function
      typename Entity::EntityPointer macro_entity_pointer(macro_grid_entity);
      cell_problem_id[i] = cp_num_manager_.get_number_of_cell_problem(macro_entity_pointer, i);

      if (!DSC_CONFIG_GET("hmm.petrov_galerkin", true)) {
        corrector_Phi[i] = std::unique_ptr<PeriodicDiscreteFunction>(
            new PeriodicDiscreteFunction("Corrector Function of Phi_j", periodicDiscreteFunctionSpace_));
        corrector_Phi[i]->clear();
        discrete_function_reader_baseSet.read(cell_problem_id[i], *(corrector_Phi[i]));
      }
    }
    // the multiplication with jacobian inverse is delegated
    macro_grid_baseSet.jacobianAll(center_local, gradient_Phi_new);
    assert(gradient_Phi == gradient_Phi_new);

    for (unsigned int i = 0; i < numMacroBaseFunctions; ++i) {
      // D_Q(\Phi_i,u_H^{n-1})
      // the jacobian of the corrector operator applied to u_H^{(n-1)} in direction of gradient \Phi_i
      PeriodicDiscreteFunction jacobian_corrector_old_u_H_Phi_i("Jacobian Corrector Function of u_H^(n-1) and Phi_i",
                                                                periodicDiscreteFunctionSpace_);
      jacobian_corrector_old_u_H_Phi_i.clear();

      discrete_function_reader_jac_cor.read(cell_problem_id[i], jacobian_corrector_old_u_H_Phi_i);

      for (unsigned int j = 0; j < numMacroBaseFunctions; ++j) {
        RangeType fine_scale_average = 0.0;

        for (const Entity& micro_grid_entity : periodicDiscreteFunctionSpace_) {
          const auto& micro_grid_geometry = micro_grid_entity.geometry();
          assert(micro_grid_entity.partitionType() == InteriorEntity);
          const auto loc_corrector_old_u_H = corrector_old_u_H.localFunction(micro_grid_entity);
          const auto loc_D_Q_old_u_H_Phi_i = jacobian_corrector_old_u_H_Phi_i.localFunction(micro_grid_entity);

          // higher order quadrature, since A^{\epsilon} is highly variable
          const auto micro_grid_quadrature = make_quadrature(micro_grid_entity, periodicDiscreteFunctionSpace_);
          const auto numQuadraturePoints = micro_grid_quadrature.nop();

          for (size_t microQuadraturePoint = 0; microQuadraturePoint < numQuadraturePoints; ++microQuadraturePoint) {
            // local (barycentric) coordinates (with respect to entity)
            const auto& local_micro_point = micro_grid_quadrature.point(microQuadraturePoint);
            auto global_point_in_Y = micro_grid_geometry.global(local_micro_point);
            const double weight_micro_quadrature = micro_grid_quadrature.weight(microQuadraturePoint) *
                                                   micro_grid_geometry.integrationElement(local_micro_point);

            JacobianRangeType grad_corrector_old_u_H, grad_D_Q_old_u_H_Phi_i;
            loc_corrector_old_u_H.jacobian(micro_grid_quadrature[microQuadraturePoint], grad_corrector_old_u_H);
            loc_D_Q_old_u_H_Phi_i.jacobian(micro_grid_quadrature[microQuadraturePoint], grad_D_Q_old_u_H_Phi_i);

            JacobianRangeType grad_corrector_Phi_j;
            if (!DSC_CONFIG_GET("hmm.petrov_galerkin", true)) {
              auto loc_corrector_Phi_j = corrector_Phi[j]->localFunction(micro_grid_entity);
              loc_corrector_Phi_j.jacobian(micro_grid_quadrature[microQuadraturePoint], grad_corrector_Phi_j);
            }

            // x_T + (delta * y)
            DomainType current_point_in_macro_grid;
            for (int k = 0; k < dimension; ++k)
              current_point_in_macro_grid[k] = macro_entity_barycenter[k] + (delta * global_point_in_Y[k]);

            // evaluate jacobian matrix of diffusion operator in 'position_vector' in direction 'direction_vector':

            JacobianRangeType position_vector;
            for (int k = 0; k < dimension; ++k)
              position_vector[0][k] = grad_old_u_H[0][k] + grad_corrector_old_u_H[0][k];

            JacobianRangeType direction_vector;
            for (int k = 0; k < dimension; ++k)
              direction_vector[0][k] = gradient_Phi[i][0][k] + grad_D_Q_old_u_H_Phi_i[0][k];

            JacobianRangeType jac_diffusion_flux;
            diffusion_operator_.jacobianDiffusiveFlux(current_point_in_macro_grid, position_vector, direction_vector,
                                                      jac_diffusion_flux);

            double cutting_function = 1.0;
            for (int k = 0; k < dimension; ++k) {
              // is the current quadrature point in the relevant cell?
              if (fabs(global_point_in_Y[k]) > (0.5 * (epsilon_estimated / delta))) {
                cutting_function *= 0.0;
              }
            }

            // if test function reconstruction
            if (!DSC_CONFIG_GET("hmm.petrov_galerkin", true)) {
              JacobianRangeType grad_reconstruction_Phi_j;
              for (int k = 0; k < dimension; ++k)
                grad_reconstruction_Phi_j[0][k] = gradient_Phi[j][0][k] + grad_corrector_Phi_j[0][k];

              fine_scale_average +=
                  cutting_function * weight_micro_quadrature * (jac_diffusion_flux[0] * grad_reconstruction_Phi_j[0]);
            } else {
              fine_scale_average +=
                  cutting_function * weight_micro_quadrature * (jac_diffusion_flux[0] * gradient_Phi[j][0]);
            }
          }
        }

        // add |T| * (delta/epsilon)^N \int_Y ...
        local_matrix.add(j, i, pow(delta / epsilon_estimated, dimension) * macro_entity_volume * fine_scale_average);
      }
    }
    number_of_macro_entity += 1;
  }

  boundary_treatment(global_matrix);
} // assemble_jacobian_matrix

//! dummy implementation of "operator()"
//! 'w' = effect of the discrete operator on 'u'
void DiscreteEllipticHMMOperator::operator()(const DiscreteFunction& /*u*/, DiscreteFunction& /*w*/) const {
  DUNE_THROW(Dune::NotImplemented,
             "the ()-operator of the DiscreteEllipticHMMOperator class is not yet implemented and still a dummy.");
}

} // namespace HMM {
} // namespace Multiscale {
} // namespace Dune {
