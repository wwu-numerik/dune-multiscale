#ifndef DUNE_MULTISCALE_COMMON_NEWTON_RHS_HH
#define DUNE_MULTISCALE_COMMON_NEWTON_RHS_HH

#include <config.h>

#include <dune/multiscale/common/traits.hh>

#include <dune/fem/quadrature/quadrature.hh>
#include <dune/multiscale/hmm/cell_problem_numbering.hh>
#include <dune/multiscale/hmm/cell_problem_solver.hh>

#include <dune/multiscale/tools/misc.hh>
#include <dune/multiscale/common/traits.hh>
#include <dune/multiscale/msfem/localproblems/subgrid-list.hh>
#include <dune/multiscale/problems/selector.hh>
#include <dune/stuff/fem/functions/checks.hh>
#include <dune/multiscale/msfem/msfem_traits.hh>
#include <dune/multiscale/msfem/localproblems/localsolutionmanager.hh>
#include <dune/multiscale/common/dirichletconstraints.hh>


#include <dune/stuff/common/logging.hh>
#include <dune/stuff/functions/interfaces.hh>

namespace Dune {
namespace Multiscale {

//! Assembler for right rand side
//! We assemble the right hand side in a LSE, i.e. f \cdot \Phi_H + G \cdot \nabala \Phi_H
//! we call f the first Source and G the second Source

class NewtonRightHandSide {
  typedef CommonTraits::DiscreteFunctionType DiscreteFunctionType;
  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

  typedef typename DiscreteFunctionSpaceType::BasisFunctionSetType BasisFunctionSetType;
  typedef typename DiscreteFunctionSpaceType::RangeType RangeType;
  typedef typename DiscreteFunctionSpaceType::DomainFieldType DomainFieldType;
  typedef DomainFieldType TimeType;
  typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
  typedef typename GridPartType::GridType GridType;
  typedef typename DiscreteFunctionSpaceType::JacobianRangeType JacobianRangeType;
  typedef typename DiscreteFunctionSpaceType::DomainType DomainType;
  typedef typename GridType::Codim<0>::Entity EntityType;
  typedef typename EntityType::Geometry GeometryType;

  typedef MsFEM::LocalSolutionManager LocalSolutionManagerType;

  static const int dimension = GridType::dimension;
  static const int polynomialOrder = DiscreteFunctionSpaceType::polynomialOrder;
  static const int quadratureOrder = 2*polynomialOrder + 2;
public:
/**
 * The rhs-assemble()-methods for non-linear elliptic problems:
 * discreteFunction is an output parameter (kind of return value)
 **/
static void assemble_for_Newton_method(const CommonTraits::FirstSourceType& f, const CommonTraits::DiffusionType& A,
                                       const DiscreteFunctionType& old_u_H, // old_u_H from the last iteration step
                                       DiscreteFunctionType& rhsVector); // end method


/**
 * The rhs-assemble()-methods for non-linear elliptic problems
 * if there is a first source f and a lower order term F:
 * discreteFunction is an output parameter (kind of return value)
 **/
static void assemble_for_Newton_method(const CommonTraits::FirstSourceType& f, const CommonTraits::DiffusionType& A,
                                       const CommonTraits::LowerOrderTermType& F,
                                       const DiscreteFunctionType& old_u_H, // old_u_H from the last iteration step
                                       DiscreteFunctionType& rhsVector); // end method

/**
 * The rhs-assemble()-methods for non-linear elliptic problems
 * if there is a first source f, a lower order term F
 * and Dirichlet and Neumann boundary conditions
 * discreteFunction is an output parameter (kind of return value)
 * \param old_u_H from the last iteration step
 * \param dirichlet_extension discrete function describing dirichlet extension
 **/
static void assemble_for_Newton_method(
    const CommonTraits::FirstSourceType& f, const CommonTraits::DiffusionType& A, const CommonTraits::LowerOrderTermType& F,
    const DiscreteFunctionType& old_u_H,
    const DiscreteFunctionType& dirichlet_extension,
    const CommonTraits::NeumannBCType& neumann_bc, DiscreteFunctionType& rhsVector); // end method

//! The rhs-assemble()-methods for non-linear elliptic problems, solved with the heterogenous multiscale method
// ( requires reconstruction of old_u_H and local fine scale averages )

static void
assemble_for_HMM_Newton_method(const CommonTraits::FirstSourceType& f, const CommonTraits::DiffusionType& A,
                               const DiscreteFunctionType& old_u_H, // old_u_H from the last iteration step
                               // to obtain some information about the periodic discrete function space (space
                               // for the cell problems)
                               const HMM::CellProblemNumberingManager& cp_num_manager,
                               const HMM::HMMTraits::PeriodicDiscreteFunctionType& dummy_func, DiscreteFunctionType& rhsVector) {
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
} // end method
};

} // end namespace Multiscale
} // end namespace Dune
#endif // DUNE_MULTISCALE_COMMON_NEWTON_RHS_HH
