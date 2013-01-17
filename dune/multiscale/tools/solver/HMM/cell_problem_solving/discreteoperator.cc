#include <dune/multiscale/problems/elliptic_problems/selector.hh>
#include <dune/stuff/common/parameter/configcontainer.hh>

namespace Dune {

template< class DiscreteFunctionImp, class DiffusionImp >
void DiscreteCellProblemOperator< DiscreteFunctionImp, DiffusionImp >::operator()(const DiscreteFunction&/* u*/,
                                                                                  DiscreteFunction& /*w*/) const {
  std::cout << "the ()-operator of the DiscreteCellProblemOperator class is not yet implemented and still a dummy."
            << std::endl;
  std::abort();
}

template< class PeriodicDiscreteFunctionImp, class DiffusionImp >
template< class MatrixType >
void DiscreteCellProblemOperator< PeriodicDiscreteFunctionImp, DiffusionImp >::assemble_matrix(
  const DomainType& x_T,
  MatrixType&
  global_matrix) const
{
  // x_T is the barycenter of the macro grid element T
  typedef typename MatrixType::LocalMatrixType LocalMatrix;

  const double delta = DSC_CONFIG_GET("problem.delta", 1.0f);

  global_matrix.reserve();
  global_matrix.clear();

  // micro scale base function:
  std::vector< RangeType > phi( periodicDiscreteFunctionSpace_.mapper().maxNumDofs() );

  // gradient of micro scale base function:
  std::vector< typename BaseFunctionSet::JacobianRangeType > gradient_phi(
    periodicDiscreteFunctionSpace_.mapper().maxNumDofs() );

  const Iterator end = periodicDiscreteFunctionSpace_.end();
  for (Iterator it = periodicDiscreteFunctionSpace_.begin(); it != end; ++it)
  {
    const Entity& cell_grid_entity = *it;
    const Geometry& cell_grid_geometry = cell_grid_entity.geometry();
    assert(cell_grid_entity.partitionType() == InteriorEntity);

    LocalMatrix local_matrix = global_matrix.localMatrix(cell_grid_entity, cell_grid_entity);

    const BaseFunctionSet& baseSet = local_matrix.domainBaseFunctionSet();
    const unsigned int numBaseFunctions = baseSet.size();

    // for constant diffusion "2*discreteFunctionSpace_.order()" is sufficient, for the general case, it is better to
    // use a higher order quadrature:
    const Quadrature quadrature(cell_grid_entity, 2 * periodicDiscreteFunctionSpace_.order() + 2);
    const size_t numQuadraturePoints = quadrature.nop();
    for (size_t quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint)
    {
      // local (barycentric) coordinates (with respect to cell grid entity)
      const typename Quadrature::CoordinateType& local_point = quadrature.point(quadraturePoint);

      // global point in the unit cell Y
      DomainType global_point = cell_grid_geometry.global(local_point);

      // x_T + (delta * global_point)
      DomainType x_T_delta_global_point;
      for (int k = 0; k < dimension; ++k)
      {
        x_T_delta_global_point[k] = x_T[k] + (delta * global_point[k]);
      }

      const double weight = quadrature.weight(quadraturePoint)
                            * cell_grid_geometry.integrationElement(local_point);

      // transposed of the the inverse jacobian
      const auto& inverse_jac = cell_grid_geometry.jacobianInverseTransposed(local_point);
      baseSet.jacobianAll(quadrature[quadraturePoint], inverse_jac, gradient_phi);
      baseSet.evaluateAll(quadrature[quadraturePoint], phi);

      for (unsigned int i = 0; i < numBaseFunctions; ++i)
      {
        // A( x_T + \delta y, \nabla \phi )
        // diffusion operator evaluated in (x_T + \delta y , \nabla \phi)
        typename LocalFunction::JacobianRangeType diffusion_in_gradient_phi;
        diffusion_operator_.diffusiveFlux(x_T_delta_global_point, gradient_phi[i], diffusion_in_gradient_phi);
        for (unsigned int j = 0; j < numBaseFunctions; ++j)
        {
          // stiffness contribution
          local_matrix.add( j, i, weight * (diffusion_in_gradient_phi[0] * gradient_phi[j][0]) );

          // mass contribution
          local_matrix.add( j, i, CELL_MASS_WEIGHT * weight * (phi[i][0] * phi[j][0]) );
        }
      }
    }
  }
} // assemble_matrix

template< class PeriodicDiscreteFunctionImp, class DiffusionImp >
template< class MatrixType >
void DiscreteCellProblemOperator< PeriodicDiscreteFunctionImp, DiffusionImp >::assemble_jacobian_matrix(
  const DomainType& x_T,
  const JacobianRangeType& grad_coarse_function,
  const PeriodicDiscreteFunctionImp& old_fine_function,
  MatrixType& global_matrix) const
{
  const double delta = DSC_CONFIG_GET("problem.delta", 1.0f);

  typedef typename MatrixType::LocalMatrixType LocalMatrix;
  typedef typename PeriodicDiscreteFunctionImp::LocalFunctionType
  LocalFunction;

  global_matrix.reserve();
  global_matrix.clear();

  // micro scale base function:
  std::vector< RangeType > phi( periodicDiscreteFunctionSpace_.mapper().maxNumDofs() );

  // gradient of micro scale base function:
  std::vector< typename BaseFunctionSet::JacobianRangeType > gradient_phi(
    periodicDiscreteFunctionSpace_.mapper().maxNumDofs() );

  const Iterator end = periodicDiscreteFunctionSpace_.end();
  for (Iterator it = periodicDiscreteFunctionSpace_.begin(); it != end; ++it)
  {
    const Entity& cell_grid_entity = *it;
    const Geometry& cell_grid_geometry = cell_grid_entity.geometry();
    assert(cell_grid_entity.partitionType() == InteriorEntity);

    LocalMatrix local_matrix = global_matrix.localMatrix(cell_grid_entity, cell_grid_entity);
    LocalFunction local_fine_function = old_fine_function.localFunction(cell_grid_entity);

    const BaseFunctionSet& baseSet = local_matrix.domainBaseFunctionSet();
    const unsigned int numBaseFunctions = baseSet.size();

    // for constant diffusion "2*periodicDiscreteFunctionSpace_.order()" is sufficient, for the general case, it is
    // better to use a higher order quadrature:
    const Quadrature quadrature(cell_grid_entity, 2 * periodicDiscreteFunctionSpace_.order() + 2);
    const size_t numQuadraturePoints = quadrature.nop();
    for (size_t quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint)
    {
      // local (barycentric) coordinates (with respect to entity)
      const typename Quadrature::CoordinateType& local_point = quadrature.point(quadraturePoint);

      // global point in the unit cell Y
      const DomainType global_point = cell_grid_geometry.global(local_point);

      // x_T + (delta * global_point)
      DomainType x_T_delta_global_point;
      for (int k = 0; k < dimension; ++k)
      {
        x_T_delta_global_point[k] = x_T[k] + (delta * global_point[k]);
      }

      const double weight = quadrature.weight(quadraturePoint) * cell_grid_geometry.integrationElement(local_point);

      // transposed of the the inverse jacobian
      const FieldMatrix< double, dimension, dimension >& inverse_jac
        = cell_grid_geometry.jacobianInverseTransposed(local_point);

      for (unsigned int i = 0; i < numBaseFunctions; ++i)
      {
        // jacobian of the base functions, with respect to the reference element
        typename BaseFunctionSet::JacobianRangeType gradient_phi_ref_element;
        baseSet.jacobian(i, quadrature[quadraturePoint], gradient_phi_ref_element);

        // multiply it with transpose of jacobian inverse to obtain the jacobian with respect to the real entity
        inverse_jac.mv(gradient_phi_ref_element[0], gradient_phi[i][0]);

        baseSet.evaluate(i, quadrature[quadraturePoint], phi[i]);
      }

      for (unsigned int i = 0; i < numBaseFunctions; ++i)
      {
        // JA( x_T + \delta y, \nabla_x PHI_H(x_T) + \nabla_y old_fine_function ) \nabla \phi
        // jacobian matrix of the diffusion operator evaluated in (x_T + \delta y , \nabla_x PHI_H(x_T) + \nabla_y
        // old_fine_function ) in direction \nabla \phi
        typename BaseFunctionSet::JacobianRangeType grad_local_fine_function;
        local_fine_function.jacobian(quadrature[quadraturePoint], grad_local_fine_function);
        // here: no multiplication with jacobian inverse transposed required!

        typename BaseFunctionSet::JacobianRangeType position_vector;
        for (int k = 0; k < dimension; ++k)
        {
          position_vector[0][k] = grad_coarse_function[0][k] + grad_local_fine_function[0][k];
        }

        // jacobian of diffusion operator evaluated in (x,grad coarse + grad fine) in direction of the gradient of the
        // current base function
        typename LocalFunction::JacobianRangeType jac_diffusion_flux;
        diffusion_operator_.jacobianDiffusiveFlux(x_T_delta_global_point,
                                                  position_vector,
                                                  gradient_phi[i], jac_diffusion_flux);

        for (unsigned int j = 0; j < numBaseFunctions; ++j)
        {
          // stiffness contribution
          local_matrix.add( j, i, weight * (jac_diffusion_flux[0] * gradient_phi[j][0]) );

          // mass contribution
          local_matrix.add( j, i, CELL_MASS_WEIGHT * weight * (phi[i][0] * phi[j][0]) );
        }
      }
    }
  }
} // assemble_jacobian_matrix

template< class DiscreteFunctionImp, class DiffusionImp >
void DiscreteCellProblemOperator< DiscreteFunctionImp, DiffusionImp >::printCellRHS(const DiscreteFunctionImp& rhs) const {
  typedef typename DiscreteFunctionImp::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionSpaceType::IteratorType        IteratorType;
  typedef typename DiscreteFunctionImp::LocalFunctionType         LocalFunctionType;

  const DiscreteFunctionSpaceType& discreteFunctionSpace
    = rhs.space();

  const IteratorType endit = discreteFunctionSpace.end();
  for (IteratorType it = discreteFunctionSpace.begin(); it != endit; ++it)
  {
    const LocalFunctionType elementOfRHS = rhs.localFunction(*it);

    const int numDofs = elementOfRHS.numDofs();
    for (int i = 0; i < numDofs; ++i)
    {
      std::cout << "Number of Dof: " << i << " ; " << rhs.name() << " : " << elementOfRHS[i] << std::endl;
    }
  }
}  // end method

template< class DiscreteFunctionImp, class DiffusionImp >
double DiscreteCellProblemOperator< DiscreteFunctionImp, DiffusionImp >::normRHS(const DiscreteFunctionImp& rhs) const {
  double norm = 0.0;

  typedef typename DiscreteFunctionImp::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionSpaceType::IteratorType        IteratorType;
  typedef typename IteratorType::Entity                           EntityType;
  typedef typename DiscreteFunctionImp::LocalFunctionType         LocalFunctionType;
  typedef typename DiscreteFunctionSpaceType::GridPartType        GridPartType;
  typedef typename DiscreteFunctionSpaceType::GridType            GridType;
  typedef typename GridType::template Codim< 0 >::Geometry
  EnGeometryType;

  const DiscreteFunctionSpaceType& discreteFunctionSpace
    = rhs.space();

  const IteratorType endit = discreteFunctionSpace.end();
  for (IteratorType it = discreteFunctionSpace.begin(); it != endit; ++it)
  {
    // entity
    const EntityType& entity = *it;

    // create quadrature for given geometry type
    const CachingQuadrature< GridPartType, 0 > quadrature(entity, 2 * discreteFunctionSpace.order() + 2);

    // get geoemetry of entity
    const EnGeometryType& geo = entity.geometry();

    const LocalFunctionType localRHS = rhs.localFunction(*it);

    // integrate
    const int quadratureNop = quadrature.nop();
    for (int quadraturePoint = 0; quadraturePoint < quadratureNop; ++quadraturePoint)
    {
      const double weight = quadrature.weight(quadraturePoint)
                            * geo.integrationElement( quadrature.point(quadraturePoint) );

      RangeType value(0.0);
      localRHS.evaluate(quadrature[quadraturePoint], value);

      norm += weight * value * value;
    }
  }
  return norm;
}  // end method

template< class DiscreteFunctionImp, class DiffusionImp >
// template< class MatrixType >
void DiscreteCellProblemOperator< DiscreteFunctionImp, DiffusionImp >::assembleCellRHS_linear(const DomainType& x_T,
  const JacobianRangeType& gradient_PHI_H,
  DiscreteFunctionImp& cell_problem_RHS) const {
  typedef typename DiscreteFunction::DiscreteFunctionSpaceType DiscreteFunctionSpace;
  typedef typename DiscreteFunction::LocalFunctionType         LocalFunction;

  typedef typename DiscreteFunctionSpace::BaseFunctionSetType BaseFunctionSet;
  typedef typename DiscreteFunctionSpace::IteratorType        Iterator;
  typedef typename Iterator::Entity                           Entity;
  typedef typename Entity::Geometry                           Geometry;

  typedef typename DiscreteFunctionSpace::GridPartType GridPart;
  typedef CachingQuadrature< GridPart, 0 >             Quadrature;

  const DiscreteFunctionSpace& discreteFunctionSpace = cell_problem_RHS.space();

  // set entries to zero:
  cell_problem_RHS.clear();

  // get edge length of cell:
  const double delta = DSC_CONFIG_GET("problem.delta", 1.0f);

  // gradient of micro scale base function:
  std::vector< JacobianRangeType > gradient_phi( discreteFunctionSpace.mapper().maxNumDofs() );

//  RangeType rhs_L2_Norm = 0.0;

  const Iterator end = discreteFunctionSpace.end();
  for (Iterator it = discreteFunctionSpace.begin(); it != end; ++it)
  {
    const Entity& cell_grid_entity = *it;
    const Geometry& geometry = cell_grid_entity.geometry();
    assert(cell_grid_entity.partitionType() == InteriorEntity);

    LocalFunction elementOfRHS = cell_problem_RHS.localFunction(cell_grid_entity);

    const BaseFunctionSet& baseSet = elementOfRHS.baseFunctionSet();
    const unsigned int numBaseFunctions = baseSet.size();

    const Quadrature quadrature(cell_grid_entity, 2 * discreteFunctionSpace.order() + 2);
    const size_t numQuadraturePoints = quadrature.nop();
    for (size_t quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint)
    {
      const typename Quadrature::CoordinateType& local_point = quadrature.point(quadraturePoint);

      // global point in the unit cell Y:
      const DomainType global_point = geometry.global(local_point);

      // x_T + (delta * global_point)
      DomainType x_T_delta_global_point;
      for (int k = 0; k < dimension; ++k)
      {
        x_T_delta_global_point[k] = x_T[k] + (delta * global_point[k]);
      }

      const double weight = quadrature.weight(quadraturePoint) * geometry.integrationElement(local_point);

      // transposed of the the inverse jacobian
      const FieldMatrix< double, dimension, dimension >& inverse_jac
        = geometry.jacobianInverseTransposed(local_point);

      // A^{\eps}( x_T + \delta y) \nabla_x PHI_H(x_T)
      // diffusion operator evaluated in (x_T + \delta y) multiplied with \nabla_x PHI_H(x_T)
      JacobianRangeType diffusion_in_gradient_PHI_H;
      diffusion_operator_.diffusiveFlux(x_T_delta_global_point, gradient_PHI_H, diffusion_in_gradient_PHI_H);

      for (unsigned int i = 0; i < numBaseFunctions; ++i)
      {
        // jacobian of the base functions, with respect to the reference element
        JacobianRangeType gradient_phi_ref_element;
        baseSet.jacobian(i, quadrature[quadraturePoint], gradient_phi_ref_element);

        // multiply it with transpose of jacobian inverse to obtain the jacobian with respect to the real entity
        inverse_jac.mv(gradient_phi_ref_element[0], gradient_phi[i][0]);
      }

      for (unsigned int i = 0; i < numBaseFunctions; ++i)
      {
        elementOfRHS[i] -= weight * (diffusion_in_gradient_PHI_H[0] * gradient_phi[i][0]);
      }
    }
  }
} // assembleCellRHS_linear

template< class DiscreteFunctionImp, class DiffusionImp >
void DiscreteCellProblemOperator< DiscreteFunctionImp, DiffusionImp >::assembleCellRHS_nonlinear(const DomainType& x_T,
  const JacobianRangeType& grad_coarse_function,
  const DiscreteFunctionImp& old_fine_function,
  DiscreteFunctionImp& cell_problem_RHS) const {
  typedef typename DiscreteFunction::DiscreteFunctionSpaceType DiscreteFunctionSpace;
  typedef typename DiscreteFunction::LocalFunctionType         LocalFunction;

  typedef typename DiscreteFunctionSpace::BaseFunctionSetType BaseFunctionSet;
  typedef typename DiscreteFunctionSpace::IteratorType        Iterator;
  typedef typename Iterator::Entity                           Entity;
  typedef typename Entity::Geometry                           Geometry;

  typedef typename DiscreteFunctionSpace::GridPartType GridPart;
  typedef CachingQuadrature< GridPart, 0 >             Quadrature;

  const DiscreteFunctionSpace& discreteFunctionSpace = cell_problem_RHS.space();

  // set entries to zero:
  cell_problem_RHS.clear();

  // get edge length of cell:
  const double delta = DSC_CONFIG_GET("problem.delta", 1.0f);

  // gradient of micro scale base function:
  std::vector< JacobianRangeType > gradient_phi( discreteFunctionSpace.mapper().maxNumDofs() );

  const Iterator end = discreteFunctionSpace.end();
  for (Iterator it = discreteFunctionSpace.begin(); it != end; ++it)
  {
    const Entity& cell_grid_entity = *it;
    const Geometry& geometry = cell_grid_entity.geometry();
    assert(cell_grid_entity.partitionType() == InteriorEntity);

    LocalFunction local_old_fine_function = old_fine_function.localFunction(cell_grid_entity);
    LocalFunction elementOfRHS = cell_problem_RHS.localFunction(cell_grid_entity);

    const BaseFunctionSet& baseSet = elementOfRHS.baseFunctionSet();
    const unsigned int numBaseFunctions = baseSet.size();

    Quadrature quadrature(cell_grid_entity, 2 * discreteFunctionSpace.order() + 2);
    const size_t numQuadraturePoints = quadrature.nop();
    for (size_t quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint)
    {
      const typename Quadrature::CoordinateType& local_point = quadrature.point(quadraturePoint);

      // global point in the unit cell Y:
      DomainType global_point = geometry.global(local_point);

      // x_T + (delta * global_point)
      DomainType x_T_delta_global_point;
      for (int k = 0; k < dimension; ++k)
      {
        x_T_delta_global_point[k] = x_T[k] + (delta * global_point[k]);
      }

      JacobianRangeType grad_old_fine_function;
      local_old_fine_function.jacobian(quadrature[quadraturePoint], grad_old_fine_function);
      // here: no multiplication with jacobian inverse transposed required!

      JacobianRangeType position_vector;
      for (int k = 0; k < dimension; ++k)
      {
        position_vector[0][k] = grad_coarse_function[0][k] + grad_old_fine_function[0][k];
      }

      // A^{\eps}( x_T + \delta y, \nabla_x grad_coarse_function(x_T) + \nabla_y old_fine_function )
      JacobianRangeType diffusive_flux;
      diffusion_operator_.diffusiveFlux(x_T_delta_global_point, position_vector, diffusive_flux);

      const double weight = quadrature.weight(quadraturePoint) * geometry.integrationElement(local_point);

      // transposed of the the inverse jacobian
      const FieldMatrix< double, dimension, dimension >& inverse_jac
        = geometry.jacobianInverseTransposed(local_point);

      for (unsigned int i = 0; i < numBaseFunctions; ++i)
      {
        // jacobian of the base functions, with respect to the reference element
        JacobianRangeType gradient_phi_ref_element;
        baseSet.jacobian(i, quadrature[quadraturePoint], gradient_phi_ref_element);

        // multiply it with transpose of jacobian inverse to obtain the jacobian with respect to the real entity
        inverse_jac.mv(gradient_phi_ref_element[0], gradient_phi[i][0]);
      }

      for (unsigned int i = 0; i < numBaseFunctions; ++i)
      {
        elementOfRHS[i] -= weight * (diffusive_flux[0] * gradient_phi[i][0]);
      }
    }
  }
} // assembleCellRHS_nonlinear


template< class DiscreteFunctionImp, class DiffusionImp >
void DiscreteCellProblemOperator< DiscreteFunctionImp, DiffusionImp >::assemble_jacobian_corrector_cell_prob_RHS
  ( // the global quadrature point in the macro grid element T
  const DomainType& x_T, // barycenter of macro entity T
  // gradient of the old coarse function (old means last iteration step)
  const JacobianRangeType& grad_old_coarse_function, // \nabla_x u_H^{(n-1)}(x_T)
  // gradient of the corrector of the old coarse function
  const DiscreteFunctionImp& corrector_of_old_coarse_function, /*Q_h(u_H^{(n-1)})*/
  // gradient of the current macroscopic base function
  const JacobianRangeType& grad_coarse_base_function, // \nabla_x \Phi_H(x_T)
  // rhs cell problem:
  DiscreteFunctionImp& jac_corrector_cell_problem_RHS) const {
  // ! typedefs for the (discrete) periodic micro space:

  typedef typename DiscreteFunction::DiscreteFunctionSpaceType DiscreteFunctionSpace;
  typedef typename DiscreteFunction::LocalFunctionType         LocalFunction;

  typedef typename DiscreteFunctionSpace::BaseFunctionSetType BaseFunctionSet;
  typedef typename DiscreteFunctionSpace::IteratorType        Iterator;
  typedef typename Iterator::Entity                           Entity;
  typedef typename Entity::Geometry                           Geometry;

  // this is a periodic grid partition:
  typedef typename DiscreteFunctionSpace::GridPartType GridPart;
  typedef CachingQuadrature< GridPart, 0 >             Quadrature;

  const DiscreteFunctionSpace& discreteFunctionSpace = corrector_of_old_coarse_function.space();

  // set entries of right hand side to zero:
  jac_corrector_cell_problem_RHS.clear();

  const double delta = DSC_CONFIG_GET("problem.delta", 1.0f);

  // gradient of micro scale base function:
  std::vector< JacobianRangeType > gradient_phi( discreteFunctionSpace.mapper().maxNumDofs() );

  // iterator of micro (or cell) grid elements
  const Iterator end = discreteFunctionSpace.end();
  for (Iterator it = discreteFunctionSpace.begin(); it != end; ++it)
  {
    const Entity& cell_grid_entity = *it;
    const Geometry& geometry = cell_grid_entity.geometry();
    assert(cell_grid_entity.partitionType() == InteriorEntity);

    // local Q_h(u_H^{(n-1)}):
    const LocalFunction local_Q_old_u_H = corrector_of_old_coarse_function.localFunction(cell_grid_entity);
    LocalFunction elementOfRHS = jac_corrector_cell_problem_RHS.localFunction(cell_grid_entity);

    const BaseFunctionSet& baseSet = elementOfRHS.baseFunctionSet();
    const unsigned int numBaseFunctions = baseSet.size();

    const Quadrature quadrature(cell_grid_entity, 2 * discreteFunctionSpace.order() + 2);
    const size_t numQuadraturePoints = quadrature.nop();
    for (size_t quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint)
    {
      const typename Quadrature::CoordinateType& local_point = quadrature.point(quadraturePoint);

      // global point in the unit cell Y:
      const DomainType global_point_in_Y = geometry.global(local_point);

      // x_T + (delta * global_point_in_Y)
      DomainType x_T_plus_delta_y;
      for (int k = 0; k < dimension; ++k)
      {
        x_T_plus_delta_y[k] = x_T[k] + (delta * global_point_in_Y[k]);
      }

      // grad_y Q_h( u_H^{(n-1)})
      JacobianRangeType grad_Q_old_u_H;
      local_Q_old_u_H.jacobian(quadrature[quadraturePoint], grad_Q_old_u_H);
      // here: no multiplication with jacobian inverse transposed required!

      JacobianRangeType position_vector;
      for (int k = 0; k < dimension; ++k)
      {
        position_vector[0][k] = grad_old_coarse_function[0][k] + grad_Q_old_u_H[0][k];
      }

      JacobianRangeType direction_vector;
      for (int k = 0; k < dimension; ++k)
      {
        direction_vector[0][k] = grad_coarse_base_function[0][k];
      }

      // DA^{\eps}( x_T + \delta y, \nabla_x u_H^{(n-1)})(x_T) + \nabla_y Q_h( u_H^{(n-1)})(y) )( \nabla_x \Phi_H(x_T) )
      JacobianRangeType jacobian_diffusive_flux;
      diffusion_operator_.jacobianDiffusiveFlux(x_T_plus_delta_y,
                                                position_vector,
                                                direction_vector,
                                                jacobian_diffusive_flux);

      const double weight = quadrature.weight(quadraturePoint) * geometry.integrationElement(local_point);

      // transposed of the the inverse jacobian
      const auto& inverse_jac = geometry.jacobianInverseTransposed(local_point);
      baseSet.jacobianAll(quadrature[quadraturePoint], inverse_jac, gradient_phi);

      for (unsigned int i = 0; i < numBaseFunctions; ++i)
      {
        elementOfRHS[i] -= weight * (jacobian_diffusive_flux[0] * gradient_phi[i][0]);
      }
    }
  }
} // assemble_jacobian_corrector_cell_prob_RHS

} // namespace Dune

