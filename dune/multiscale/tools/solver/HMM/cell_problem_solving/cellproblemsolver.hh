#ifndef DiscreteEllipticCellProblem_HH
#define DiscreteEllipticCellProblem_HH

#include <dune/common/fmatrix.hh>

#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/operator/common/operator.hh>

// artificical mass coefficient to guarantee uniqueness and existence of the cell problem solution
// (should be as small as possible)
#define CELL_MASS_WEIGHT 0.0000001

// CELLSOLVER_VERBOSE: 0 = false, 1 = true
#define CELLSOLVER_VERBOSE false

#include <dune/fem/operator/2order/lagrangematrixsetup.hh>

#include <dune/multiscale/tools/disc_func_writer/discretefunctionwriter.hh>

namespace Dune {
// define output traits
struct CellProblemDataOutputParameters
  : public DataOutputParameters
{
public:
  std::string my_prefix_;
  std::string my_path_;

  void set_prefix(std::string my_prefix) {
    my_prefix_ = my_prefix;
    // std :: cout << "Set prefix. my_prefix_ = " << my_prefix_ << std :: endl;
  }

  void set_path(std::string my_path) {
    my_path_ = my_path;
  }

  // base of file name for data file
  std::string prefix() const {
    if (my_prefix_ == "")
      return "solutions";
    else
      return my_prefix_;
  }

  // path where the data is stored
  std::string path() const {
    if (my_path_ == "")
      return "data_output_hmm";
    else
      return my_path_;
  }

  // format of output:
  int outputformat() const {
    // return 0; // GRAPE (lossless format)
    return 1; // VTK
    // return 2; // VTK vertex data
    // return 3; // gnuplot
  }
};

// Imp stands for Implementation
template< class PeriodicDiscreteFunctionImp, class DiffusionImp >
class DiscreteCellProblemOperator
  : public Operator< typename PeriodicDiscreteFunctionImp::RangeFieldType,
                     typename PeriodicDiscreteFunctionImp::RangeFieldType, PeriodicDiscreteFunctionImp,
                     PeriodicDiscreteFunctionImp >
{
  typedef DiscreteCellProblemOperator< PeriodicDiscreteFunctionImp, DiffusionImp > This;

public:
  typedef PeriodicDiscreteFunctionImp DiscreteFunction;
  typedef DiffusionImp                DiffusionModel;

  typedef typename DiscreteFunction::DiscreteFunctionSpaceType DiscreteFunctionSpace;

  typedef typename DiscreteFunctionSpace::GridPartType   GridPart;
  typedef typename DiscreteFunctionSpace::GridType       GridType;
  typedef typename DiscreteFunctionSpace::RangeFieldType RangeFieldType;

  typedef typename DiscreteFunctionSpace::DomainType DomainType;
  typedef typename DiscreteFunctionSpace::RangeType  RangeType;
  typedef typename DiscreteFunctionSpace::JacobianRangeType
  JacobianRangeType;

protected:
  static const int dimension = GridPart::GridType::dimension;
  static const int polynomialOrder = DiscreteFunctionSpace::polynomialOrder;

  typedef typename DiscreteFunction::LocalFunctionType LocalFunction;

  typedef typename DiscreteFunctionSpace::BaseFunctionSetType                   BaseFunctionSet;
  typedef typename DiscreteFunctionSpace::LagrangePointSetType                  LagrangePointSet;
  typedef typename LagrangePointSet::template Codim< 1 >::SubEntityIteratorType FaceDofIterator;

  typedef typename DiscreteFunctionSpace::IteratorType Iterator;
  typedef typename Iterator::Entity                    Entity;
  typedef typename Entity::Geometry                    Geometry;

  typedef typename GridPart::IntersectionIteratorType IntersectionIterator;
  typedef typename IntersectionIterator::Intersection Intersection;

  typedef CachingQuadrature< GridPart, 0 > Quadrature;

public:
  DiscreteCellProblemOperator(const DiscreteFunctionSpace& periodicDiscreteFunctionSpace,
                              const DiffusionModel& diffusion_op)
    : periodicDiscreteFunctionSpace_(periodicDiscreteFunctionSpace)
      , diffusion_operator_(diffusion_op)
  {}

private:
  DiscreteCellProblemOperator(const This&);

public:
  // dummy operator
  virtual void operator()(const DiscreteFunction& u, DiscreteFunction& w) const;

  template< class MatrixType >
  void assemble_matrix(const DomainType& x_T, MatrixType& global_matrix) const;

  template< class MatrixType >
  void assemble_jacobian_matrix(const DomainType& x_T,
                                JacobianRangeType& grad_coarse_function,
                                DiscreteFunction& old_fine_function,
                                MatrixType& global_matrix) const;

  // the right hand side assembler methods
  void assembleCellRHS_linear(  // the global quadrature point in the macro grid element T
    const DomainType& x_T,
    // \nabla_x \Phi_H(x_T) (the coarse function to reconstruct):
    JacobianRangeType& grad_coarse_function,
    // rhs cell problem:
    DiscreteFunction& cell_problem_RHS) const;

  void assembleCellRHS_nonlinear(  // the global quadrature point in the macro grid element T
    const DomainType& x_T,
    // \nabla_x \Phi_H(x_T) :
    JacobianRangeType& grad_coarse_function,
    // old solution from the last iteration step
    DiscreteFunction& old_fine_function,
    // rhs cell problem:
    DiscreteFunction& cell_problem_RHS) const;

  // assemble method for the right hand side of the jacobian corrector cell problem
  void assemble_jacobian_corrector_cell_prob_RHS( // the global quadrature point in the macro grid element T
    const DomainType& x_T,
    // gradient of the old coarse function (old means last iteration step)
    JacobianRangeType& grad_old_coarse_function,
    // gradient of the corrector of the old coarse function
    DiscreteFunction& corrector_of_old_coarse_function,
    // gradient of the current macroscopic base function
    JacobianRangeType& grad_coarse_base_function,
    // rhs cell problem:
    DiscreteFunction& jac_corrector_cell_problem_RHS) const;

  void printCellRHS(DiscreteFunction& rhs) const;

  double normRHS(DiscreteFunction& rhs) const;

private:
  const DiscreteFunctionSpace& periodicDiscreteFunctionSpace_;
  const DiffusionModel& diffusion_operator_;
};

// dummy implementation of "operator()"
// 'w' = effect of the discrete operator on 'u'
template< class DiscreteFunctionImp, class DiffusionImp >
void DiscreteCellProblemOperator< DiscreteFunctionImp, DiffusionImp >::operator()(const DiscreteFunction& u,
                                                                                  DiscreteFunction& w) const {
  std::cout << "the ()-operator of the DiscreteCellProblemOperator class is not yet implemented and still a dummy."
            << std::endl;
  std::abort();
}

// ! stiffness matrix for a linear elliptic diffusion operator
// we obtain entries of the following kind
// (cell problem for the macro grid element 'T' and for the base-function '\Phi_H',
// x_T denotes the barycenter of T, \delta denotes the cell size )
// \int_Y A_h^{\eps}(t,x_T + \delta*y) \nabla phi_h_i(y) \cdot \nabla phi_h_j(y)
// + CELL_MASS_WEIGHT * \int_Y phi_h_i(y) \phi_h_j(y)
// (the second summand yields an artificical mass term to guarantee uniqueness and existence
// for the problem with periodic boundary condition.
// This is an alternative to the 'average zero' condition.)
template< class PeriodicDiscreteFunctionImp, class DiffusionImp >
template< class MatrixType >
void DiscreteCellProblemOperator< PeriodicDiscreteFunctionImp, DiffusionImp >::assemble_matrix(
  const DomainType& x_T,
  MatrixType&
  global_matrix) const
// x_T is the barycenter of the macro grid element T
{
  typedef typename MatrixType::LocalMatrixType LocalMatrix;

  Problem::ModelProblemData model_info;
  const double delta = model_info.getDelta();

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
    const unsigned int numBaseFunctions = baseSet.numBaseFunctions();

    // for constant diffusion "2*discreteFunctionSpace_.order()" is sufficient, for the general case, it is better to
    // use a higher order quadrature:
    Quadrature quadrature(cell_grid_entity, 2 * periodicDiscreteFunctionSpace_.order() + 2);
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

#if 1

// ! stiffness matrix for a non linear elliptic diffusion operator

// here we obtain the contribution of the jacobian matrix of the diffusion operator evaluated in the gradient of a
// certain discrete function (in case of the Newton method, it is the preceeding iterate v_h^{(n-1)} )

// Let PHI_H denote a macrocopic discrete function (that we want to reconstruct)
// we obtain entries of the following kind:
// (cell problem for the macro grid element 'T')

// \int_Y JA^{\eps}( x_T + \delta*y, \nabla_x PHI_H(x_T) + \nabla_y v_h^{(n-1)} ) \nabla phi_h_i(y) \cdot \nabla
// phi_h_j(y)
// + CELL_MASS_WEIGHT * \int_Y phi_h_i(y) \phi_h_j(y)

// (here, JA^{\eps} denotes the jacobian matrix of the diffusion operator A^{\eps},
// x_T denotes the barycenter of T, \delta denotes the cell size )

template< class PeriodicDiscreteFunctionImp, class DiffusionImp >
template< class MatrixType >
void DiscreteCellProblemOperator< PeriodicDiscreteFunctionImp, DiffusionImp >::assemble_jacobian_matrix(
  const DomainType& x_T,
  // macroscopic quadrature point
  JacobianRangeType
  & grad_coarse_function,
  // the gradient of the macroscopic function (that we want to reconstruct) evaluated in x_T
  PeriodicDiscreteFunctionImp
  & old_fine_function,
  // the microscopic function (fine-scale correction) from the last iteration step of the Newton
  // method
  MatrixType&
  global_matrix) const {
  Problem::ModelProblemData model_info;
  const double delta = model_info.getDelta();

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
    Quadrature quadrature(cell_grid_entity, 2 * periodicDiscreteFunctionSpace_.order() + 2);
    const size_t numQuadraturePoints = quadrature.nop();
    for (size_t quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint)
    {
      // local (barycentric) coordinates (with respect to entity)
      const typename Quadrature::CoordinateType& local_point = quadrature.point(quadraturePoint);

      // global point in the unit cell Y
      DomainType global_point = cell_grid_geometry.global(local_point);

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

#endif // if 1

#if 1
template< class DiscreteFunctionImp, class DiffusionImp >
void DiscreteCellProblemOperator< DiscreteFunctionImp, DiffusionImp >::printCellRHS(DiscreteFunctionImp& rhs) const {
  typedef typename DiscreteFunctionImp::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionSpaceType::IteratorType        IteratorType;
  typedef typename DiscreteFunctionImp::LocalFunctionType         LocalFunctionType;

  const DiscreteFunctionSpaceType& discreteFunctionSpace
    = rhs.space();

  const IteratorType endit = discreteFunctionSpace.end();
  for (IteratorType it = discreteFunctionSpace.begin(); it != endit; ++it)
  {
    LocalFunctionType elementOfRHS = rhs.localFunction(*it);

    const int numDofs = elementOfRHS.numDofs();
    for (int i = 0; i < numDofs; ++i)
    {
      std::cout << "Number of Dof: " << i << " ; " << rhs.name() << " : " << elementOfRHS[i] << std::endl;
    }
  }
}  // end method

template< class DiscreteFunctionImp, class DiffusionImp >
double DiscreteCellProblemOperator< DiscreteFunctionImp, DiffusionImp >::normRHS(DiscreteFunctionImp& rhs) const {
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
    CachingQuadrature< GridPartType, 0 > quadrature(entity, 2 * discreteFunctionSpace.order() + 2);

    // get geoemetry of entity
    const EnGeometryType& geo = entity.geometry();

    LocalFunctionType localRHS = rhs.localFunction(*it);

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

#endif // if 1

#if 1
// assemble the right hand side of a cell problem
// ----------------------------------------------

// assemble method for the case of a linear diffusion operator
// (in this case, no Newton method is required, which is why there is no dependency on an old fine-scale discrete
// function / old iteration step )

// we compute the following entries for each fine-scale base function phi_h_i:
// - \int_Y A^{\eps}( x_T + \delta*y ) \nabla_x PHI_H(x_T) \cdot \nabla_y phi_h_i(y)
template< class DiscreteFunctionImp, class DiffusionImp >
// template< class MatrixType >
void DiscreteCellProblemOperator< DiscreteFunctionImp, DiffusionImp >::assembleCellRHS_linear
  ( // the global quadrature point in the macro grid element T
  const DomainType& x_T,
  // \nabla_x \Phi_H(x_T):
  JacobianRangeType& gradient_PHI_H,
  // rhs cell problem:
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

  // model problem data:
  Problem::ModelProblemData problem_info;

  // get edge length of cell:
  const double delta = problem_info.getDelta();

  // gradient of micro scale base function:
  std::vector< JacobianRangeType > gradient_phi( discreteFunctionSpace.mapper().maxNumDofs() );

  RangeType rhs_L2_Norm = 0.0;

  const Iterator end = discreteFunctionSpace.end();
  for (Iterator it = discreteFunctionSpace.begin(); it != end; ++it)
  {
    const Entity& cell_grid_entity = *it;
    const Geometry& geometry = cell_grid_entity.geometry();
    assert(cell_grid_entity.partitionType() == InteriorEntity);

    LocalFunction elementOfRHS = cell_problem_RHS.localFunction(cell_grid_entity);

    const BaseFunctionSet& baseSet = elementOfRHS.baseFunctionSet();
    const unsigned int numBaseFunctions = baseSet.numBaseFunctions();

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

#endif // if 1

#if 1

// assemble method for the case of a nonlinear diffusion operator
// (in this case a Newton method is required, which is why there is an additional dependency on an old fine-scale
// discrete function v_h^{(n-1)} / old iteration step )

// we compute the following entries for each fine-scale base function phi_h_i:
// - \int_Y A^{\eps}( x_T + \delta*y , \nabla_x PHI_H(x_T) + \nabla_y \nabla_y v_h^{(n-1)} ) \cdot \nabla_y phi_h_i(y)
template< class DiscreteFunctionImp, class DiffusionImp >
void DiscreteCellProblemOperator< DiscreteFunctionImp, DiffusionImp >::assembleCellRHS_nonlinear
  ( // the global quadrature point in the macro grid element T
  const DomainType& x_T,
  // in the linear setting we typically reconstruct macroscopic base functions, in the non-linear setting we are in a
  // more general setting.
  // gradient of the coarse function, that we want to reconstruct:
  JacobianRangeType& grad_coarse_function,
  // old solution from the last iteration step
  DiscreteFunctionImp& old_fine_function,
  // rhs cell problem:
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

  // model problem data:
  Problem::ModelProblemData problem_info;

  // get edge length of cell:
  const double delta = problem_info.getDelta();

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

#endif // if 1

#if 1

// assemble method for the right hand side of the jacobian corrector cell problem
// (details, see 'solve_jacobiancorrector_cellproblem' below)
template< class DiscreteFunctionImp, class DiffusionImp >
void DiscreteCellProblemOperator< DiscreteFunctionImp, DiffusionImp >::assemble_jacobian_corrector_cell_prob_RHS
  ( // the global quadrature point in the macro grid element T
  const DomainType& x_T, // barycenter of macro entity T
  // gradient of the old coarse function (old means last iteration step)
  JacobianRangeType& grad_old_coarse_function, // \nabla_x u_H^{(n-1)}(x_T)
  // gradient of the corrector of the old coarse function
  DiscreteFunctionImp& corrector_of_old_coarse_function, /*Q_h(u_H^{(n-1)})*/
  // gradient of the current macroscopic base function
  JacobianRangeType& grad_coarse_base_function, // \nabla_x \Phi_H(x_T)
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

  // model problem data:
  Problem::ModelProblemData problem_info;

  // get edge length of cell:
  const double delta = problem_info.getDelta();

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
    LocalFunction local_Q_old_u_H = corrector_of_old_coarse_function.localFunction(cell_grid_entity);
    LocalFunction elementOfRHS = jac_corrector_cell_problem_RHS.localFunction(cell_grid_entity);

    const BaseFunctionSet& baseSet = elementOfRHS.baseFunctionSet();
    const unsigned int numBaseFunctions = baseSet.size();

    Quadrature quadrature(cell_grid_entity, 2 * discreteFunctionSpace.order() + 2);
    const size_t numQuadraturePoints = quadrature.nop();
    for (size_t quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint)
    {
      const typename Quadrature::CoordinateType& local_point = quadrature.point(quadraturePoint);

      // global point in the unit cell Y:
      DomainType global_point_in_Y = geometry.global(local_point);

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
        elementOfRHS[i] -= weight * (jacobian_diffusive_flux[0] * gradient_phi[i][0]);
      }
    }
  }
} // assemble_jacobian_corrector_cell_prob_RHS

#endif // if 1

// ! ------------------------------------------------------------------------------------------------
// ! ------------------------------------------------------------------------------------------------

// ! ------------------------------------------------------------------------------------------------
// ! ---------------- the cell problem numbering manager classes ------------------------------------

// comparison class for the CellProblemNumberingManager:
template< class GridPartType, class DomainType, class EntityPointerType >
struct classcomp
{
  bool operator()(const std::pair< EntityPointerType, int >& left_entity_pair,
                  const std::pair< EntityPointerType, int >& right_entity_pair) const {
    // compare the barycenteres of the entities with the lexicographic order, than compare the int's (number of local
    // base function)

    typedef CachingQuadrature< GridPartType, 0 > Quadrature;

    // ------ right element

    const typename EntityPointerType::Entity::Geometry& geometry_right = ( *(right_entity_pair.first) ).geometry();

    Quadrature quadrature_right( ( *(right_entity_pair.first) ), 0 );

    // local barycenter (with respect to entity)
    const typename Quadrature::CoordinateType& local_point_right = quadrature_right.point(0);

    DomainType barycenter_right_entity = geometry_right.global(local_point_right);

    // ------ left element

    const typename EntityPointerType::Entity::Geometry& geometry_left = ( *(left_entity_pair.first) ).geometry();

    Quadrature quadrature_left( ( *(left_entity_pair.first) ), 0 );

    // local barycenter (with respect to entity)
    const typename Quadrature::CoordinateType& local_point_left = quadrature_left.point(0);

    DomainType barycenter_left_entity = geometry_left.global(local_point_left);

    enum { dimension = GridPartType::GridType::dimension };

    int current_axis = dimension - 1;

    while (current_axis >= 0)
    {
      if (barycenter_left_entity[current_axis] < barycenter_right_entity[current_axis])
      { return true; } else if (barycenter_left_entity[current_axis] > barycenter_right_entity[current_axis])
      { return false; }

      current_axis -= 1;
    }

    if (left_entity_pair.second < right_entity_pair.second)
    {
      return true;
    } else
    { return false; }

    return true;
  } // ()
};

// comparison class for the CellProblemNumberingManager (just comparison of two entities!)
template< class GridPartType, class DomainType, class EntityPointerType >
struct entity_compare
{
  bool operator()(EntityPointerType left_entity,
                  EntityPointerType right_entity) const {
    // compare the barycenteres of the entities with the lexicographic order

    typedef CachingQuadrature< GridPartType, 0 > Quadrature;

    // ------ right element

    const typename EntityPointerType::Entity::Geometry& geometry_right = (*right_entity).geometry();

    Quadrature quadrature_right(*right_entity, 0);

    // local barycenter (with respect to entity)
    const typename Quadrature::CoordinateType& local_point_right = quadrature_right.point(0);

    DomainType barycenter_right_entity = geometry_right.global(local_point_right);

    // ------ left element

    const typename EntityPointerType::Entity::Geometry& geometry_left = (*left_entity).geometry();

    Quadrature quadrature_left(*left_entity, 0);

    // local barycenter (with respect to entity)
    const typename Quadrature::CoordinateType& local_point_left = quadrature_left.point(0);

    DomainType barycenter_left_entity = geometry_left.global(local_point_left);

    enum { dimension = GridPartType::GridType::dimension };

    int current_axis = dimension - 1;

    while (current_axis >= 0)
    {
      if (barycenter_left_entity[current_axis] < barycenter_right_entity[current_axis])
      { return true; } else if (barycenter_left_entity[current_axis] > barycenter_right_entity[current_axis])
      { return false; }

      current_axis -= 1;
    }

    return false;
  } // ()
};

// only for the combination entity + number of local base function on entity
template< class DiscreteFunctionSpaceType >
class CellProblemNumberingManager
{
public:
  typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
  typedef typename GridPartType::GridType                  GridType;

  typedef typename GridType::template Codim< 0 >::Entity         EntityType;
  typedef typename GridType::template  Codim< 0 >::EntityPointer EntityPointerType;

  typedef typename DiscreteFunctionSpaceType::BaseFunctionSetType BaseFunctionSetType;
  typedef typename DiscreteFunctionSpaceType::IteratorType        IteratorType;

  typedef typename DiscreteFunctionSpaceType::DomainType DomainType;

  typedef classcomp< GridPartType, DomainType, EntityPointerType > CompClass;

  typedef std::map< std::pair< EntityPointerType, int >, int, CompClass > CellNumMapType;

  // for the comparison of two entities:
  typedef entity_compare< GridPartType, DomainType, EntityPointerType > CompEntityClass;

  typedef std::map< EntityPointerType, int, CompEntityClass > CellNumMapNLType;

  CellNumMapType* cell_numbering_map_;
  CellNumMapNLType* cell_numbering_map_NL_;

  // simpliefied: in general we need CellNumMapType for the cell problem numering in the linear setting (entity and
  // local number of base function) and in the nonlinear case we need CellNumMapNLType (NL stands for nonlinear).
  // CellNumMapType is also required in the nonlinear case if we use test function reconstruction (TFR)

  inline explicit CellProblemNumberingManager(DiscreteFunctionSpaceType& discreteFunctionSpace) {
    cell_numbering_map_ = new CellNumMapType;
    cell_numbering_map_NL_ = new CellNumMapNLType;

    int counter = 0;
    int number_of_entity = 0;

    IteratorType endit = discreteFunctionSpace.end();
    for (IteratorType it = discreteFunctionSpace.begin(); it != endit; ++it)
    {
      cell_numbering_map_NL_->insert( std::make_pair(EntityPointerType(*it), number_of_entity) );

      const BaseFunctionSetType baseSet
        = discreteFunctionSpace.baseFunctionSet(*it);

      // number of base functions on entity
      const int numBaseFunctions = baseSet.size();

      for (int i = 0; i < numBaseFunctions; ++i)
      {
        std::pair< EntityPointerType, int > idPair(EntityPointerType(*it), i);
        cell_numbering_map_->insert( std::make_pair(idPair, counter) );
        counter++;
      }

      number_of_entity++;
    }
  }

  // use 'cp_num_manager.get_number_of_cell_problem( it, i )'
  inline int get_number_of_cell_problem(EntityPointerType& ent, const int& numOfBaseFunction) const {
    std::pair< EntityPointerType, int > idPair(ent, numOfBaseFunction);
    // !TODO this can create elements
    return (*cell_numbering_map_)[idPair];
  }

  // use 'cp_num_manager.get_number_of_cell_problem( it )'
  // Note: 'get_number_of_cell_problem( it )' is NOT equal to 'get_number_of_cell_problem( it , 0 )'!
  inline int get_number_of_cell_problem(EntityPointerType& ent) const {
    // !TODO this can create elements
    return (*cell_numbering_map_NL_)[ent];
  }
};

// ! ------------------ end of the cell problem numbering manager classes ---------------------------
// ! ------------------------------------------------------------------------------------------------

// ! ------------------------------------------------------------------------------------------------

// ! ------------------------------------------------------------------------------------------------
// ! --------------------- the essential cell problem solver class ----------------------------------

template< class PeriodicDiscreteFunctionImp, class DiffusionOperatorImp >
class CellProblemSolver
{
public:
  // ! type of discrete functions
  typedef PeriodicDiscreteFunctionImp PeriodicDiscreteFunctionType;

  // ! type of discrete function space
  typedef typename PeriodicDiscreteFunctionType::DiscreteFunctionSpaceType
  PeriodicDiscreteFunctionSpaceType;

  // ! type of grid partition
  typedef typename PeriodicDiscreteFunctionSpaceType::GridPartType PeriodicGridPartType;

  // ! type of grid
  typedef typename PeriodicDiscreteFunctionSpaceType::GridType PeriodicGridType;

  // ! type of range vectors
  typedef typename PeriodicDiscreteFunctionSpaceType::RangeType RangeType;

  // ! type of range vectors
  typedef typename PeriodicDiscreteFunctionSpaceType::DomainType DomainType;

  // ! polynomial order of base functions
  enum { polynomialOrder = PeriodicDiscreteFunctionSpaceType::polynomialOrder };

  // ! type of the (possibly non-linear) diffusion operator
  typedef DiffusionOperatorImp DiffusionType;

  struct CellMatrixTraits
  {
    typedef PeriodicDiscreteFunctionSpaceType                          RowSpaceType;
    typedef PeriodicDiscreteFunctionSpaceType                          ColumnSpaceType;
    typedef LagrangeMatrixSetup< false >                               StencilType;
    typedef ParallelScalarProduct< PeriodicDiscreteFunctionSpaceType > ParallelScalarProductType;

    template< class M >
    struct Adapter
    {
      typedef LagrangeParallelMatrixAdapter< M > MatrixAdapterType;
    };
  };

  typedef SparseRowMatrixOperator< PeriodicDiscreteFunctionType, PeriodicDiscreteFunctionType,
                                   CellMatrixTraits > CellFEMMatrix;

  // OEMGMRESOp //OEMBICGSQOp // OEMBICGSTABOp
  typedef OEMBICGSTABOp< PeriodicDiscreteFunctionType, CellFEMMatrix > InverseCellFEMMatrix;

  // discrete elliptic operator describing the elliptic cell problems
  typedef DiscreteCellProblemOperator< PeriodicDiscreteFunctionType, DiffusionType > CellProblemOperatorType;

private:
  const PeriodicDiscreteFunctionSpaceType& periodicDiscreteFunctionSpace_; // Referenz &, wenn & verwendet, dann unten:
  DiffusionType& diffusion_;

  std::ofstream* data_file_;

public:
  // ! constructor - with diffusion operator A^{\epsilon}(x)
  CellProblemSolver(const PeriodicDiscreteFunctionSpaceType& periodicDiscreteFunctionSpace,
                    DiffusionType& diffusion_operator)
    : periodicDiscreteFunctionSpace_(periodicDiscreteFunctionSpace)
      , diffusion_(diffusion_operator)
      , data_file_(NULL)
  {}

  // ! constructor - with diffusion operator A^{\epsilon}(x)
  CellProblemSolver(const PeriodicDiscreteFunctionSpaceType& periodicDiscreteFunctionSpace,
                    DiffusionType& diffusion_operator,
                    std::ofstream& data_file)
    : periodicDiscreteFunctionSpace_(periodicDiscreteFunctionSpace)
      , diffusion_(diffusion_operator)
      , data_file_(&data_file)
  {}

  // ! ----------- method: solve cell problem ------------------------------------------

  template< class JacobianRangeImp >
  void solvecellproblem(JacobianRangeImp& gradient_PHI_H,
                        // the barycenter x_T of a macro grid element 'T'
                        const DomainType& globalQuadPoint,
                        PeriodicDiscreteFunctionType& cell_problem_solution) {
    // set solution equal to zero:
    cell_problem_solution.clear();

    // ! the matrix in our linear system of equations
    // in the non-linear case, it is the matrix for each iteration step
    CellFEMMatrix cell_system_matrix("Cell Problem System Matrix",
                                     periodicDiscreteFunctionSpace_,
                                     periodicDiscreteFunctionSpace_);

    // ! define the discrete (elliptic) cell problem operator
    // ( effect of the discretized differential operator on a certain discrete function )
    CellProblemOperatorType cell_problem_op(periodicDiscreteFunctionSpace_, diffusion_);

    // ! right hand side vector of the algebraic cell problem
    // (in the non-linear setting it changes for every iteration step)
    PeriodicDiscreteFunctionType cell_problem_rhs("rhs of cell problem", periodicDiscreteFunctionSpace_);
    cell_problem_rhs.clear();

    // NOTE:
    // is the right hand side of the cell problem equal to zero or almost identical to zero?
    // if yes, the solution of the cell problem is also identical to zero. The solver is getting a problem with this
    // situation, which is why we do not solve cell problems for zero-right-hand-side, since we already know the result.

    #ifdef LINEAR_PROBLEM

    // assemble the stiffness matrix
    cell_problem_op.assemble_matrix(globalQuadPoint, cell_system_matrix);

    // assemble right hand side of algebraic cell problem
    cell_problem_op.assembleCellRHS_linear(globalQuadPoint, gradient_PHI_H, cell_problem_rhs);

    const double norm_rhs = cell_problem_op.normRHS(cell_problem_rhs);

    if ( !( cell_problem_rhs.dofsValid() ) )
    {
      std::cout << "Cell Problem RHS invalid." << std::endl;
      abort();
    }

    if (norm_rhs < /*1e-06*/ 1e-10)
    {
      cell_problem_solution.clear();
      // std :: cout << "Cell problem with solution zero." << std :: endl;
    } else {
      InverseCellFEMMatrix cell_fem_biCGStab(cell_system_matrix, 1e-8, 1e-8, 20000, CELLSOLVER_VERBOSE);
      cell_fem_biCGStab(cell_problem_rhs, cell_problem_solution);
    }

    // end linear case.
    #else // ifdef LINEAR_PROBLEM
    // nonlinear case:

    // starting value for the Newton method
    PeriodicDiscreteFunctionType zero_func(" constant zero function ", periodicDiscreteFunctionSpace_);
    zero_func.clear();

    // residual vector (current residual)
    PeriodicDiscreteFunctionType cell_problem_residual("Cell Problem Residual", periodicDiscreteFunctionSpace_);
    cell_problem_residual.clear();

    L2Error< PeriodicDiscreteFunctionType > l2error;
    RangeType relative_newton_error = 10000.0;

    int iteration_step = 0;

    // the Newton step for for solving the current cell problem (solved with Newton Method):
    // L2-Norm of residual < tolerance ?
    #ifdef STOCHASTIC_PERTURBATION
    double tolerance = 1e-01 * VARIANCE;
    #else
    double tolerance = 1e-06;
    #endif // ifdef STOCHASTIC_PERTURBATION
    while (relative_newton_error > tolerance)
    {
      // (here: cellproblem_solution = solution from the last iteration step)

      // assemble the stiffness matrix
      cell_problem_op.assemble_jacobian_matrix(globalQuadPoint,
                                               gradient_PHI_H,
                                               cell_problem_solution,
                                               cell_system_matrix);

      // assemble right hand side of algebraic cell problem (for the current iteration step)
      cell_problem_op.assembleCellRHS_nonlinear(globalQuadPoint,
                                                gradient_PHI_H,
                                                cell_problem_solution,
                                                cell_problem_rhs);

      const double norm_rhs = cell_problem_op.normRHS(cell_problem_rhs);

      if ( !( cell_problem_rhs.dofsValid() ) )
      {
        std::cout << "Cell Problem RHS invalid." << std::endl;
        abort();
      }

      if ( (norm_rhs < /*1e-06*/ 1e-10) /*|| ( (norm_rhs > 1e-06) && (norm_rhs < 1.5e-06) )*/ )
      {
        // residual solution almost identical to zero: break
        break;
      }

      double biCG_tolerance = 1e-8;
      bool cell_solution_convenient = false;
      while (cell_solution_convenient == false)
      {
        cell_problem_residual.clear();
        InverseCellFEMMatrix cell_fem_newton_biCGStab(cell_system_matrix,
                                                      1e-8, biCG_tolerance, 20000, CELLSOLVER_VERBOSE);

        cell_fem_newton_biCGStab(cell_problem_rhs, cell_problem_residual);

        if ( cell_problem_residual.dofsValid() )
        { cell_solution_convenient = true; }

        if (biCG_tolerance > 1e-4)
        {
          std::cout << "WARNING! Iteration step " << iteration_step
                    << ". Invalid dofs in 'cell_problem_residual', but '" << relative_newton_error
                    <<
          " = relative_newton_error > 1e-01' and 'biCG_tolerance > 1e-4'. L^2-Norm of right hand side of cell problem: "
                    << norm_rhs << ". Therefore possibly inaccurate solution." << std::endl;
          std::cout << "Information:" << std::endl;
          std::cout << "x_T = globalQuadPoint = " << globalQuadPoint << "." << std::endl;
          std::cout << "nabla u_H^{(n-1)} = gradient_PHI_H = " << gradient_PHI_H[0] << "." << std::endl;
          std::cout << "Print right hand side? y/n: ";
          char answer;
          std::cin >> answer;
          if ( !(answer == 'n') )
          { cell_problem_op.printCellRHS(cell_problem_rhs); }
          std::abort();
        }

        biCG_tolerance *= 10.0;
      }

      cell_problem_solution += cell_problem_residual;

      relative_newton_error = l2error.template norm2< (2* polynomialOrder) + 2 >(cell_problem_residual, zero_func);
      RangeType norm_cell_solution = l2error.template norm2< (2* polynomialOrder) + 2 >(cell_problem_solution,
                                                                                        zero_func);
      relative_newton_error = relative_newton_error / norm_cell_solution;

      // std :: cout << "L2-Norm of cell problem residual = " << residual_L2_norm << std :: endl;

      cell_problem_residual.clear();

      if (iteration_step > 10)
      {
        std::cout << "Warning! Algorithm already reached Newton-iteration step " << iteration_step
                  << " for computing the nonlinear cellproblem." << std::endl;
        std::cout << "relative_newton_error = " << relative_newton_error << std::endl << std::endl;
        #ifdef FORCE
        residual_L2_norm = 0.0;
        #endif
      }

      iteration_step += 1;
    }

    #endif // ifdef LINEAR_PROBLEM
    // end nonlinear case.

    if ( !( cell_problem_solution.dofsValid() ) )
    {
      std::cout << "Current solution of the cell problem invalid!" << std::endl;
      std::abort();
    }
  } // solvecellproblem

  // ! ----------- end method: solve cell problem ------------------------------------------

  // ! -------- method: solve the jacobian corrector cell problem ----------------------------
  // Problem to determine the Jacobian operator of the correction operator
  // ( the correction operator Q is defined via cell problems, here we determine the corresponding derivative DQ )
  // (see paper what it means and where it comes from. Note that it is only required for the nonlinear case)
  template< class JacobianRangeImp >
  void solve_jacobiancorrector_cellproblem(
    // gradient of macroscopic base function
    JacobianRangeImp& gradient_PHI_H,
    // gradient of the macroscopic function from the last iteration step
    JacobianRangeImp& grad_old_coarse_function,
    // gradient_y of the corrector of the macroscopic function from the last iteration step
    PeriodicDiscreteFunctionType& corrector_of_old_coarse_function,
    // the barycenter x_T of a macro grid element 'T'
    const DomainType& globalQuadPoint,
    PeriodicDiscreteFunctionType& jac_cor_cell_problem_solution) {
    // set solution equal to zero:
    jac_cor_cell_problem_solution.clear();

    // ! the matrix in our linear system of equations
    // system matrix for the jacobian corrector cell problem to solve:
    // entries:
    // - int_Y DA^{\eps}( x_T + \delta y, \nabla_x u_H^{(n-1)})(x_T) + \nabla_y Q_h( u_H^{(n-1)})(y) ) \nablay_y
    // \phi_h_i(y) \nablay_y \phi_h_j(y) dy
    // ( \phi_h_i and \phi_h_j denote microscopic base functions.)
    CellFEMMatrix jac_cor_cell_system_matrix("Jacobian Corrector Cell Problem System Matrix",
                                             periodicDiscreteFunctionSpace_,
                                             periodicDiscreteFunctionSpace_);

    // ! define the discrete (elliptic) cell problem operator
    // ( effect of the discretized differential operator on a certain discrete function )
    CellProblemOperatorType cell_problem_op(periodicDiscreteFunctionSpace_, diffusion_);
    // we are looking for the derivative of the operator in \nabla_x u_H^{(n-1)})(x_T) + \nabla_y Q_h( u_H^{(n-1)})(y)
    // and in direction of the macroscopic base function

    // ! right hand side vector of the algebraic jacobian corrector cell problem
    // entries of the right hand side vector:
    // - int_Y DA^{\eps}( x_T + \delta y, \nabla_x u_H^{(n-1)})(x_T) + \nabla_y Q_h( u_H^{(n-1)})(y) )( \nabla_x
    // \Phi_H(x_T) ) \nablay_y \phi_h(y) dy
    // \Phi_H denotes the current macroscopic base function and
    // \phi_h denotes the current microscopic base function.
    PeriodicDiscreteFunctionType jac_cor_cell_problem_rhs("rhs of jacobian corrector cell problem",
                                                          periodicDiscreteFunctionSpace_);
    jac_cor_cell_problem_rhs.clear();

    // assemble the stiffness matrix
    cell_problem_op.assemble_jacobian_matrix(globalQuadPoint,
                                             grad_old_coarse_function,
                                             corrector_of_old_coarse_function,
                                             jac_cor_cell_system_matrix);

    // assemble right hand side of algebraic jacobian corrector cell problem
    cell_problem_op.assemble_jacobian_corrector_cell_prob_RHS
      (globalQuadPoint,
      grad_old_coarse_function,
      corrector_of_old_coarse_function,
      gradient_PHI_H,
      jac_cor_cell_problem_rhs);

    const double norm_rhs = cell_problem_op.normRHS(jac_cor_cell_problem_rhs);

    if ( !( jac_cor_cell_problem_rhs.dofsValid() ) )
    {
      std::cout << "Jacobian Corrector Cell Problem RHS invalid." << std::endl;
      std::abort();
    }

    // is the right hand side of the jacobian corrector cell problem equal to zero or almost identical to zero?
    // if yes, the solution of the cell problem is also identical to zero. The solver is getting a problem with this
    // situation, which is why we do not solve cell problems for zero-right-hand-side, since we already know the result.
    if (norm_rhs < 1e-10)
    {
      jac_cor_cell_problem_solution.clear();
      // std :: cout << "Jacobian Corrector Cell Problem with solution zero." << std :: endl;
    } else {
      InverseCellFEMMatrix jac_cor_cell_fem_biCGStab(jac_cor_cell_system_matrix, 1e-8, 1e-8, 20000, CELLSOLVER_VERBOSE);
      jac_cor_cell_fem_biCGStab(jac_cor_cell_problem_rhs, jac_cor_cell_problem_solution);
    }
  } // solve_jacobiancorrector_cellproblem

  // ! ---------- method: solve the jacobian corrector cell problem -------------------------

  // two methods for solving and saving the solutions of the cell problems.
  // 1. save solutions for the whole set of macroscopic base function
  // (in general for the case of a linear diffusion operator)
  // 2. save the solutions for a fixed discrete function
  // (in general for the case of a nonlinear diffusion operator)

  // ! ---- method: solve and save the cell problems for the set of macroscopic base functions -----

  // here we need a 'cell problem numbering manager' to determine the number of the cell problem
  // (a combination of number of entity and number of local base function)
  // Structure:
  // Struktur der Indizierung fuer das Abspeichern der Loesungen der Zellprobleme:
  // wir loesen Zellprobleme fuer jede Entity des Makro-Grids und jede Basisfunktion, die einen nichtleeren support auf
  // dieser Entity besitzt, also schematisch:
  // Sei n=0,...,N einer Durchnumerierung der Entitys und i=0,...I_n eine zu einer festen Entity gehoerende Nummerierung
  // der Basisfunktionen mit nicht-leeren support.
  // Die Durchnummerierung der Loesungen der Zellprobleme k=0,...,K ist dann gegeben durch: k(n,i_n) = ( sum_(l=0)^(n-1)
  // ( I_l + 1) ) + i_n
  // NOTE: es verhaelt sich NICHT wie die vorhandene Methode mapToGlobal(entity,i) ! (die gibt die globale Nummer der
  // Basisfunktion zurueck, es gibt aber  deutlich mehr Zellprobleme zum Loesen!
  // (das wird aber alles im Hintergrund vom 'cell problem numbering manager')

  // compute and save solutions of the cell problems for the base function set of the 'discreteFunctionSpace'
  // requires cell problem numbering manager
  template< class DiscreteFunctionImp, class CellProblemNumberingManagerImp >
  void saveTheSolutions_baseSet(
    const typename DiscreteFunctionImp::DiscreteFunctionSpaceType& discreteFunctionSpace,
    const CellProblemNumberingManagerImp& cp_num_manager,   // just to check, if we use the correct numeration
    const std::string& filename) {
    typedef DiscreteFunctionImp DiscreteFunctionType;

    typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType
    DiscreteFunctionSpaceType;

    typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;

    typedef typename DiscreteFunctionSpaceType::GridType GridType;

    typedef typename DiscreteFunctionSpaceType::JacobianRangeType
    JacobianRangeType;

    typedef typename DiscreteFunctionSpaceType::BaseFunctionSetType
    BaseFunctionSetType;

    typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;

    typedef typename DiscreteFunctionType::LocalFunctionType
    LocalFunctionType;

    typedef typename GridType::template Codim< 0 >::Entity EntityType;
    typedef typename EntityType::Geometry                  EntityGeometryType;

    typedef CachingQuadrature< GridPartType, 0 > EntityQuadratureType;

    enum { dimension = GridType::dimension };
    enum { maxnumOfBaseFct = 100 };

    bool writer_is_open = false;

    std::string cell_solution_location = "data/HMM/" + filename + "_cellSolutions_baseSet";
    DiscreteFunctionWriter dfw( (cell_solution_location).c_str() );

    writer_is_open = dfw.open();

    long double starting_time = clock();

    // we want to determine minimum, average and maxiumum time for solving a cell problem in the current method
    double minimum_time_c_p = 1000000;
    double average_time_c_p = 0;
    double maximum_time_c_p = 0;

    int number_of_cell_problem = 0;

    if (writer_is_open)
    {
      IteratorType endit = discreteFunctionSpace.end();
      for (IteratorType it = discreteFunctionSpace.begin(); it != endit; ++it)
      {
        // gradients of the macroscopic base functions
        JacobianRangeType gradientPhi[maxnumOfBaseFct];

        // entity
        const EntityType& entity = *it;

        const BaseFunctionSetType baseSet
          = discreteFunctionSpace.baseFunctionSet(entity);

        const EntityGeometryType& geometry = entity.geometry();

        EntityQuadratureType quadrature(entity, 0);

        DomainType barycenter_of_entity = geometry.global( quadrature.point(0) );

        // number of base functions on entity
        const int numBaseFunctions = baseSet.numBaseFunctions();

        // calc Jacobian inverse before volume is evaluated
        const FieldMatrix< double, dimension, dimension >& inv
          = geometry.jacobianInverseTransposed( quadrature.point(0 /*=quadraturePoint*/) );

        PeriodicDiscreteFunctionType correctorPhi_i("corrector Phi_i", periodicDiscreteFunctionSpace_);

        for (int i = 0; i < numBaseFunctions; ++i)
        {
          baseSet.jacobian(i, quadrature[0 /*=quadraturePoint*/], gradientPhi[i]);
          // multiply with transpose of jacobian inverse
          gradientPhi[i][0] = FMatrixHelp::mult(inv, gradientPhi[i][0]);
        }

        for (int i = 0; i < numBaseFunctions; ++i)
        {
          correctorPhi_i.clear();

          // take time
          long double time_now = clock();

          solvecellproblem< JacobianRangeType >
            (gradientPhi[i], barycenter_of_entity, correctorPhi_i);

          // min/max time
          if ( (clock() - time_now) / CLOCKS_PER_SEC > maximum_time_c_p )
          { maximum_time_c_p = (clock() - time_now) / CLOCKS_PER_SEC; }
          if ( (clock() - time_now) / CLOCKS_PER_SEC < minimum_time_c_p )
          { minimum_time_c_p = (clock() - time_now) / CLOCKS_PER_SEC; }

          dfw.append(correctorPhi_i);

          // check if we use a correct numeration of the cell problems:
          if ( !(cp_num_manager.get_number_of_cell_problem(it, i) == number_of_cell_problem) )
          {
            std::cout << "Numeration of cell problems incorrect." << std::endl;
            std::abort();
          }

          number_of_cell_problem++;
        }
      } // end: for-loop: IteratorType it
    } // end: 'if ( writer_is_open )'

    if (data_file_)
    {
      if ( data_file_->is_open() )
      {
        (*data_file_) << std::endl;
        (*data_file_) << "In method: saveTheSolutions_baseSet." << std::endl << std::endl;
        (*data_file_) << "Cell problems solved for " << discreteFunctionSpace.grid().size(0) << " leaf entities."
                      << std::endl;
        (*data_file_) << "Minimum time for solving a cell problem = " << minimum_time_c_p << "s." << std::endl;
        (*data_file_) << "Maximum time for solving a cell problem = " << maximum_time_c_p << "s." << std::endl;
        (*data_file_) << "Average time for solving a cell problem = "
                      << ( (clock() - starting_time) / CLOCKS_PER_SEC ) / number_of_cell_problem << "s." << std::endl;
        (*data_file_) << "Total time for computing and saving the cell problems = "
                      << ( (clock() - starting_time) / CLOCKS_PER_SEC ) << "s," << std::endl << std::endl;
      }
    }
  } // saveTheSolutions_baseSet

  // ! ---- method: solve and save the cell problems for a fixed macroscopic discrete function -----

  // this method does not require a cell problem numbering manager
  // (it uses the standard counter for entities, provided by the corrsponding iterator)

  // compute and save solutions of the cell problems for a fixed macroscopic discrete function
  // (in gerneral it is the macro solution from the last iteration step)
  template< class DiscreteFunctionImp >
  void saveTheSolutions_discFunc(
    const DiscreteFunctionImp& macro_discrete_function,
    const std::string& filename) {
    typedef DiscreteFunctionImp DiscreteFunctionType;

    typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType
    DiscreteFunctionSpaceType;

    typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;

    typedef typename DiscreteFunctionSpaceType::GridType GridType;

    typedef typename DiscreteFunctionSpaceType::JacobianRangeType
    JacobianRangeType;

    typedef typename DiscreteFunctionSpaceType::BaseFunctionSetType
    BaseFunctionSetType;

    typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;

    typedef typename DiscreteFunctionType::LocalFunctionType
    LocalFunctionType;

    typedef typename GridType::template Codim< 0 >::Entity EntityType;
    typedef typename EntityType::Geometry                  EntityGeometryType;

    typedef CachingQuadrature< GridPartType, 0 > EntityQuadratureType;

    enum { dimension = GridType::dimension };
    enum { maxnumOfBaseFct = 100 };

    bool writer_is_open = false;

    std::string cell_solution_location = "data/HMM/" + filename + "_cellSolutions_discFunc";
    DiscreteFunctionWriter dfw( (cell_solution_location).c_str() );

    writer_is_open = dfw.open();

    long double starting_time = clock();

    // we want to determine minimum, average and maxiumum time for solving a cell problem in the current method
    double minimum_time_c_p = 1000000;
    double DUNE_UNUSED(average_time_c_p) = 0;
    double maximum_time_c_p = 0;

    const DiscreteFunctionSpaceType& discreteFunctionSpace = macro_discrete_function.space();

    if (writer_is_open)
    {
      int number_of_entity = 0;
      IteratorType endit = discreteFunctionSpace.end();
      for (IteratorType it = discreteFunctionSpace.begin(); it != endit; ++it)
      {
        // entity
        const EntityType& entity = *it;

        const EntityGeometryType& geometry = entity.geometry();

        EntityQuadratureType quadrature(entity, 0);

        DomainType barycenter_of_entity = geometry.global( quadrature.point(0) );

        LocalFunctionType local_macro_disc = macro_discrete_function.localFunction(entity);
        JacobianRangeType grad_macro_discrete_function;
        local_macro_disc.jacobian(quadrature[0], grad_macro_discrete_function);

        PeriodicDiscreteFunctionType cell_solution_on_entity("corrector of macro discrete function",
                                                             periodicDiscreteFunctionSpace_);

        // take time
        long double time_now = clock();

        solvecellproblem< JacobianRangeType >
          (grad_macro_discrete_function, barycenter_of_entity, cell_solution_on_entity);

        // min/max time
        if ( (clock() - time_now) / CLOCKS_PER_SEC > maximum_time_c_p )
        { maximum_time_c_p = (clock() - time_now) / CLOCKS_PER_SEC; }
        if ( (clock() - time_now) / CLOCKS_PER_SEC < minimum_time_c_p )
        { minimum_time_c_p = (clock() - time_now) / CLOCKS_PER_SEC; }

        dfw.append(cell_solution_on_entity);

        #if 0
        // !LOESCHEN:
        if (number_of_entity == 341)
        {
          // in case you want to save the solutions of the two cell problems:
          #if 1
          typedef Tuple< PeriodicDiscreteFunctionImp* > IOTupleType;
          typedef DataOutput< GridType, IOTupleType >   DataOutputType;

          // general output parameters
          CellProblemDataOutputParameters outputparam;

          // sequence stamp
          std::stringstream outstring;

          // ------- cell problem -------------

          // create and initialize output class
          IOTupleType cellproblem_tuple(&cell_solution_on_entity);

          outputparam.set_prefix("cellSolution_saved_");
          outputparam.set_path("data/HMM/");
          DataOutputType cellSolution_dataoutput(periodicDiscreteFunctionSpace_.grid(), cellproblem_tuple, outputparam);

          // write data
          outstring << "cellSolution_saved_";
          cellSolution_dataoutput.writeData( 1.0 /*dummy*/, outstring.str() );
          // clear the std::stringstream:
          outstring.str( std::string() );

          #endif // if 1
        }
        // !-------------------------------------
        #endif // if 0

        number_of_entity += 1;
      } // end: for-loop: IteratorType it
    } // end: 'if ( writer_is_open )'

    if (data_file_)
    {
      if ( data_file_->is_open() )
      {
        (*data_file_) << std::endl;
        (*data_file_) << "In method: saveTheSolutions_discFunc." << std::endl << std::endl;
        (*data_file_) << "Cell problems solved for " << discreteFunctionSpace.grid().size(0) << " leaf entities."
                      << std::endl;
        (*data_file_) << "Minimum time for solving a cell problem = " << minimum_time_c_p << "s." << std::endl;
        (*data_file_) << "Maximum time for solving a cell problem = " << maximum_time_c_p << "s." << std::endl;
        (*data_file_) << "Average time for solving a cell problem = "
                      << ( (clock()
              - starting_time) / CLOCKS_PER_SEC ) / discreteFunctionSpace.grid().size(0) << "s." << std::endl;
        (*data_file_) << "Total time for computing and saving the cell problems = "
                      << ( (clock() - starting_time) / CLOCKS_PER_SEC ) << "s," << std::endl << std::endl;
      }
    }
  } // saveTheSolutions_discFunc

  #if 1

  // compute and save solutions of the jacobian corrector cell problems for the base function set of the
  // 'discreteFunctionSpace' and for a fixed macroscopic discrete function
  // requires cell problem numbering manager
  template< class DiscreteFunctionImp, class CellProblemNumberingManagerImp >
  void saveTheJacCorSolutions_baseSet_discFunc(
    const DiscreteFunctionImp& macro_discrete_function,
    const CellProblemNumberingManagerImp& cp_num_manager,   // just to check, if we use the correct numeration
    const std::string& filename) {
    typedef DiscreteFunctionImp DiscreteFunctionType;

    typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType
    DiscreteFunctionSpaceType;

    typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;

    typedef typename DiscreteFunctionSpaceType::GridType GridType;

    typedef typename DiscreteFunctionSpaceType::JacobianRangeType
    JacobianRangeType;

    typedef typename DiscreteFunctionSpaceType::BaseFunctionSetType
    BaseFunctionSetType;

    typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;

    typedef typename DiscreteFunctionType::LocalFunctionType
    LocalFunctionType;

    typedef typename GridType::template Codim< 0 >::Entity EntityType;
    typedef typename EntityType::EntityPointer
    EntityPointerType;
    typedef typename EntityType::Geometry EntityGeometryType;

    typedef CachingQuadrature< GridPartType, 0 > EntityQuadratureType;

    enum { dimension = GridType::dimension };
    enum { maxnumOfBaseFct = 100 };

    bool writer_is_open = false;

    // where we save the solutions:
    std::string cell_solution_location = "data/HMM/" + filename + "_JacCorCellSolutions_baseSet_discFunc";
    DiscreteFunctionWriter dfw( (cell_solution_location).c_str() );
    writer_is_open = dfw.open();

    // where we saved the solutions for the discrete function
    // NOTE: they already need to be assembled, i.e. we already applied the method saveSolutions_discFunc!
    std::string cell_solution_discFunc_location = "data/HMM/" + filename + "_cellSolutions_discFunc";

    bool reader_is_open = false;

    // reader for the cell problem data file (discrete functions):
    DiscreteFunctionReader discrete_function_reader( (cell_solution_discFunc_location).c_str() );
    reader_is_open = discrete_function_reader.open();

    long double starting_time = clock();

    // we want to determine minimum, average and maxiumum time for solving a cell problem in the current method
    double minimum_time_c_p = 1000000;
    double average_time_c_p = 0;
    double maximum_time_c_p = 0;

    int number_of_cell_problem = 0;

    const DiscreteFunctionSpaceType& discreteFunctionSpace = macro_discrete_function.space();

    if (writer_is_open)
    {
      int number_of_entity = 0;
      IteratorType endit = discreteFunctionSpace.end();
      for (IteratorType it = discreteFunctionSpace.begin(); it != endit; ++it)
      {
        // gradients of the macroscopic base functions
        JacobianRangeType gradientPhi[maxnumOfBaseFct];

        // entity
        const EntityType& entity = *it;

        const BaseFunctionSetType baseSet
          = discreteFunctionSpace.baseFunctionSet(entity);

        const EntityGeometryType& geometry = entity.geometry();

        EntityQuadratureType quadrature(entity, 0);

        DomainType barycenter_of_entity = geometry.global( quadrature.point(0) );

        // number of base functions on entity
        const int numBaseFunctions = baseSet.size();

        // calc Jacobian inverse before volume is evaluated
        const FieldMatrix< double, dimension, dimension >& inv
          = geometry.jacobianInverseTransposed( quadrature.point(0 /*=quadraturePoint*/) );

        LocalFunctionType local_macro_disc = macro_discrete_function.localFunction(entity);
        JacobianRangeType grad_macro_discrete_function;
        local_macro_disc.jacobian(quadrature[0], grad_macro_discrete_function);

        PeriodicDiscreteFunctionType corrector_macro_discrete_function("corrector of macro discrete function",
                                                                       periodicDiscreteFunctionSpace_);
        if (reader_is_open)
        { discrete_function_reader.read(number_of_entity, corrector_macro_discrete_function); }

        // the solution that we want to save to the data file
        PeriodicDiscreteFunctionType jac_corrector_Phi_i("jacobian corrector of Phi_i", periodicDiscreteFunctionSpace_);

        for (int i = 0; i < numBaseFunctions; ++i)
        {
          baseSet.jacobian(i, quadrature[0 /*=quadraturePoint*/], gradientPhi[i]);
          // multiply with transpose of jacobian inverse
          gradientPhi[i][0] = FMatrixHelp::mult(inv, gradientPhi[i][0]);
        }

        for (int i = 0; i < numBaseFunctions; ++i)
        {
          jac_corrector_Phi_i.clear();

          // take time
          long double time_now = clock();

          solve_jacobiancorrector_cellproblem< JacobianRangeType >
            (gradientPhi[i],
            grad_macro_discrete_function,
            corrector_macro_discrete_function,
            barycenter_of_entity,
            jac_corrector_Phi_i);

          // min/max time
          if ( (clock() - time_now) / CLOCKS_PER_SEC > maximum_time_c_p )
          { maximum_time_c_p = (clock() - time_now) / CLOCKS_PER_SEC; }
          if ( (clock() - time_now) / CLOCKS_PER_SEC < minimum_time_c_p )
          { minimum_time_c_p = (clock() - time_now) / CLOCKS_PER_SEC; }

          dfw.append(jac_corrector_Phi_i);

          // check if we use a correct numeration of the cell problems:
          EntityPointerType entity_pointer(*it);
          if ( !(cp_num_manager.get_number_of_cell_problem(entity_pointer, i) == number_of_cell_problem) )
          {
            std::cout << "Numeration of cell problems incorrect." << std::endl;
            std::abort();
          }

          number_of_cell_problem++;
        }

        number_of_entity += 1;
      } // end: for-loop: IteratorType it
    } // end: 'if ( writer_is_open )'

    if (data_file_)
    {
      if ( data_file_->is_open() )
      {
        (*data_file_) << std::endl;
        (*data_file_) << "In method: saveTheJacCorSolutions_baseSet_discFunc." << std::endl << std::endl;
        (*data_file_) << "Cell problems solved for " << discreteFunctionSpace.grid().size(0) << " leaf entities."
                      << std::endl;
        (*data_file_) << "Minimum time for solving a cell problem = " << minimum_time_c_p << "s." << std::endl;
        (*data_file_) << "Maximum time for solving a cell problem = " << maximum_time_c_p << "s." << std::endl;
        (*data_file_) << "Average time for solving a cell problem = "
                      << ( (clock() - starting_time) / CLOCKS_PER_SEC ) / number_of_cell_problem << "s." << std::endl;
        (*data_file_) << "Total time for computing and saving the cell problems = "
                      << ( (clock() - starting_time) / CLOCKS_PER_SEC ) << "s," << std::endl << std::endl;
      }
    }
  } // saveTheJacCorSolutions_baseSet_discFunc

  #endif // if 1
}; // end class
}

#endif // #ifndef DiscreteElliptic_HH
