#ifndef DiscreteEllipticMsFEMLocalProblem_HH
#define DiscreteEllipticMsFEMLocalProblem_HH

#include <dune/fem/operator/matrix/spmatrix.hh>
#include <dune/common/fmatrix.hh>

#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/operator/common/operator.hh>

// LOCPROBLEMSOLVER_VERBOSE: 0 = false, 1 = true
#define LOCPROBLEMSOLVER_VERBOSE false

// write solutions of the local problems (vtk)
// #define LOCALDATAOUTPUT

#include <dune/fem/operator/2order/lagrangematrixsetup.hh>

#include <dune/multiscale/tools/disc_func_writer/discretefunctionwriter.hh>

// the local problem:
// Let 'T' denote a coarse grid element and
// let 'T_0' denote the reference element in 2D
// ( the rectangular triangle with corners (0,0), (1,0) and (0,1) )
// let 'Fx = Ax + b' denote the corresponding affine transformation with ' F(T) = T_0 '
// Now, for fixed 'T' and for fixed coarse grid basis function PHI_H,
// we define the operaor Q^eps(\Phi_H) \in H^1(T) and zero boundary by:
// \int_T A^eps(x) ∇Q^eps(\Phi_H)(x) · ∇ \phi(x) + \int_T A^eps(x) ∇ \Phi_H(x_T) · ∇ \phi(x) = 0
// for all \phi \in \mathring{H}^1(T).
// Applying the transformation formula, we get:
// \int_{T_0} (A^eps ○ F)(x) ( ∇Q^eps(\Phi_H)) ○ F )(x) · ∇ \phi(x)
// + \int_{T_0} (A^eps ○ F)(x) ∇ \Phi_H(x_T) · ∇ \phi(x) = 0
// This yields that (Q^eps(\Phi_H) ○ F) is the solution of
// \int_{T_0} (A^eps ○ F)(x) (A^{-1})^T ∇( Q^eps(\Phi_H) ○ F )(x) · ( (A^{-1})^T ∇ \phi(x) )
// + \int_{T_0} (A^eps ○ F)(x) ∇ \Phi_H(x_T) · ( (A^{-1})^T ∇ \phi(x) ) = 0
//
// Here, (A^{-1})^T denotes the transposed of the inverse of the matrix 'A' of F(x)=Ax+b

// Let (a_0,a_1,a_2) deonte the corners of the 2-simplex T, then the matrix A in the affine transformation
// F(x) = Ax + a_0, F : T_0 -> T is given by
// A_11 = a_1(1) - a_0(1)     A_12 = a_2(1) - a_0(1)
// A_21 = a_1(2) - a_0(2)     A_22 = a_2(2) - a_0(2)

// define 'c := (a_1(1) - a_0(1))·(a_2(2) - a_0(2)) - (a_1(2) - a_0(2))·(a_2(1) - a_0(1))
// then the inverse A^{-1} is given by:
// A^{-1}_11 = (1/c) (a_2(2) - a_0(2))     A^{-1}_12 = (1/c) (a_0(1) - a_2(1))
// A^{-1}_21 = (1/c) (a_0(2) - a_1(2))     A^{-1}_22 = (1/c) (a_1(1) - a_0(1))

// jacobian-inverse-transposed

namespace Dune {
// define output traits
struct LocalProblemDataOutputParameters
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
    return 1;   // VTK
    // return 2; // VTK vertex data
    // return 3; // gnuplot
  }
};

// Imp stands for Implementation
template< class GlobalEntityDiscreteFunctionImp, class DiffusionImp >
class LocalProblemOperator
  : public Operator< typename GlobalEntityDiscreteFunctionImp::RangeFieldType,
                     typename GlobalEntityDiscreteFunctionImp::RangeFieldType, GlobalEntityDiscreteFunctionImp,
                     GlobalEntityDiscreteFunctionImp >
{
  typedef LocalProblemOperator< GlobalEntityDiscreteFunctionImp, DiffusionImp > This;

public:
  typedef GlobalEntityDiscreteFunctionImp DiscreteFunction;
  typedef DiffusionImp                    DiffusionModel;

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
  LocalProblemOperator(const DiscreteFunctionSpace& globalEntityDiscreteFunctionSpace,
                       const DiffusionModel& diffusion_op)
    : globalEntityDiscreteFunctionSpace_(globalEntityDiscreteFunctionSpace)
      , diffusion_operator_(diffusion_op)
  {}

private:
  LocalProblemOperator(const This&);

public:
  // dummy operator
  virtual void operator()(const DiscreteFunction& u, DiscreteFunction& w) const;

  template< class MatrixType >
  void assemble_matrix(Iterator& pointer_on_entity, MatrixType& global_matrix) const;

  // the right hand side assembler methods
  void assemble_local_RHS(    // pointer on the coarse grid element T
    Iterator& pointer_on_entity,
    // \nabla_x \Phi_H(x_T) (the coarse function to reconstruct):
    JacobianRangeType& grad_coarse_function,
    // rhs local msfem problem:
    DiscreteFunction& local_problem_RHS) const;

  void printLocalRHS(DiscreteFunction& rhs) const;

  double normRHS(DiscreteFunction& rhs) const;

private:
  const DiscreteFunctionSpace& globalEntityDiscreteFunctionSpace_;
  const DiffusionModel& diffusion_operator_;
};

// dummy implementation of "operator()"
// 'w' = effect of the discrete operator on 'u'
template< class DiscreteFunctionImp, class DiffusionImp >
void LocalProblemOperator< DiscreteFunctionImp, DiffusionImp >::operator()(const DiscreteFunction& u,
                                                                           DiscreteFunction& w) const {
  std::cout << "the ()-operator of the LocalProblemOperator class is not yet implemented and still a dummy."
            << std::endl;
  std::abort();
}

// ! stiffness matrix for a linear elliptic diffusion operator
// we obtain entries of the following kind
// (local msfem problem for the macro grid element 'T' and for the base-function '\Phi_H',
// T_0 denotes the "reference element" in 2D
// x_T denotes the barycenter of T
// We solve for (Q^eps(\Phi_H) ○ F), which is the solution of
// \int_{T_0} (A^eps ○ F)(x) (A^{-1})^T ∇( Q^eps(\Phi_H) ○ F )(x) · ∇ \phi(x)
// + \int_{T_0} (A^eps ○ F)(x) ∇ \Phi_H(x_T) · ∇ \phi(x) = 0
// In the method "assemble_matrix", we are only concerened with
// \int_{T_0} (A^eps ○ F)(x) (A^{-1})^T ∇( Q^eps(\Phi_H) ○ F )(x) · ∇ \phi(x)
template< class GlobalEntityDiscreteFunctionImp, class DiffusionImp >
template< class MatrixType >
void LocalProblemOperator< GlobalEntityDiscreteFunctionImp, DiffusionImp >::assemble_matrix(Iterator& pointer_on_entity,
                                                                                            MatrixType& global_matrix)
const
// x_T is the barycenter of the macro grid element T
{
  // the coarse grid element T:
  const Entity& entity_T = *pointer_on_entity;
  const Geometry& geometry_of_T = entity_T.geometry();

  // transformation F : T_0 -> T
  // describe the mapping F(x) = Ax + b with F(T_0)=T for an entity T and the reference element T_0:
  // arguments: entity T, point in T_0, point in T.

  // Let (a_0,a_1,a_2) deonte the corners of the 2-simplex T, then the matrix A in the affine transformation
  // F(x) = Ax + a_0, F : T_0 -> T is given by
  // A_11 = a_1( 1 ) - a_0( 1 )     A_12 = a_2( 1 ) - a_0( 1 )
  // A_21 = a_1( 2 ) - a_0( 2 )     A_22 = a_2( 2 ) - a_0( 2 )

  // corners of the reference element:
  typename Quadrature::CoordinateType ref_corner_0, ref_corner_1, ref_corner_2;

  ref_corner_0[0] = 0.0;
  ref_corner_0[1] = 0.0;

  ref_corner_1[0] = 1.0;
  ref_corner_1[1] = 0.0;

  ref_corner_2[0] = 0.0;
  ref_corner_2[1] = 1.0;

  // corner of the global element:
  DomainType corner_0_of_T = geometry_of_T.global(ref_corner_0);
  DomainType corner_1_of_T = geometry_of_T.global(ref_corner_1);
  DomainType corner_2_of_T = geometry_of_T.global(ref_corner_2);

  // std :: cout << "corner_0_of_T = " << corner_0_of_T << std :: endl;
  // std :: cout << "corner_1_of_T = " << corner_1_of_T << std :: endl;
  // std :: cout << "corner_2_of_T = " << corner_2_of_T << std :: endl;

  // value of the matrix A (in F(x) = Ax + a_0)
  double val_A[dimension][dimension];
  val_A[0][0] = corner_1_of_T[0] - corner_0_of_T[0];
  val_A[0][1] = corner_2_of_T[0] - corner_0_of_T[0];
  val_A[1][0] = corner_1_of_T[1] - corner_0_of_T[1];
  val_A[1][1] = corner_2_of_T[1] - corner_0_of_T[1];

  // define 'c := (a_1(1) - a_0(1))·(a_2(2) - a_0(2)) - (a_1(2) - a_0(2))·(a_2(1) - a_0(1))
  double c = 1.0 / ( (val_A[0][0] * val_A[1][1]) - (val_A[0][1] * val_A[1][0]) );
  // then the inverse A^{-1} is given by:
  // A^{-1}_11 = (1/c) (a_2(2) - a_0(2))     A^{-1}_12 = (1/c) (a_0(1) - a_2(1))
  // A^{-1}_21 = (1/c) (a_0(2) - a_1(2))     A^{-1}_22 = (1/c) (a_1(1) - a_0(1))

  // (A^{-1})^T:
  double val_A_inverse_transposed[dimension][dimension];
  val_A_inverse_transposed[0][0] = c * val_A[1][1];
  val_A_inverse_transposed[1][0] = c * (-1.0) * val_A[0][1];
  val_A_inverse_transposed[0][1] = c * (-1.0) * val_A[1][0];
  val_A_inverse_transposed[1][1] = c * val_A[0][0];

  typedef typename MatrixType::LocalMatrixType LocalMatrix;

  Problem::ModelProblemData model_info;

  global_matrix.reserve();
  global_matrix.clear();

  // micro scale base function:
  std::vector< RangeType > phi( globalEntityDiscreteFunctionSpace_.mapper().maxNumDofs() );

  // gradient of micro scale base function:
  std::vector< typename BaseFunctionSet::JacobianRangeType > gradient_phi(
    globalEntityDiscreteFunctionSpace_.mapper().maxNumDofs() );

  const Iterator end = globalEntityDiscreteFunctionSpace_.end();
  for (Iterator it = globalEntityDiscreteFunctionSpace_.begin(); it != end; ++it)
  {
    const Entity& local_grid_entity = *it;
    const Geometry& local_grid_geometry = local_grid_entity.geometry();
    assert(local_grid_entity.partitionType() == InteriorEntity);

    LocalMatrix local_matrix = global_matrix.localMatrix(local_grid_entity, local_grid_entity);

    const BaseFunctionSet& baseSet = local_matrix.domainBaseFunctionSet();
    const unsigned int numBaseFunctions = baseSet.numBaseFunctions();

    // for constant diffusion "2*discreteFunctionSpace_.order()" is sufficient, for the general case, it is better to
    // use a higher order quadrature:
    Quadrature quadrature(local_grid_entity, 2 * globalEntityDiscreteFunctionSpace_.order() + 2);
    const size_t numQuadraturePoints = quadrature.nop();
    for (size_t quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint)
    {
      // local (barycentric) coordinates (with respect to local grid entity)
      const typename Quadrature::CoordinateType& local_point = quadrature.point(quadraturePoint);
      // remember, we are concerned with: \int_{T_0} (A^eps ○ F)(x) (A^{-1})^T ∇( Q^eps(\Phi_H) ○ F )(x) · ∇
      // \phi(x)

      // global point in the reference element T_0
      DomainType global_point = local_grid_geometry.global(local_point);

      // F ( global point in the reference element T_0 )
      // (the transformation of the global point in T_0 to its position in T)
      DomainType global_point_transformed(0.0);

      for (int k = 0; k < dimension; ++k)
        for (int l = 0; l < dimension; ++l)
          global_point_transformed[k] += (val_A[k][l] * global_point[l]);

      global_point_transformed += corner_0_of_T;   // = global_point;

      // F(x) = Ax + a_0, F : T_0 -> T is given by
      // A_11 = a_1(1) - a_0(1)     A_12 = a_2(1) - a_0(1)
      // A_21 = a_1(2) - a_0(2)     A_22 = a_2(2) - a_0(2)

      const double weight = quadrature.weight(quadraturePoint)
                            * local_grid_geometry.integrationElement(local_point);

      // transposed of the the inverse jacobian
      const FieldMatrix< double, dimension, dimension >& inverse_jac
        = local_grid_geometry.jacobianInverseTransposed(local_point);

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
        typename BaseFunctionSet::JacobianRangeType val_A_inverse_transposed_gradient_phi_i(0.0);
        for (int k = 0; k < dimension; ++k)
          for (int l = 0; l < dimension; ++l)
            val_A_inverse_transposed_gradient_phi_i[0][k] += val_A_inverse_transposed[k][l] * gradient_phi[i][0][l];

        // A( x_T + \delta y, \nabla \phi )
        // diffusion operator evaluated in (x_T + \delta y , \nabla \phi)
        typename LocalFunction::JacobianRangeType diffusion_in_A_inv_T_gradient_phi;
        diffusion_operator_.diffusiveFlux(global_point_transformed,
                                          val_A_inverse_transposed_gradient_phi_i,
                                          diffusion_in_A_inv_T_gradient_phi);
        for (unsigned int j = 0; j < numBaseFunctions; ++j)
        {
          typename BaseFunctionSet::JacobianRangeType val_A_inverse_transposed_gradient_phi_j(0.0);
          for (int k = 0; k < dimension; ++k)
            for (int l = 0; l < dimension; ++l)
              val_A_inverse_transposed_gradient_phi_j[0][k] += val_A_inverse_transposed[k][l] * gradient_phi[j][0][l];

          // stiffness contribution
          local_matrix.add( j, i, weight
                            * (diffusion_in_A_inv_T_gradient_phi[0] * val_A_inverse_transposed_gradient_phi_j[0]) );

          // mass contribution (just for stabilization!)
          // local_matrix.add( j, i, 0.00000001 * weight * (phi[ i ][ 0 ] * phi[ j ][ 0 ]) );
        }
      }
    }
  }
} // assemble_matrix

template< class DiscreteFunctionImp, class DiffusionImp >
void LocalProblemOperator< DiscreteFunctionImp, DiffusionImp >::printLocalRHS(DiscreteFunctionImp& rhs) const {
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
}      // end method

template< class DiscreteFunctionImp, class DiffusionImp >
double LocalProblemOperator< DiscreteFunctionImp, DiffusionImp >::normRHS(DiscreteFunctionImp& rhs) const {
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

// std :: cout << "value rhs = " << value << std::endl;

      norm += weight * value * value;
    }
  }

  return norm;
}      // end method

// assemble the right hand side of a local problem (reconstruction problem on entity)
// ----------------------------------------------

// assemble method for the case of a linear diffusion operator

// we compute the following entries for each fine-scale base function phi_h_i:
// - \int_{T_0} (A^eps ○ F)(x) ∇ \Phi_H(x_T) · ∇ \phi_h_i(x)
template< class DiscreteFunctionImp, class DiffusionImp >
// template< class MatrixType >
void LocalProblemOperator< DiscreteFunctionImp, DiffusionImp >::assemble_local_RHS
  (      // pointer in the macro grid element T
  Iterator& pointer_on_entity,
  // \nabla_x \Phi_H(x_T):
  JacobianRangeType& gradient_PHI_H,
  // rhs local msfem problem:
  DiscreteFunctionImp& local_problem_RHS) const {
  typedef typename DiscreteFunction::DiscreteFunctionSpaceType DiscreteFunctionSpace;
  typedef typename DiscreteFunction::LocalFunctionType         LocalFunction;

  typedef typename DiscreteFunctionSpace::BaseFunctionSetType BaseFunctionSet;
  typedef typename DiscreteFunctionSpace::IteratorType        Iterator;
  typedef typename Iterator::Entity                           Entity;
  typedef typename Entity::Geometry                           Geometry;

  typedef typename DiscreteFunctionSpace::GridPartType GridPart;
  typedef CachingQuadrature< GridPart, 0 >             Quadrature;

  const DiscreteFunctionSpace& discreteFunctionSpace = local_problem_RHS.space();

  // the coarse grid element T:
  const Entity& entity_T = *pointer_on_entity;
  const Geometry& geometry_of_T = entity_T.geometry();

  // transformation F : T_0 -> T
  // describe the mapping F(x) = Ax + b with F(T_0)=T for an entity T and the reference element T_0:
  // arguments: entity T, point in T_0, point in T.

  // Let (a_0,a_1,a_2) deonte the corners of the 2-simplex T, then the matrix A in the affine transformation
  // F(x) = Ax + a_0, F : T_0 -> T is given by
  // A_11 = a_1( 1 ) - a_0( 1 )     A_12 = a_2( 1 ) - a_0( 1 )
  // A_21 = a_1( 2 ) - a_0( 2 )     A_22 = a_2( 2 ) - a_0( 2 )

  // corners of the reference element:
  typename Quadrature::CoordinateType ref_corner_0, ref_corner_1, ref_corner_2;

  ref_corner_0[0] = 0.0;
  ref_corner_0[1] = 0.0;

  ref_corner_1[0] = 1.0;
  ref_corner_1[1] = 0.0;

  ref_corner_2[0] = 0.0;
  ref_corner_2[1] = 1.0;

  // corner of the global element:
  DomainType corner_0_of_T = geometry_of_T.global(ref_corner_0);
  DomainType corner_1_of_T = geometry_of_T.global(ref_corner_1);
  DomainType corner_2_of_T = geometry_of_T.global(ref_corner_2);

  // std :: cout << "corner_0_of_T = " << corner_0_of_T << std :: endl;
  // std :: cout << "corner_1_of_T = " << corner_1_of_T << std :: endl;
  // std :: cout << "corner_2_of_T = " << corner_2_of_T << std :: endl << std :: endl;

  // value of the matrix A (in F(x) = Ax + a_0)
  double val_A[dimension][dimension];
  val_A[0][0] = corner_1_of_T[0] - corner_0_of_T[0];
  val_A[0][1] = corner_2_of_T[0] - corner_0_of_T[0];
  val_A[1][0] = corner_1_of_T[1] - corner_0_of_T[1];
  val_A[1][1] = corner_2_of_T[1] - corner_0_of_T[1];

  // define 'c := (a_1(1) - a_0(1))·(a_2(2) - a_0(2)) - (a_1(2) - a_0(2))·(a_2(1) - a_0(1))
  double c = 1.0 / ( (val_A[0][0] * val_A[1][1]) - (val_A[0][1] * val_A[1][0]) );
  // then the inverse A^{-1} is given by:
  // A^{-1}_11 = (1/c) (a_2(2) - a_0(2))     A^{-1}_12 = (1/c) (a_0(1) - a_2(1))
  // A^{-1}_21 = (1/c) (a_0(2) - a_1(2))     A^{-1}_22 = (1/c) (a_1(1) - a_0(1))

  // (A^{-1})^T:
  double val_A_inverse_transposed[dimension][dimension];
  val_A_inverse_transposed[0][0] = c * val_A[1][1];
  val_A_inverse_transposed[1][0] = c * (-1.0) * val_A[0][1];
  val_A_inverse_transposed[0][1] = c * (-1.0) * val_A[1][0];
  val_A_inverse_transposed[1][1] = c * val_A[0][0];

  // set entries to zero:
  local_problem_RHS.clear();

  // model problem data:
  Problem::ModelProblemData problem_info;

  // gradient of micro scale base function:
  std::vector< JacobianRangeType > gradient_phi( discreteFunctionSpace.mapper().maxNumDofs() );

  RangeType rhs_L2_Norm = 0.0;

  const Iterator end = discreteFunctionSpace.end();
  for (Iterator it = discreteFunctionSpace.begin(); it != end; ++it)
  {
    const Entity& local_grid_entity = *it;
    const Geometry& geometry = local_grid_entity.geometry();
    assert(local_grid_entity.partitionType() == InteriorEntity);

    LocalFunction elementOfRHS = local_problem_RHS.localFunction(local_grid_entity);

    const BaseFunctionSet& baseSet = elementOfRHS.baseFunctionSet();
    const unsigned int numBaseFunctions = baseSet.numBaseFunctions();

    Quadrature quadrature(local_grid_entity, 2 * discreteFunctionSpace.order() + 2);
    const size_t numQuadraturePoints = quadrature.nop();
    for (size_t quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint)
    {
      const typename Quadrature::CoordinateType& local_point = quadrature.point(quadraturePoint);
      // remember, we are concerned with: - \int_{T_0} (A^eps ○ F)(x) ∇ \Phi_H(x_T) · ∇ \phi(x)
      // global point in the reference element T_0
      DomainType global_point = geometry.global(local_point);

      // F ( global point in the reference element T_0 )
      // (the transformation of the global point in T_0 to its position in T)
      DomainType global_point_transformed(0.0);

      for (int k = 0; k < dimension; ++k)
        for (int l = 0; l < dimension; ++l)
          global_point_transformed[k] += (val_A[k][l] * global_point[l]);

      global_point_transformed += corner_0_of_T;

      // F(x) = Ax + a_0, F : T_0 -> T is given by
      // A_11 = a_1(1) - a_0(1)     A_12 = a_2(1) - a_0(1)
      // A_21 = a_1(2) - a_0(2)     A_22 = a_2(2) - a_0(2)

      const double weight = quadrature.weight(quadraturePoint) * geometry.integrationElement(local_point);
      // transposed of the the inverse jacobian
      const FieldMatrix< double, dimension, dimension >& inverse_jac
        = geometry.jacobianInverseTransposed(local_point);

      // (A^eps ○ F)(x) ∇ \Phi_H(x_T)
      // diffusion operator evaluated in 'F(x)' multiplied with ∇ PHI_H(x_T)
      JacobianRangeType diffusion_in_gradient_PHI_H;
      diffusion_operator_.diffusiveFlux(global_point_transformed, gradient_PHI_H, diffusion_in_gradient_PHI_H);

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
        typename BaseFunctionSet::JacobianRangeType val_A_inverse_transposed_gradient_phi_i(0.0);
        for (int k = 0; k < dimension; ++k)
          for (int l = 0; l < dimension; ++l)
            val_A_inverse_transposed_gradient_phi_i[0][k] += val_A_inverse_transposed[k][l] * gradient_phi[i][0][l];

        elementOfRHS[i] -= weight * (diffusion_in_gradient_PHI_H[0] * val_A_inverse_transposed_gradient_phi_i[0]);
      }
    }
  }
} // assemble_local_RHS


// ! ------------------------------------------------------------------------------------------------
// ! ------------------------------------------------------------------------------------------------

// ! ------------------------------------------------------------------------------------------------
// ! ---------------- the local msfem problem numbering manager classes -----------------------------

// comparison class for the LocalProblemNumberingManager:
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

// comparison class for the LocalProblemNumberingManager (just comparison of two entities!)
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
class LocalProblemNumberingManager
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

  typedef std::map< std::pair< EntityPointerType, int >, int, CompClass > LocProbNumMapType;

  // for the comparison of two entities:
  typedef entity_compare< GridPartType, DomainType, EntityPointerType > CompEntityClass;

  typedef std::map< EntityPointerType, int, CompEntityClass > LocProbNumMapNLType;

private:
  // ! location of the solutions of the local problems
  std::string* location_;

  LocProbNumMapType* loc_prob_numbering_map_;
  LocProbNumMapNLType* loc_prob_numbering_map_NL_;

public:
  // space of the solutions of the local problems
  const DiscreteFunctionSpaceType& localDiscreteFunctionSpace_;

  // simpliefied: in general we need LocProbNumMapType for the local msfem problem numering in the linear setting
  // (entity and local number of base function) and in the nonlinear case we need LocProbNumMapNLType (NL stands for
  // nonlinear).
  // LocProbNumMapType is also required in the nonlinear case if we use the standard MsFEM formulation (no PGF)

public:
  inline explicit LocalProblemNumberingManager(const DiscreteFunctionSpaceType& discreteFunctionSpace) {
// : discreteFunctionSpace_( discreteFunctionSpace )
    loc_prob_numbering_map_ = new LocProbNumMapType;
    loc_prob_numbering_map_NL_ = new LocProbNumMapNLType;

    // discreteFunctionSpace_ = &discreteFunctionSpace;

    int counter = 0;
    int number_of_entity = 0;

    IteratorType endit = discreteFunctionSpace.end();
    for (IteratorType it = discreteFunctionSpace.begin(); it != endit; ++it)
    {
      loc_prob_numbering_map_NL_->insert( std::make_pair(it, number_of_entity) );

      const BaseFunctionSetType baseSet
        = discreteFunctionSpace.baseFunctionSet(*it);

      // number of base functions on entity
      const int numBaseFunctions = baseSet.numBaseFunctions();

      for (int i = 0; i < numBaseFunctions; ++i)
      {
        std::pair< EntityPointerType, int > idPair(it, i);
        loc_prob_numbering_map_->insert( std::make_pair(idPair, counter) );
        counter++;
      }

      number_of_entity++;
    }
  }

  inline explicit LocalProblemNumberingManager(const DiscreteFunctionSpaceType& discreteFunctionSpace,
                                               std::string& filename) {
// : discreteFunctionSpace_( discreteFunctionSpace )
    location_ = new std::string;
    *location_ = "data/MsFEM/" + filename + "/local_problems/_localProblemSolutions_baseSet";

    loc_prob_numbering_map_ = new LocProbNumMapType;
    loc_prob_numbering_map_NL_ = new LocProbNumMapNLType;

    // discreteFunctionSpace_ = &discreteFunctionSpace;

    int counter = 0;
    int number_of_entity = 0;

    IteratorType endit = discreteFunctionSpace.end();
    for (IteratorType it = discreteFunctionSpace.begin(); it != endit; ++it)
    {
      loc_prob_numbering_map_NL_->insert( std::make_pair(it, number_of_entity) );

      const BaseFunctionSetType baseSet
        = discreteFunctionSpace.baseFunctionSet(*it);

      // number of base functions on entity
      const int numBaseFunctions = baseSet.numBaseFunctions();

      for (int i = 0; i < numBaseFunctions; ++i)
      {
        std::pair< EntityPointerType, int > idPair(it, i);
        loc_prob_numbering_map_->insert( std::make_pair(idPair, counter) );
        counter++;
      }

      number_of_entity++;
    }
  }

  inline explicit LocalProblemNumberingManager(const DiscreteFunctionSpaceType& discreteFunctionSpace,
                                               const DiscreteFunctionSpaceType& localDiscreteFunctionSpace,
                                               std::string& filename)
    : localDiscreteFunctionSpace_(localDiscreteFunctionSpace) {
    location_ = new std::string;
    *location_ = "data/MsFEM/" + filename + "/local_problems/_localProblemSolutions_baseSet";

    loc_prob_numbering_map_ = new LocProbNumMapType;
    loc_prob_numbering_map_NL_ = new LocProbNumMapNLType;

    // discreteFunctionSpace_ = &discreteFunctionSpace;

    int counter = 0;
    int number_of_entity = 0;

    IteratorType endit = discreteFunctionSpace.end();
    for (IteratorType it = discreteFunctionSpace.begin(); it != endit; ++it)
    {
      loc_prob_numbering_map_NL_->insert( std::make_pair(it, number_of_entity) );

      const BaseFunctionSetType baseSet
        = discreteFunctionSpace.baseFunctionSet(*it);

      // number of base functions on entity
      const int numBaseFunctions = baseSet.numBaseFunctions();

      for (int i = 0; i < numBaseFunctions; ++i)
      {
        std::pair< EntityPointerType, int > idPair(it, i);
        loc_prob_numbering_map_->insert( std::make_pair(idPair, counter) );
        counter++;
      }

      number_of_entity++;
    }
  }

  // use 'lp_num_manager.get_number_of_local_problem( it, i )'
  inline int get_number_of_local_problem(EntityPointerType& ent, const int& numOfBaseFunction) const {
    std::pair< EntityPointerType, int > idPair(ent, numOfBaseFunction);
    return (*loc_prob_numbering_map_)[idPair];
  }

  // use 'lp_num_manager.get_number_of_local_problem( it )'
  // Note: 'get_number_of_local_problem( it )' is NOT equal to 'get_number_of_local_problem( it , 0 )'!
  inline int get_number_of_local_problem(EntityPointerType& ent) const {
    return (*loc_prob_numbering_map_NL_)[ent];
  }

  inline const DiscreteFunctionSpaceType& get_local_discrete_function_space() const {
    return localDiscreteFunctionSpace_;
  }

  inline std::string get_location() const {
    if (location_)
    { return *location_; } else
    { std::cout << "Warning! Location of the solutions of the local problems not specified!" << std::endl; }
  }
};

// ! ------------------ end of the local msfem problem numbering manager classes --------------------
// ! ------------------------------------------------------------------------------------------------

// ! ------------------------------------------------------------------------------------------------

// ! ------------------------------------------------------------------------------------------------
// ! --------------------- the essential local msfem problem solver class ---------------------------

template< class GlobalEntityDiscreteFunctionImp, class DiffusionOperatorImp >
class MsFEMLocalProblemSolver
{
public:
  // ! type of discrete functions
  typedef GlobalEntityDiscreteFunctionImp GlobalEntityDiscreteFunctionType;

  typedef typename GlobalEntityDiscreteFunctionType::LocalFunctionType LocalFunctionType;

  // ! type of discrete function space
  typedef typename GlobalEntityDiscreteFunctionType::DiscreteFunctionSpaceType
  GlobalEntityDiscreteFunctionSpaceType;

  // ! type of grid partition
  typedef typename GlobalEntityDiscreteFunctionSpaceType::GridPartType GridPartType;

  // ! type of grid
  typedef typename GlobalEntityDiscreteFunctionSpaceType::GridType GridType;

  // ! type of range vectors
  typedef typename GlobalEntityDiscreteFunctionSpaceType::RangeType RangeType;

  // ! type of range vectors
  typedef typename GlobalEntityDiscreteFunctionSpaceType::DomainType DomainType;

  typedef typename GlobalEntityDiscreteFunctionSpaceType::LagrangePointSetType LagrangePointSetType;

  typedef typename GlobalEntityDiscreteFunctionSpaceType::IteratorType GIteratorType;
  typedef typename GIteratorType::Entity                               GEntityType;

  enum { faceCodim = 1 };
  typedef typename LagrangePointSetType::template Codim< faceCodim >::SubEntityIteratorType FaceDofIteratorType;

  // ! polynomial order of base functions
  enum { polynomialOrder = GlobalEntityDiscreteFunctionSpaceType::polynomialOrder };

  // ! type of the (possibly non-linear) diffusion operator
  typedef DiffusionOperatorImp DiffusionType;

  struct LocProbMatrixTraits
  {
    typedef GlobalEntityDiscreteFunctionSpaceType                          RowSpaceType;
    typedef GlobalEntityDiscreteFunctionSpaceType                          ColumnSpaceType;
    typedef LagrangeMatrixSetup< false >                                   StencilType;
    typedef ParallelScalarProduct< GlobalEntityDiscreteFunctionSpaceType > ParallelScalarProductType;

    template< class M >
    struct Adapter
    {
      typedef LagrangeParallelMatrixAdapter< M > MatrixAdapterType;
    };
  };

  typedef SparseRowMatrixOperator< GlobalEntityDiscreteFunctionType, GlobalEntityDiscreteFunctionType,
                                   LocProbMatrixTraits > LocProbFEMMatrix;

  // OEMGMRESOp //OEMBICGSQOp // OEMBICGSTABOp
  typedef OEMBICGSTABOp< GlobalEntityDiscreteFunctionType, LocProbFEMMatrix > InverseLocProbFEMMatrix;

  // discrete elliptic operator describing the elliptic local msfem problems
  typedef LocalProblemOperator< GlobalEntityDiscreteFunctionType, DiffusionType > LocalProblemOperatorType;

private:
  const GlobalEntityDiscreteFunctionSpaceType& globalEntityDiscreteFunctionSpace_;   // Referenz &, wenn & verwendet,
                                                                                     // dann unten:
  DiffusionType& diffusion_;

  std::ofstream* data_file_;

public:
  // ! constructor - with diffusion operator A^{\epsilon}(x)
  MsFEMLocalProblemSolver(const GlobalEntityDiscreteFunctionSpaceType& globalEntityDiscreteFunctionSpace,
                          DiffusionType& diffusion_operator)
    : globalEntityDiscreteFunctionSpace_(globalEntityDiscreteFunctionSpace)
      , diffusion_(diffusion_operator)
      , data_file_(NULL)
  {}

  // ! constructor - with diffusion operator A^{\epsilon}(x)
  MsFEMLocalProblemSolver(const GlobalEntityDiscreteFunctionSpaceType& globalEntityDiscreteFunctionSpace,
                          DiffusionType& diffusion_operator,
                          std::ofstream& data_file)
    : globalEntityDiscreteFunctionSpace_(globalEntityDiscreteFunctionSpace)
      , diffusion_(diffusion_operator)
      , data_file_(&data_file)
  {}

  // ! ----------- method: solve the local MsFEM problem ------------------------------------------

  template< class JacobianRangeImp, class EntityPointerImp >
  void solvelocalproblem(JacobianRangeImp& gradient_PHI_H,
                         // the barycenter x_T of a macro grid element 'T'
                         EntityPointerImp& entityPointer,    // pointer on 'T'
                         GlobalEntityDiscreteFunctionType& local_problem_solution) {
    // set solution equal to zero:
    local_problem_solution.clear();

    // ! the matrix in our linear system of equations
    // in the non-linear case, it is the matrix for each iteration step
    LocProbFEMMatrix locprob_system_matrix("Local Problem System Matrix",
                                           globalEntityDiscreteFunctionSpace_,
                                           globalEntityDiscreteFunctionSpace_);

    // ! define the discrete (elliptic) local MsFEM problem operator
    // ( effect of the discretized differential operator on a certain discrete function )
    LocalProblemOperatorType local_problem_op(globalEntityDiscreteFunctionSpace_, diffusion_);

    const GridPartType& gridPart = globalEntityDiscreteFunctionSpace_.gridPart();
    typedef typename GlobalEntityDiscreteFunctionSpaceType::IteratorType GEIteratorType;
    typedef typename GridPartType::IntersectionIteratorType              GEIntersectionIteratorType;
    GEIteratorType ge_endit = globalEntityDiscreteFunctionSpace_.end();

    // ! right hand side vector of the algebraic local MsFEM problem
    // (in the non-linear setting it changes for every iteration step)
    GlobalEntityDiscreteFunctionType local_problem_rhs("rhs of local MsFEM problem", globalEntityDiscreteFunctionSpace_);
    local_problem_rhs.clear();

    // NOTE:
    // is the right hand side of the local MsFEM problem equal to zero or almost identical to zero?
    // if yes, the solution of the local MsFEM problem is also identical to zero. The solver is getting a problem with
    // this situation, which is why we do not solve local msfem problems for zero-right-hand-side, since we already know
    // the result.

    // assemble the stiffness matrix
    local_problem_op.assemble_matrix(entityPointer, locprob_system_matrix);

    // ! boundary treatment:
    typedef typename LocProbFEMMatrix::LocalMatrixType        LocalMatrix;
    typedef typename GEIntersectionIteratorType::Intersection Intersection;

    for (GEIteratorType it = globalEntityDiscreteFunctionSpace_.begin(); it != ge_endit; ++it)
    {
      const GEntityType& local_grid_entity = *it;
      if ( !local_grid_entity.hasBoundaryIntersections() )
        continue;

      LocalMatrix local_matrix = locprob_system_matrix.localMatrix(local_grid_entity, local_grid_entity);

      const LagrangePointSetType& lagrangePointSet = globalEntityDiscreteFunctionSpace_.lagrangePointSet(
        local_grid_entity);

      const GEIntersectionIteratorType iend = gridPart.iend(local_grid_entity);
      for (GEIntersectionIteratorType iit = gridPart.ibegin(local_grid_entity); iit != iend; ++iit)
      {
        const Intersection& intersection = *iit;
        if ( !intersection.boundary() )
          continue;

        const int face = intersection.indexInInside();
        const FaceDofIteratorType fdend = lagrangePointSet.template endSubEntity< 1 >(face);
        for (FaceDofIteratorType fdit = lagrangePointSet.template beginSubEntity< 1 >(face); fdit != fdend; ++fdit)
          local_matrix.unitRow(*fdit);
      }
    }
    // assemble right hand side of algebraic local msfem problem
    local_problem_op.assemble_local_RHS(entityPointer, gradient_PHI_H, local_problem_rhs);

    // zero boundary condition for 'cell problems':
    // set Dirichlet Boundary to zero
    for (GEIteratorType it = globalEntityDiscreteFunctionSpace_.begin(); it != ge_endit; ++it)
    {
      GEIntersectionIteratorType iit = gridPart.ibegin(*it);
      const GEIntersectionIteratorType ge_iendit = gridPart.iend(*it);

      for ( ; iit != ge_iendit; ++iit)
      {
        if ( !(*iit).boundary() )
          continue;

        LocalFunctionType rhs_on_entity = local_problem_rhs.localFunction(*it);

        const LagrangePointSetType& lagrangePointSet
          = globalEntityDiscreteFunctionSpace_.lagrangePointSet(*it);

        const int face = (*iit).indexInInside();

        FaceDofIteratorType faceIterator
          = lagrangePointSet.template beginSubEntity< faceCodim >(face);

        const FaceDofIteratorType faceEndIterator
          = lagrangePointSet.template endSubEntity< faceCodim >(face);

        for ( ; faceIterator != faceEndIterator; ++faceIterator)
          rhs_on_entity[*faceIterator] = 0;
      }
    }

    const double norm_rhs = local_problem_op.normRHS(local_problem_rhs);
    if ( !( local_problem_rhs.dofsValid() ) )
    {
      std::cout << "Local MsFEM Problem RHS invalid." << std::endl;
      abort();
    }

    if (norm_rhs < /*1e-06*/ 1e-10)
    {
      local_problem_solution.clear();
      // std :: cout << "Local MsFEM problem with solution zero." << std :: endl;
    } else {
      InverseLocProbFEMMatrix locprob_fem_biCGStab(locprob_system_matrix, 1e-8, 1e-8, 20000, LOCPROBLEMSOLVER_VERBOSE);
      locprob_fem_biCGStab(local_problem_rhs, local_problem_solution);
    }

    if ( !( local_problem_solution.dofsValid() ) )
    {
      std::cout << "Current solution of the local msfem problem invalid!" << std::endl;
      std::abort();
    }
  } // solvelocalproblem

  // ! ----------- end method: solve local MsFEM problem ------------------------------------------

  // method for solving and saving the solutions of the local msfem problems
  // for the whole set of macroscopic base function

  // ! ---- method: solve and save the local msfem problems for the set of macroscopic base functions -----

  // here we need a 'local msfem problem numbering manager' to determine the number of the local msfem problem
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
  // (das wird aber alles im Hintergrund vom 'local msfem problem numbering manager')

  // compute and save solutions of the local msfem problems for the base function set of the 'discreteFunctionSpace'
  // requires local msfem problem numbering manager
  template< class DiscreteFunctionImp, class LocalProblemNumberingManagerImp >
  void saveTheSolutions_baseSet(
    const typename DiscreteFunctionImp::DiscreteFunctionSpaceType& discreteFunctionSpace,
    const LocalProblemNumberingManagerImp& lp_num_manager,       // just to check, if we use the correct numeration
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

    std::string locprob_solution_location = "data/MsFEM/" + filename + "_localProblemSolutions_baseSet";
    DiscreteFunctionWriter dfw( (locprob_solution_location).c_str() );

    writer_is_open = dfw.open();

    long double starting_time = clock();

    // we want to determine minimum, average and maxiumum time for solving a local msfem problem in the current method
    double minimum_time_c_p = 1000000;
    double average_time_c_p = 0;
    double maximum_time_c_p = 0;

    int number_of_local_problem = 0;

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

        GlobalEntityDiscreteFunctionType correctorPhi_i("corrector Phi_i", globalEntityDiscreteFunctionSpace_);

        for (int i = 0; i < numBaseFunctions; ++i)
        {
          baseSet.jacobian(i, quadrature[0 /*=quadraturePoint*/], gradientPhi[i]);
          // multiply with transpose of jacobian inverse
          gradientPhi[i][0] = FMatrixHelp::mult(inv, gradientPhi[i][0]);
        }

        for (int i = 0; i < numBaseFunctions; ++i)
        {
          std::cout << "Number of the local problem: " << number_of_local_problem
                    << " (of " << numBaseFunctions* discreteFunctionSpace.grid().size(0) << " problems in total)"
                    << std::endl;

          correctorPhi_i.clear();

          // take time
          long double time_now = clock();
          gradientPhi[i][0][0] = 1.0;
          gradientPhi[i][0][1] = 0.0;

          solvelocalproblem< JacobianRangeType, IteratorType >
            (gradientPhi[i], it, correctorPhi_i);

          // min/max time
          if ( (clock() - time_now) / CLOCKS_PER_SEC > maximum_time_c_p )
          { maximum_time_c_p = (clock() - time_now) / CLOCKS_PER_SEC; }
          if ( (clock() - time_now) / CLOCKS_PER_SEC < minimum_time_c_p )
          { minimum_time_c_p = (clock() - time_now) / CLOCKS_PER_SEC; }

          dfw.append(correctorPhi_i);

          // check if we use a correct numeration of the local msfem problems:
          if ( !(lp_num_manager.get_number_of_local_problem(it, i) == number_of_local_problem) )
          {
            std::cout << "Numeration of local problems incorrect." << std::endl;
            std::abort();
          }

          #ifdef LOCALDATAOUTPUT
          // --------------- writing data output ---------------------

// if ( number_of_local_problem == 7 )
          {
            typedef Tuple< DiscreteFunctionType* >      IOTupleType;
            typedef DataOutput< GridType, IOTupleType > DataOutputType;

            // general output parameters
            LocalProblemDataOutputParameters outputparam;
            outputparam.set_path("data/MsFEM/" + filename);

            // sequence stamp
            std::stringstream outstring;

            // --------- data output local solution --------------

            // create and initialize output class
            IOTupleType local_solution_series(&correctorPhi_i);

            char ls_name[50];
            sprintf(ls_name, "a_local_problem_solution_%d", number_of_local_problem);
            std::string ls_name_s(ls_name);

            outputparam.set_prefix(ls_name_s);
            DataOutputType localsol_dataoutput(
              globalEntityDiscreteFunctionSpace_.gridPart().grid(), local_solution_series, outputparam);

            // write data
            outstring << "a-local-problem-solution";
            localsol_dataoutput.writeData( 1.0 /*dummy*/, outstring.str() );
            // clear the std::stringstream:
            outstring.str( std::string() );

            // -------------------------------------------------------

            // abort();
          }

          #endif // ifdef LOCALDATAOUTPUT

          number_of_local_problem++;
        }
      }    // end: for-loop: IteratorType it
    }   // end: 'if ( writer_is_open )'

    if (data_file_)
    {
      if ( data_file_->is_open() )
      {
        (*data_file_) << std::endl;
        (*data_file_) << "In method: saveTheSolutions_baseSet." << std::endl << std::endl;
        (*data_file_) << "Local MsFEM problems solved for " << discreteFunctionSpace.grid().size(0)
                      << " leaf entities." << std::endl;
        (*data_file_) << "Minimum time for solving a local problem = " << minimum_time_c_p << "s." << std::endl;
        (*data_file_) << "Maximum time for solving a localproblem = " << maximum_time_c_p << "s." << std::endl;
        (*data_file_) << "Average time for solving a localproblem = "
                      << ( (clock() - starting_time) / CLOCKS_PER_SEC ) / number_of_local_problem << "s." << std::endl;
        (*data_file_) << "Total time for computing and saving the localproblems = "
                      << ( (clock() - starting_time) / CLOCKS_PER_SEC ) << "s," << std::endl << std::endl;
      }
    }
  } // saveTheSolutions_baseSet
};  // end class
}

#endif // #ifndef DiscreteEllipticMsFEMLocalProblem_HH
