#ifndef DiscreteElliptic_HH
#define DiscreteElliptic_HH

#include <dune/common/fmatrix.hh>

#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/solver/oemsolver/oemsolver.hh>
#include <dune/fem/solver/inverseoperators.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/operator/discreteoperatorimp.hh>
#include <dune/fem/operator/2order/lagrangematrixsetup.hh>
#include <dune/fem/operator/matrix/spmatrix.hh>

#include <dune/multiscale/tools/assembler/righthandside_assembler.hh>

#include <dune/stuff/common/ranges.hh>

namespace Dune {
// Imp stands for Implementation
template< class DiscreteFunctionImp, class DiffusionImp, class ReactionImp >
class DiscreteEllipticOperator
  : public Operator< typename DiscreteFunctionImp::RangeFieldType, typename DiscreteFunctionImp::RangeFieldType,
                     DiscreteFunctionImp, DiscreteFunctionImp >
{
  typedef DiscreteEllipticOperator< DiscreteFunctionImp, DiffusionImp, ReactionImp > This;

public:
  typedef DiscreteFunctionImp DiscreteFunction;
  typedef DiffusionImp        DiffusionModel;
  typedef ReactionImp         Reaction;

  typedef typename DiscreteFunction::DiscreteFunctionSpaceType DiscreteFunctionSpace;

  typedef typename DiscreteFunctionSpace::GridPartType   GridPart;
  typedef typename GridPart::GridType                    Grid;
  typedef typename DiscreteFunctionSpace::RangeFieldType RangeFieldType;

  typedef typename DiscreteFunctionSpace::DomainType DomainType;
  typedef typename DiscreteFunctionSpace::RangeType  RangeType;

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
  DiscreteEllipticOperator(const DiscreteFunctionSpace& discreteFunctionSpace, const DiffusionModel& diffusion_op)
    : discreteFunctionSpace_(discreteFunctionSpace)
      , diffusion_operator_(diffusion_op)
      , reaction_coefficient_(nullptr)
  {}

  DiscreteEllipticOperator(const DiscreteFunctionSpace& discreteFunctionSpace,
                           const DiffusionModel& diffusion_op,
                           const Reaction& reaction_coefficient)
    : discreteFunctionSpace_(discreteFunctionSpace)
      , diffusion_operator_(diffusion_op)
      , reaction_coefficient_(&reaction_coefficient)
  {}

private:
  DiscreteEllipticOperator(const This&);

public:
  // dummy operator
  virtual void operator()(const DiscreteFunction& u, DiscreteFunction& w) const;

  template< class MatrixType >
  void assemble_matrix(MatrixType& global_matrix, bool boundary_treatment = true) const;

  // Matrix Assembler for local problems on a Subgrid of the Hostgrid:
  template< class MatrixType, class HostDiscreteFunctionSpaceType >
  void assemble_matrix(MatrixType& global_matrix, HostDiscreteFunctionSpaceType& hostSpace,
                       bool boundary_treatment = true) const;

  template< class MatrixType >
  void assemble_jacobian_matrix(DiscreteFunction& disc_func, MatrixType& global_matrix, bool boundary_treatment = true) const;

private:
  const DiscreteFunctionSpace& discreteFunctionSpace_;
  const DiffusionModel& diffusion_operator_;
  const Reaction* const  reaction_coefficient_;
};

// dummy implementation of "operator()"
// 'w' = effect of the discrete operator on 'u'
template< class DiscreteFunctionImp, class DiffusionImp, class ReactionImp >
void DiscreteEllipticOperator< DiscreteFunctionImp, DiffusionImp, ReactionImp >::operator()(const DiscreteFunction& /*u*/,
                                                                                            DiscreteFunction& /*w*/)
const {
  DUNE_THROW(Dune::NotImplemented,"the ()-operator of the DiscreteEllipticOperator class is not yet implemented and still a dummy.");
} // ()

template< class DiscreteFunctionImp, class DiffusionImp, class ReactionImp >
template< class MatrixType >
void DiscreteEllipticOperator< DiscreteFunctionImp, DiffusionImp, ReactionImp >::assemble_matrix(
  MatrixType& global_matrix,
  bool boundary_treatment ) const {
  typedef typename MatrixType::LocalMatrixType LocalMatrix;

  global_matrix.reserve();
  global_matrix.clear();

  std::vector< typename BaseFunctionSet::JacobianRangeType > gradient_phi( discreteFunctionSpace_.mapper().maxNumDofs() );

  // micro scale base function:
  std::vector< RangeType > phi( discreteFunctionSpace_.mapper().maxNumDofs() );

  const Iterator end = discreteFunctionSpace_.end();
  for (Iterator it = discreteFunctionSpace_.begin(); it != end; ++it)
  {
    const Entity& entity = *it;
    const Geometry& geometry = entity.geometry();
    assert(entity.partitionType() == InteriorEntity);

    LocalMatrix local_matrix = global_matrix.localMatrix(entity, entity);

    const BaseFunctionSet& baseSet = local_matrix.domainBaseFunctionSet();
    const unsigned int numBaseFunctions = baseSet.size();

    // for constant diffusion "2*discreteFunctionSpace_.order()" is sufficient, for the general case, it is better to
    // use a higher order quadrature:
    const Quadrature quadrature(entity, 2 * discreteFunctionSpace_.order() + 2);
    const size_t numQuadraturePoints = quadrature.nop();
    for (size_t quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint)
    {
      // local (barycentric) coordinates (with respect to entity)
      const typename Quadrature::CoordinateType& local_point = quadrature.point(quadraturePoint);

      const DomainType global_point = geometry.global(local_point);

      const double weight = quadrature.weight(quadraturePoint) * geometry.integrationElement(local_point);

      // transposed of the the inverse jacobian
      const FieldMatrix< double, dimension, dimension >& inverse_jac
        = geometry.jacobianInverseTransposed(local_point);

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
        // A( \nabla \phi ) // diffusion operator evaluated in (x,\nabla \phi)
        typename LocalFunction::JacobianRangeType diffusion_in_gradient_phi;
        diffusion_operator_.diffusiveFlux(global_point, gradient_phi[i], diffusion_in_gradient_phi);
        for (unsigned int j = 0; j < numBaseFunctions; ++j)
        {
          local_matrix.add( j, i, weight * (diffusion_in_gradient_phi[0] * gradient_phi[j][0]) );

          if (reaction_coefficient_)
          {
            RangeType c;
            reaction_coefficient_->evaluate(global_point, c);
            local_matrix.add( j, i, weight * c * (phi[i][0] * phi[j][0]) );
          }
        }
      }
    }
  }

  // boundary treatment
  if (boundary_treatment)
  {
    const GridPart& gridPart = discreteFunctionSpace_.gridPart();
    for (Iterator it = discreteFunctionSpace_.begin(); it != end; ++it)
    {
      const Entity& entity = *it;
      if ( !entity.hasBoundaryIntersections() )
        continue;

      LocalMatrix local_matrix = global_matrix.localMatrix(entity, entity);

      const LagrangePointSet& lagrangePointSet = discreteFunctionSpace_.lagrangePointSet(entity);

      const IntersectionIterator iend = gridPart.iend(entity);
      for (IntersectionIterator iit = gridPart.ibegin(entity); iit != iend; ++iit)
      {
        const Intersection& intersection = *iit;
        if ( !intersection.boundary() )
          continue;

        const int face = intersection.indexInInside();
        const FaceDofIterator fdend = lagrangePointSet.template endSubEntity< 1 >(face);
        for (FaceDofIterator fdit = lagrangePointSet.template beginSubEntity< 1 >(face); fdit != fdend; ++fdit)
          local_matrix.unitRow(*fdit);
      }
    }
  }
} // assemble_matrix

// Matrix Assembler for local problems on a Subgrid of the Hostgrid:
template< class DiscreteFunctionImp, class DiffusionImp, class ReactionImp >
template< class MatrixType, class HostDiscreteFunctionSpaceType >
void DiscreteEllipticOperator< DiscreteFunctionImp,
                               DiffusionImp,
                               ReactionImp >::assemble_matrix
  (MatrixType& global_matrix, HostDiscreteFunctionSpaceType& hostSpace, bool boundary_treatment ) const {
  typedef typename MatrixType::LocalMatrixType LocalMatrix;

  global_matrix.reserve();
  global_matrix.clear();

  std::vector< typename BaseFunctionSet::JacobianRangeType > gradient_phi( discreteFunctionSpace_.mapper().maxNumDofs() );

  // micro scale base function:
  std::vector< RangeType > phi( discreteFunctionSpace_.mapper().maxNumDofs() );

  const Iterator end = discreteFunctionSpace_.end();
  for (Iterator it = discreteFunctionSpace_.begin(); it != end; ++it)
  {
    const Entity& entity = *it;
    const Geometry& geometry = entity.geometry();
    assert(entity.partitionType() == InteriorEntity);

    LocalMatrix local_matrix = global_matrix.localMatrix(entity, entity);

    const BaseFunctionSet& baseSet = local_matrix.domainBaseFunctionSet();
    const unsigned int numBaseFunctions = baseSet.numBaseFunctions();

    // for constant diffusion "2*discreteFunctionSpace_.order()" is sufficient, for the general case, it is better to
    // use a higher order quadrature:
    const Quadrature quadrature(entity, 2 * discreteFunctionSpace_.order() + 2);
    const size_t numQuadraturePoints = quadrature.nop();
    for (size_t quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint)
    {
      // local (barycentric) coordinates (with respect to entity)
      const typename Quadrature::CoordinateType& local_point = quadrature.point(quadraturePoint);

      const DomainType global_point = geometry.global(local_point);

      const double weight = quadrature.weight(quadraturePoint) * geometry.integrationElement(local_point);

      // transposed of the the inverse jacobian
      const FieldMatrix< double, dimension, dimension >& inverse_jac
        = geometry.jacobianInverseTransposed(local_point);

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
        // A( \nabla \phi ) // diffusion operator evaluated in (x,\nabla \phi)
        typename LocalFunction::JacobianRangeType diffusion_in_gradient_phi;
        diffusion_operator_.diffusiveFlux(global_point, gradient_phi[i], diffusion_in_gradient_phi);
        for (unsigned int j = 0; j < numBaseFunctions; ++j)
        {
          local_matrix.add( j, i, weight * (diffusion_in_gradient_phi[0] * gradient_phi[j][0]) );

          if (reaction_coefficient_)
          {
            RangeType c;
            reaction_coefficient_->evaluate(global_point, c);
            local_matrix.add( j, i, weight * c * (phi[i][0] * phi[j][0]) );
          }
        }
      }
    }
  }

  // boundary treatment
  if (boundary_treatment)
  {
    typedef typename HostDiscreteFunctionSpaceType::GridPartType HostGridPartType;
    typedef typename HostDiscreteFunctionSpaceType::GridType     HostGridType;

    typedef typename HostDiscreteFunctionSpaceType::IteratorType::Entity HostEntityType;
    typedef typename HostEntityType::EntityPointer                       HostEntityPointerType;

    typedef typename HostGridPartType::IntersectionIteratorType HostIntersectionIterator;

    const HostGridPartType& hostGridPart = hostSpace.gridPart();

    const GridPart& gridPart = discreteFunctionSpace_.gridPart();
    const Grid& subGrid = discreteFunctionSpace_.grid();

    for (Iterator it = discreteFunctionSpace_.begin(); it != end; ++it)
    {
      const Entity& entity = *it;

      const HostEntityPointerType host_entity_pointer = subGrid.template getHostEntity< 0 >(entity);
      const HostEntityType& host_entity = *host_entity_pointer;

      LocalMatrix local_matrix = global_matrix.localMatrix(entity, entity);

      const LagrangePointSet& lagrangePointSet = discreteFunctionSpace_.lagrangePointSet(entity);

      const HostIntersectionIterator iend = hostGridPart.iend(host_entity);
      for (HostIntersectionIterator iit = hostGridPart.ibegin(host_entity); iit != iend; ++iit)
      {
        if ( iit->neighbor() ) // if there is a neighbor entity
        {
          // check if the neighbor entity is in the subgrid
          const HostEntityPointerType neighborHostEntityPointer = iit->outside();
          const HostEntityType& neighborHostEntity = *neighborHostEntityPointer;
          if ( subGrid.template contains< 0 >(neighborHostEntity) )
          {
            continue;
          }
        }

        const int face = (*iit).indexInInside();
        const FaceDofIterator fdend = lagrangePointSet.template endSubEntity< 1 >(face);
        for (FaceDofIterator fdit = lagrangePointSet.template beginSubEntity< 1 >(face); fdit != fdend; ++fdit)
          local_matrix.unitRow(*fdit);
      }
    }
  }
} // assemble_matrix

// assemble stiffness matrix for the jacobian matrix of the diffusion operator evaluated in the gradient of a certain
// discrete function (in case of the Newton method, it is the preceeding iterate u_H^{(n-1)} )
// stiffness matrix with entries
// \int JA(\nabla disc_func) \nabla phi_i \nabla phi_j
// (here, JA denotes the jacobian matrix of the diffusion operator A)
template< class DiscreteFunctionImp, class DiffusionImp, class ReactionImp >
template< class MatrixType >
void DiscreteEllipticOperator< DiscreteFunctionImp, DiffusionImp, ReactionImp >::assemble_jacobian_matrix(
  DiscreteFunction& disc_func,
  MatrixType& global_matrix,
  bool boundary_treatment ) const {
  typedef typename MatrixType::LocalMatrixType LocalMatrix;

  typedef typename DiscreteFunction::LocalFunctionType
  LocalFunction;

  global_matrix.reserve();
  global_matrix.clear();

  std::vector< typename BaseFunctionSet::JacobianRangeType > gradient_phi( discreteFunctionSpace_.mapper().maxNumDofs() );

  // micro scale base function:
  std::vector< RangeType > phi( discreteFunctionSpace_.mapper().maxNumDofs() );

  const Iterator end = discreteFunctionSpace_.end();
  for (Iterator it = discreteFunctionSpace_.begin(); it != end; ++it)
  {
    const Entity& entity = *it;
    const Geometry& geometry = entity.geometry();
    assert(entity.partitionType() == InteriorEntity);

    LocalMatrix local_matrix = global_matrix.localMatrix(entity, entity);
    LocalFunction local_disc_function = disc_func.localFunction(entity);

    const BaseFunctionSet& baseSet = local_matrix.domainBaseFunctionSet();
    const unsigned int numBaseFunctions = baseSet.size();

    // for constant diffusion "2*discreteFunctionSpace_.order()" is sufficient, for the general case, it is better to
    // use a higher order quadrature:
    const Quadrature quadrature(entity, 2 * discreteFunctionSpace_.order() + 2);
    const size_t numQuadraturePoints = quadrature.nop();
    for (size_t quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint)
    {
      // local (barycentric) coordinates (with respect to entity)
      const typename Quadrature::CoordinateType& local_point = quadrature.point(quadraturePoint);

      const DomainType global_point = geometry.global(local_point);

      const double weight = quadrature.weight(quadraturePoint) * geometry.integrationElement(local_point);

      // transposed of the the inverse jacobian
      const FieldMatrix< double, dimension, dimension >& inverse_jac
        = geometry.jacobianInverseTransposed(local_point);

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
        typename BaseFunctionSet::JacobianRangeType grad_local_disc_function;
        local_disc_function.jacobian(quadrature[quadraturePoint], grad_local_disc_function);
        // here: no multiplication with jacobian inverse transposed required!

        // JA( \nabla u_H ) \nabla phi_i // jacobian of diffusion operator evaluated in (x,grad_local_disc_function) in
        // direction of the gradient of the current base function
        typename LocalFunction::JacobianRangeType jac_diffusion_flux;
        diffusion_operator_.jacobianDiffusiveFlux(global_point,
                                                  grad_local_disc_function,
                                                  gradient_phi[i],
                                                  jac_diffusion_flux);

        for (unsigned int j = 0; j < numBaseFunctions; ++j)
        {
          local_matrix.add( j, i, weight * (jac_diffusion_flux[0] * gradient_phi[j][0]) );

          if (reaction_coefficient_)
          {
            RangeType c;
            reaction_coefficient_->evaluate(global_point, c);
            local_matrix.add( j, i, weight * c * (phi[i][0] * phi[j][0]) );
          }
        }
      }
    }
  }

  // boundary treatment
  if (boundary_treatment)
  {
    const GridPart& gridPart = discreteFunctionSpace_.gridPart();
    for (Iterator it = discreteFunctionSpace_.begin(); it != end; ++it)
    {
      const Entity& entity = *it;
      if ( !entity.hasBoundaryIntersections() )
        continue;

      LocalMatrix local_matrix = global_matrix.localMatrix(entity, entity);

      const LagrangePointSet& lagrangePointSet = discreteFunctionSpace_.lagrangePointSet(entity);

      const IntersectionIterator iend = gridPart.iend(entity);
      for (IntersectionIterator iit = gridPart.ibegin(entity); iit != iend; ++iit)
      {
        const Intersection& intersection = *iit;
        if ( !intersection.boundary() )
          continue;

        const int face = intersection.indexInInside();
        const FaceDofIterator fdend = lagrangePointSet.template endSubEntity< 1 >(face);
        for (FaceDofIterator fdit = lagrangePointSet.template beginSubEntity< 1 >(face); fdit != fdend; ++fdit)
          local_matrix.unitRow(*fdit);
      }
    }
  }
} // assemble_jacobian_matrix
}

#endif // #ifndef DiscreteElliptic_HH
