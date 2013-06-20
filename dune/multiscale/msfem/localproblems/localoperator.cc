// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifdef HAVE_CMAKE_CONFIG
 #include "cmake_config.h"
#elif defined (HAVE_CONFIG_H)
 #include <config.h>
#endif // ifdef HAVE_CMAKE_CONFIG

#include "localoperator.hh"

#include <dune/stuff/common/ranges.hh>
#include <dune/stuff/fem/localmatrix_proxy.hh>

#include <dune/multiscale/tools/misc/uzawa.hh>
#include <dune/multiscale/tools/misc/weighted-clement-operator.hh>
#include <dune/multiscale/problems/elliptic/selector.hh>

namespace Dune {
namespace Multiscale {
namespace MsFEM {

LocalProblemOperator::LocalProblemOperator(const DiscreteFunctionSpace& subDiscreteFunctionSpace,
                                           const DiffusionModel& diffusion_op)
  : subDiscreteFunctionSpace_(subDiscreteFunctionSpace)
    , diffusion_operator_(diffusion_op)
{}

// is a given 'point' in the convex hull of corner 0, corner 1 and corner 2 (which forms a codim 0 entity)
bool LocalProblemOperator::point_is_in_element( const DomainType& corner_0,
                const DomainType& corner_1,
                const DomainType& corner_2,
                const DomainType& point) const
    {
      DomainType v = corner_0 - corner_2;
      DomainType w = corner_1 - corner_2;
      DomainType p = point - corner_2;

      double lambda_0, lambda_1;

      if ( v[0] != 0.0 )
      {
    double val0 = p[1] - ( (v[1]/v[0])*p[0]);
    double val1 = w[1] - ( (v[1]/v[0])*w[0]);

    if ( val1 != 0.0 )
    {
      lambda_1 = val0 / val1;
      lambda_0 = (p[0] - (lambda_1*w[0]))/v[0];
          if ( (0.0 <= lambda_0) && (1.0 >= lambda_0) && (0.0 <= lambda_1) && (1.0 >= lambda_1) && ( lambda_0+lambda_1<= 1.0) )
        return true;
      else
        return false;
    }
    else
    {
      DUNE_THROW(Dune::InvalidStateException, "... in method 'point_is_in_element': Given corners do not span a codim 0 entity in 2D.");
    }
      }
      else
      {
    if ( (w[0] != 0.0) && (v[1] != 0.0) )
    {
      lambda_1 = p[0] / w[0];
      lambda_0 = (p[1] - (lambda_1*w[1]))/v[1];
          if ( (0.0 <= lambda_0) && (1.0 >= lambda_0) && (0.0 <= lambda_1) && (1.0 >= lambda_1) && ( lambda_0+lambda_1<= 1.0) )
        return true;
      else
        return false;
    }
    else
    {DUNE_THROW(Dune::InvalidStateException, "... in method 'point_is_in_element': Given corners do not span a codim 0 entity in 2D.");}
      }

    }



//! stiffness matrix for a linear elliptic diffusion operator
// for oversampling strategy 1 (no constraints)
void LocalProblemOperator::assemble_matrix(MsFEMLocalProblemSolver::LocProbFEMMatrix& global_matrix) const
// x_T is the barycenter of the macro grid element T
{
  global_matrix.reserve();
  global_matrix.clear();

  // local grid basis functions:
  std::vector< RangeType > phi( subDiscreteFunctionSpace_.mapper().maxNumDofs() );

  // gradient of micro scale base function:
  std::vector< typename BaseFunctionSet::JacobianRangeType > gradient_phi(
    subDiscreteFunctionSpace_.mapper().maxNumDofs() );
  typename LocalFunction::JacobianRangeType diffusion_in_gradient_phi;

  const Iterator end = subDiscreteFunctionSpace_.end();
  for (Iterator it = subDiscreteFunctionSpace_.begin(); it != end; ++it)
  {
    const Entity& sub_grid_entity = *it;
    const Geometry& sub_grid_geometry = sub_grid_entity.geometry();

    DSFe::LocalMatrixProxy<MsFEMLocalProblemSolver::LocProbFEMMatrix> local_matrix(global_matrix, sub_grid_entity, sub_grid_entity);

    const BaseFunctionSet& baseSet = local_matrix.domainBaseFunctionSet();
    const auto numBaseFunctions = baseSet.size();

    // for constant diffusion "2*discreteFunctionSpace_.order()" is sufficient, for the general case, it is better to
    // use a higher order quadrature:
    const Quadrature quadrature(sub_grid_entity, 2 * subDiscreteFunctionSpace_.order() + 2);
    const size_t numQuadraturePoints = quadrature.nop();
    for (size_t quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint)
    {
      // local (barycentric) coordinates (with respect to local grid entity)
      const typename Quadrature::CoordinateType& local_point = quadrature.point(quadraturePoint);
      const DomainType global_point = sub_grid_geometry.global(local_point);

      const double weight = quadrature.weight(quadraturePoint)
                            * sub_grid_geometry.integrationElement(local_point);

      // transposed of the the inverse jacobian
      const auto& inverse_jac = sub_grid_geometry.jacobianInverseTransposed(local_point);
      baseSet.jacobianAll(quadrature[quadraturePoint], inverse_jac, gradient_phi);
      baseSet.evaluateAll(quadrature[quadraturePoint], phi);

      for (unsigned int i = 0; i < numBaseFunctions; ++i)
      {
        // A( x, \nabla \phi(x) )
        diffusion_operator_.diffusiveFlux(global_point, gradient_phi[i], diffusion_in_gradient_phi);
        for (unsigned int j = 0; j < numBaseFunctions; ++j)
        {
          // stiffness contribution
          local_matrix.add( j, i, weight * (diffusion_in_gradient_phi[0] * gradient_phi[j][0]) );
          // mass contribution (just for stabilization!)
          // local_matrix.add( j, i, 0.00000001 * weight * (phi[ i ][ 0 ] * phi[ j ][ 0 ]) );
        }
      }
    }
  }
} // assemble_matrix


//! stiffness matrix for a linear elliptic diffusion operator
void LocalProblemOperator::assemble_matrix(MsFEMLocalProblemSolver::LocProbFEMMatrix& global_matrix,
                                           const SubGridList::CoarseNodeVectorType& coarse_node_vector ) const
// x_T is the barycenter of the macro grid element T
{
  global_matrix.reserve();
  global_matrix.clear();

  // local grid basis functions:
  std::vector< RangeType > phi( subDiscreteFunctionSpace_.mapper().maxNumDofs() );

  // gradient of micro scale base function:
  std::vector< typename BaseFunctionSet::JacobianRangeType > gradient_phi(
    subDiscreteFunctionSpace_.mapper().maxNumDofs() );

  const Iterator end = subDiscreteFunctionSpace_.end();
  for (Iterator it = subDiscreteFunctionSpace_.begin(); it != end; ++it)
  {
    const Entity& sub_grid_entity = *it;
    const Geometry& sub_grid_geometry = sub_grid_entity.geometry();
    assert(sub_grid_entity.partitionType() == InteriorEntity);

    std::vector< int > sub_grid_entity_corner_is_relevant;
    for ( int c = 0; c < sub_grid_geometry.corners(); ++c )
      {
        for ( size_t coarse_node_local_id = 0; coarse_node_local_id < coarse_node_vector.size(); ++coarse_node_local_id )
        {
     // if the subgrid corner 'c' is in the 'relevant coarse node vector' and if 'c' was not yet added to the
     // vector 'sub_grid_entity_corner_is_relevant' then add it to the vector
         if ( (coarse_node_vector[coarse_node_local_id] == sub_grid_geometry.corner(c))
         && (std::find(sub_grid_entity_corner_is_relevant.begin(), sub_grid_entity_corner_is_relevant.end(), c) == sub_grid_entity_corner_is_relevant.end()) )
         { sub_grid_entity_corner_is_relevant.push_back(c); }
        }
      }

    DSFe::LocalMatrixProxy<MsFEMLocalProblemSolver::LocProbFEMMatrix> local_matrix(global_matrix, sub_grid_entity, sub_grid_entity);

    const BaseFunctionSet& baseSet = local_matrix.domainBaseFunctionSet();
    const auto numBaseFunctions = baseSet.size();

    std::vector<RangeType> value_phi(numBaseFunctions);
    // for constant diffusion "2*discreteFunctionSpace_.order()" is sufficient, for the general case, it is better to
    // use a higher order quadrature:
    const Quadrature quadrature(sub_grid_entity, 2 * subDiscreteFunctionSpace_.order() + 2);
    const size_t numQuadraturePoints = quadrature.nop();
    for (size_t quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint)
    {
      // local (barycentric) coordinates (with respect to local grid entity)
      const typename Quadrature::CoordinateType& local_point = quadrature.point(quadraturePoint);
      const DomainType global_point = sub_grid_geometry.global(local_point);

      const double weight = quadrature.weight(quadraturePoint)
                            * sub_grid_geometry.integrationElement(local_point);

      // transposed of the the inverse jacobian
      const auto& inverse_jac = sub_grid_geometry.jacobianInverseTransposed(local_point);
      baseSet.jacobianAll(quadrature[quadraturePoint], inverse_jac, gradient_phi);
      baseSet.evaluateAll(quadrature[quadraturePoint], phi);

      for ( size_t sgec = 0; sgec < sub_grid_entity_corner_is_relevant.size(); ++sgec )
      {
        baseSet.evaluateAll(sub_grid_geometry.local(sub_grid_geometry.corner(sub_grid_entity_corner_is_relevant[sgec])), value_phi);
        for (unsigned int i = 0; i < numBaseFunctions; ++i)
        {
          if ( value_phi[i] == 1.0 )
          {
            assert( dimension == 2);
            phi[i][0] = 0.0;
            gradient_phi[i][0][0] = 0.0;
            gradient_phi[i][0][1] = 0.0;
          }
        }
      }

      for (unsigned int i = 0; i < numBaseFunctions; ++i)
      {
        // A( x, \nabla \phi(x) )
        typename LocalFunction::JacobianRangeType diffusion_in_gradient_phi;
        diffusion_operator_.diffusiveFlux(global_point, gradient_phi[i], diffusion_in_gradient_phi);
        for (unsigned int j = 0; j < numBaseFunctions; ++j)
        {
          // stiffness contribution
          local_matrix.add( j, i, weight * (diffusion_in_gradient_phi[0] * gradient_phi[j][0]) );
          // mass contribution (just for stabilization!)
          // local_matrix.add( j, i, 0.00000001 * weight * (phi[ i ][ 0 ] * phi[ j ][ 0 ]) );
        }
      }
    }
  }
} // assemble_matrix

void LocalProblemOperator::printLocalRHS(const LocalProblemOperator::DiscreteFunction& rhs) const {
  typedef typename DiscreteFunction::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionSpaceType::IteratorType        IteratorType;
  typedef typename DiscreteFunction::LocalFunctionType         LocalFunctionType;

  const DiscreteFunctionSpaceType& discreteFunctionSpace
    = rhs.space();

  const IteratorType endit = discreteFunctionSpace.end();
  for (IteratorType it = discreteFunctionSpace.begin(); it != endit; ++it)
  {
    LocalFunctionType elementOfRHS = rhs.localFunction(*it);

    const int numDofs = elementOfRHS.numDofs();
    for (int i = 0; i < numDofs; ++i)
    {
      DSC_LOG_DEBUG << "Number of Dof: " << i << " ; " << rhs.name() << " : " << elementOfRHS[i] << std::endl;
    }
  }
}  // end method


void LocalProblemOperator::set_zero_boundary_condition_RHS(const HostDiscreteFunctionSpace& host_space,
                                                           LocalProblemOperator::DiscreteFunction& rhs ) const {
  typedef typename DiscreteFunction::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionSpaceType::IteratorType        IteratorType;
  typedef typename DiscreteFunction::LocalFunctionType         LocalFunctionType;

  const DiscreteFunctionSpaceType& discreteFunctionSpace
    = rhs.space();

  const GridType& subGrid = discreteFunctionSpace.grid();
  const HostGridPart& hostGridPart = host_space.gridPart();
  
  // set Dirichlet Boundary to zero
  const Iterator endit = discreteFunctionSpace.end();
  for (Iterator it = discreteFunctionSpace.begin(); it != endit; ++it)
  {
    const Entity& subgrid_entity = *it;

    HostEntityPointer host_entity_pointer = subGrid.getHostEntity< 0 >(subgrid_entity);
    const HostEntity& host_entity = *host_entity_pointer;

    HostIntersectionIterator iit = hostGridPart.ibegin(host_entity);
    const HostIntersectionIterator endiit = hostGridPart.iend(host_entity);
    for ( ; iit != endiit; ++iit)
    {
      if ( iit->neighbor() ) // if there is a neighbor entity
      {
        // check if the neighbor entity is in the subgrid
        const HostEntityPointer neighborHostEntityPointer = iit->outside();
        const HostEntity& neighborHostEntity = *neighborHostEntityPointer;

        if ( subGrid.contains< 0 >(neighborHostEntity) )
        {
          continue;
        }
      }

      const LagrangePointSet& lagrangePointSet
          = discreteFunctionSpace.lagrangePointSet(subgrid_entity);

      const int face = (*iit).indexInInside();

      FaceDofIterator faceIterator
          = lagrangePointSet.beginSubEntity< faceCodim >(face);
      const FaceDofIterator faceEndIterator
          = lagrangePointSet.endSubEntity< faceCodim >(face);

      for ( ; faceIterator != faceEndIterator; ++faceIterator)
        (rhs.localFunction( subgrid_entity ))[*faceIterator] = 0;
    }

  }

}  // end method



double LocalProblemOperator::normRHS(const LocalProblemOperator::DiscreteFunction& rhs) const {
  double norm = 0.0;

  typedef typename DiscreteFunction::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionSpaceType::IteratorType        IteratorType;
  typedef typename IteratorType::Entity                           EntityType;
  typedef typename DiscreteFunction::LocalFunctionType         LocalFunctionType;
  typedef typename DiscreteFunctionSpaceType::GridPartType        GridPartType;
  typedef typename DiscreteFunctionSpaceType::GridType            GridType;
  typedef typename GridType::Codim< 0 >::Geometry
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


void LocalProblemOperator
      ::assemble_local_RHS(const JacobianRangeType &e, // direction 'e'
                            // rhs local msfem problem:
                            LocalProblemOperator::DiscreteFunction& local_problem_RHS) const {
  typedef typename DiscreteFunction::DiscreteFunctionSpaceType DiscreteFunctionSpace;
  typedef typename DiscreteFunction::LocalFunctionType         LocalFunction;

  typedef typename DiscreteFunctionSpace::BaseFunctionSetType BaseFunctionSet;
  typedef typename DiscreteFunctionSpace::IteratorType        Iterator;
  typedef typename Iterator::Entity                           Entity;
  typedef typename Entity::Geometry                           Geometry;

  typedef typename DiscreteFunctionSpace::GridPartType GridPart;
  typedef CachingQuadrature< GridPart, 0 >             Quadrature;

  const DiscreteFunctionSpace& discreteFunctionSpace = local_problem_RHS.space();

  // set entries to zero:
  local_problem_RHS.clear();

  // gradient of micro scale base function:
  std::vector< JacobianRangeType > gradient_phi( discreteFunctionSpace.mapper().maxNumDofs() );

  const Iterator end = discreteFunctionSpace.end();
  for (Iterator it = discreteFunctionSpace.begin(); it != end; ++it)
  {
    const Entity& local_grid_entity = *it;
    const Geometry& geometry = local_grid_entity.geometry();

    LocalFunction elementOfRHS = local_problem_RHS.localFunction(local_grid_entity);

    const BaseFunctionSet& baseSet = elementOfRHS.baseFunctionSet();
    const auto numBaseFunctions = baseSet.size();

    const Quadrature quadrature(local_grid_entity, 2 * discreteFunctionSpace.order() + 2);
    const size_t numQuadraturePoints = quadrature.nop();
    for (size_t quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint)
    {
      const typename Quadrature::CoordinateType& local_point = quadrature.point(quadraturePoint);

      // remember, we are concerned with: - \int_{U(T)} (A^eps)(x) e · ∇ \phi(x)

      // global point in the subgrid
      const DomainType global_point = geometry.global(local_point);

      const double weight = quadrature.weight(quadraturePoint) * geometry.integrationElement(local_point);

      // transposed of the the inverse jacobian
      const FieldMatrix< double, dimension, dimension >& inverse_jac
        = geometry.jacobianInverseTransposed(local_point);

      // A^eps(x) e
      // diffusion operator evaluated in 'x' multiplied with e
      JacobianRangeType diffusion_in_e;
      diffusion_operator_.diffusiveFlux(global_point, e, diffusion_in_e);

      baseSet.jacobianAll(quadrature[quadraturePoint], inverse_jac, gradient_phi);
      for (unsigned int i = 0; i < numBaseFunctions; ++i)
      {
        elementOfRHS[i] -= weight * (diffusion_in_e[0] * gradient_phi[i][0]);
      }
    }
  }
} // assemble_local_RHS


void LocalProblemOperator
      ::assemble_local_RHS(const JacobianRangeType &e, // direction 'e'
                           const SubGridList::CoarseNodeVectorType& coarse_node_vector, // for constraints on the space
                           const int& oversampling_strategy,
                           // rhs local msfem problem:
                           DiscreteFunction& local_problem_RHS) const {
  typedef typename DiscreteFunction::DiscreteFunctionSpaceType DiscreteFunctionSpace;
  typedef typename DiscreteFunction::LocalFunctionType         LocalFunction;

  typedef typename DiscreteFunctionSpace::BaseFunctionSetType BaseFunctionSet;
  typedef typename DiscreteFunctionSpace::IteratorType        Iterator;
  typedef typename Iterator::Entity                           Entity;
  typedef typename Entity::Geometry                           Geometry;

  typedef typename DiscreteFunctionSpace::GridPartType GridPart;
  typedef CachingQuadrature< GridPart, 0 >             Quadrature;

  const DiscreteFunctionSpace& discreteFunctionSpace = local_problem_RHS.space();

  // set entries to zero:
  local_problem_RHS.clear();

  // gradient of micro scale base function:
  std::vector< JacobianRangeType > gradient_phi( discreteFunctionSpace.mapper().maxNumDofs() );

  const Iterator end = discreteFunctionSpace.end();
  for (Iterator it = discreteFunctionSpace.begin(); it != end; ++it)
  {
    const Entity& local_grid_entity = *it;
    const Geometry& geometry = local_grid_entity.geometry();
    assert(local_grid_entity.partitionType() == InteriorEntity);

    // for strategy 3, we only integrate over 'T' instead of 'U(T)', therefor check if 'it' belongs to 'T':
    if ( oversampling_strategy == 3 )
      {
        // the first three elements of the 'coarse_node_vector' are the corners of the relevant coarse grid entity
        // (the coarse grid entity that was the starting entity to create the current subgrid that was constructed by enrichment)
        if ( !(point_is_in_element( coarse_node_vector[0], coarse_node_vector[1], coarse_node_vector[2], geometry.center() )) )
       continue;
      }

    // 'oversampling_strategy == 3' means that we use the rigorous MsFEM
    bool clement = false;
    if (oversampling_strategy == 3)
     clement = (DSC_CONFIG_GET( "rigorous_msfem.oversampling_strategy", "Clement" ) == "Clement" );

    // if the we use oversampling strategy 2 or 3/Lagrange, we need to sort out some coarse grid nodes:
    std::vector< int > sub_grid_entity_corner_is_relevant;
    if (!clement)
    {
      for ( int c = 0; c < geometry.corners(); ++c )
      {
        for ( size_t coarse_node_local_id = 0; coarse_node_local_id < coarse_node_vector.size(); ++coarse_node_local_id )
        {
         // if the subgrid corner 'c' is in the 'relevant coarse node vector' and if 'c' was not yet added to the
         // vector 'sub_grid_entity_corner_is_relevant' then add it to the vector
         if ( (coarse_node_vector[coarse_node_local_id] == geometry.corner(c))
         && (std::find(sub_grid_entity_corner_is_relevant.begin(), sub_grid_entity_corner_is_relevant.end(), c) == sub_grid_entity_corner_is_relevant.end()) )
         { sub_grid_entity_corner_is_relevant.push_back(c);
         //std :: cout << std ::endl << "geometry.corner(" << c << ") = " << geometry.corner(c) << " is relevant." << std ::endl;

        }
        }
      }
    }
    LocalFunction elementOfRHS = local_problem_RHS.localFunction(local_grid_entity);

    const BaseFunctionSet& baseSet = elementOfRHS.baseFunctionSet();
    const auto numBaseFunctions = baseSet.size();

    const Quadrature quadrature(local_grid_entity, 2 * discreteFunctionSpace.order() + 2);
    const size_t numQuadraturePoints = quadrature.nop();
    for (size_t quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint)
    {
      const typename Quadrature::CoordinateType& local_point = quadrature.point(quadraturePoint);

      // remember, we are concerned with: - \int_{U(T)} (A^eps)(x) e · ∇ \phi(x)

      // global point in the subgrid
      const DomainType global_point = geometry.global(local_point);

      const double weight = quadrature.weight(quadraturePoint) * geometry.integrationElement(local_point);

      // transposed of the the inverse jacobian
      const auto& inverse_jac = geometry.jacobianInverseTransposed(local_point);

      // A^eps(x) e
      // diffusion operator evaluated in 'x' multiplied with e
      JacobianRangeType diffusion_in_e;
      diffusion_operator_.diffusiveFlux(global_point, e, diffusion_in_e);
      baseSet.jacobianAll(quadrature[quadraturePoint], inverse_jac, gradient_phi);
      std::vector< std::vector<RangeType> > phi_values(sub_grid_entity_corner_is_relevant.size());
      for(auto j : DSC::valueRange(sub_grid_entity_corner_is_relevant.size())) {
        baseSet.evaluateAll(geometry.local(geometry.corner(sub_grid_entity_corner_is_relevant[j])), phi_values[j]);
      }

      for (unsigned int i = 0; i < numBaseFunctions; ++i)
      {
        bool zero_entry = false;
        for ( size_t sgec = 0; sgec < sub_grid_entity_corner_is_relevant.size(); ++sgec )
        {
          const auto& value_phi_i = phi_values[sgec][i];
          if ( value_phi_i == 1.0 )
          {
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


void LocalProblemOperator
      ::assemble_local_RHS_lg_problems( const HostDiscreteFunction& coarse_basis_func, double clement_weight,
                                           DiscreteFunction& local_problem_RHS ) const {


  typedef typename DiscreteFunction::DiscreteFunctionSpaceType DiscreteFunctionSpace;
  typedef typename DiscreteFunction::LocalFunctionType         LocalFunction;

  typedef typename DiscreteFunctionSpace::BaseFunctionSetType BaseFunctionSet;
  typedef typename DiscreteFunctionSpace::IteratorType        Iterator;
  typedef typename Iterator::Entity                           Entity;
  typedef typename Entity::Geometry                           Geometry;

  typedef typename DiscreteFunctionSpace::GridPartType GridPart;
  typedef CachingQuadrature< GridPart, 0 >             Quadrature;

  const DiscreteFunctionSpace& discreteFunctionSpace = local_problem_RHS.space();

  const GridType& subGrid = discreteFunctionSpace.grid();
  
  // set entries to zero:
  local_problem_RHS.clear();

  // gradient of micro scale base function:
  std::vector< JacobianRangeType > gradient_phi( discreteFunctionSpace.mapper().maxNumDofs() );

  const Iterator end = discreteFunctionSpace.end();
  for (Iterator it = discreteFunctionSpace.begin(); it != end; ++it)
  {
    const Entity& local_grid_entity = *it;
    const Geometry& geometry = local_grid_entity.geometry();
    assert(local_grid_entity.partitionType() == InteriorEntity);

    LocalFunction elementOfRHS = local_problem_RHS.localFunction(local_grid_entity);

    const BaseFunctionSet& baseSet = elementOfRHS.baseFunctionSet();
    const auto numBaseFunctions = baseSet.size();
        
    HostEntityPointer host_entity_pointer = subGrid.getHostEntity< 0 >( local_grid_entity );
    const HostEntity& host_entity = *host_entity_pointer;

    HostLocalFunction local_coarse_basis_func = coarse_basis_func.localFunction( host_entity );

    const Quadrature quadrature(local_grid_entity, 2 * discreteFunctionSpace.order() + 2);
    const size_t numQuadraturePoints = quadrature.nop();
    for (size_t quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint)
    {
      const typename Quadrature::CoordinateType& local_point = quadrature.point(quadraturePoint);

      const double weight = clement_weight * quadrature.weight(quadraturePoint) * geometry.integrationElement(local_point);

      std::vector<RangeType> fine_phi_x;
      baseSet.evaluateAll( quadrature[quadraturePoint], fine_phi_x);
      
      RangeType value_coarse_basis_func;
      local_coarse_basis_func.evaluate( quadrature[quadraturePoint] , value_coarse_basis_func);

      for (unsigned int i = 0; i < numBaseFunctions; ++i)
        elementOfRHS[i] += weight * value_coarse_basis_func * fine_phi_x[i];

    }
  }

} // assemble_local_RHS_pre_processing


void LocalProblemOperator
      ::assemble_local_RHS_lg_problems_all( const std::vector< std::shared_ptr<HostDiscreteFunction > >& coarse_basis_func_list,
                                            std::vector< double >& clement_weights,
                                            std::vector< int >& ids_basis_functions_in_subgrid,
                                            std::vector< std::unique_ptr< DiscreteFunction > >& local_problem_RHS ) const {

  typedef typename DiscreteFunction::DiscreteFunctionSpaceType DiscreteFunctionSpace;
  typedef typename DiscreteFunction::LocalFunctionType         LocalFunction;

  typedef typename DiscreteFunctionSpace::BaseFunctionSetType BaseFunctionSet;
  typedef typename DiscreteFunctionSpace::IteratorType        Iterator;
  typedef typename Iterator::Entity                           Entity;
  typedef typename Entity::Geometry                           Geometry;

  typedef typename DiscreteFunctionSpace::GridPartType GridPart;
  typedef CachingQuadrature< GridPart, 0 >             Quadrature;

  const DiscreteFunctionSpace& discreteFunctionSpace = local_problem_RHS[0]->space();

  const GridType& subGrid = discreteFunctionSpace.grid();
  
  for (int j = 0; j < local_problem_RHS.size(); ++j)
  {
    // set entries to zero:
    local_problem_RHS[j]->clear();
  }

  // gradient of micro scale base function:
  std::vector< JacobianRangeType > gradient_phi( discreteFunctionSpace.mapper().maxNumDofs() );

  const Iterator end = discreteFunctionSpace.end();
  for (Iterator it = discreteFunctionSpace.begin(); it != end; ++it)
  {
    const Entity& local_grid_entity = *it;
    const Geometry& geometry = local_grid_entity.geometry();
    assert(local_grid_entity.partitionType() == InteriorEntity);
    
    const BaseFunctionSet& baseSet = (local_problem_RHS[0]->localFunction(local_grid_entity)).baseFunctionSet();
    const auto numBaseFunctions = baseSet.size();
        
    HostEntityPointer host_entity_pointer = subGrid.getHostEntity< 0 >( local_grid_entity );
    const HostEntity& host_entity = *host_entity_pointer;

    const Quadrature quadrature(local_grid_entity, 2 * discreteFunctionSpace.order() + 2);
    const size_t numQuadraturePoints = quadrature.nop();
    for (size_t quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint)
    {
      const typename Quadrature::CoordinateType& local_point = quadrature.point(quadraturePoint);

      const double weight = quadrature.weight(quadraturePoint) * geometry.integrationElement(local_point);

      std::vector<RangeType> fine_phi_x;
      baseSet.evaluateAll( quadrature[quadraturePoint], fine_phi_x);
      
      for (int j = 0; j < local_problem_RHS.size(); ++j)
      {
         int interior_basis_func_id = ids_basis_functions_in_subgrid[j];
         HostLocalFunction local_coarse_basis_func = coarse_basis_func_list[interior_basis_func_id]->localFunction( host_entity );
         LocalFunction elementOfRHS = local_problem_RHS[j]->localFunction( local_grid_entity );

         RangeType value_coarse_basis_func;
         local_coarse_basis_func.evaluate( quadrature[quadraturePoint] , value_coarse_basis_func);

         for (unsigned int i = 0; i < numBaseFunctions; ++i)
           elementOfRHS[i] += clement_weights[interior_basis_func_id] * weight * value_coarse_basis_func * fine_phi_x[i];

      }
    }

  }

} // assemble_local_RHS_pre_processing_all


} //namespace MsFEM {
} //namespace Multiscale {
} //namespace Dune {
