// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include <config.h>

#include "localoperator.hh"

#include <dune/stuff/common/ranges.hh>
#include <dune/stuff/fem/localmatrix_proxy.hh>
#include <dune/stuff/discretefunction/projection/heterogenous.hh>
#include <dune/stuff/fem/matrix_object.hh>

#include <dune/multiscale/tools/misc/uzawa.hh>
#include <dune/multiscale/problems/base.hh>
#include <dune/multiscale/problems/selector.hh>

namespace Dune {
namespace Multiscale {
namespace MsFEM {

LocalProblemOperator::LocalProblemOperator(const DiscreteFunctionSpaceType& subDiscreteFunctionSpace,
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
void LocalProblemOperator::assemble_matrix(MsFEMLocalProblemSolver::LocProbFEMMatrixType& global_matrix) const
// x_T is the barycenter of the macro grid element T
{
  global_matrix.reserve(DSFe::diagonalAndNeighborStencil(global_matrix));
  global_matrix.clear();

  // local grid basis functions:
  std::vector< RangeType > phi( subDiscreteFunctionSpace_.mapper().maxNumDofs() );

  // gradient of micro scale base function:
  std::vector< typename BasisFunctionSetType::JacobianRangeType > gradient_phi(
    subDiscreteFunctionSpace_.mapper().maxNumDofs() );
  typename LocalFunctionType::JacobianRangeType diffusion_in_gradient_phi;

  const EntityIteratorType end = subDiscreteFunctionSpace_.end();
  for (EntityIteratorType it = subDiscreteFunctionSpace_.begin(); it != end; ++it)
  {
    const EntityType& sub_grid_entity = *it;
    const GeometryType& sub_grid_geometry = sub_grid_entity.geometry();

    DSFe::LocalMatrixProxy<MsFEMLocalProblemSolver::LocProbFEMMatrixType> local_matrix(global_matrix, sub_grid_entity, sub_grid_entity);

    const BasisFunctionSetType& baseSet = local_matrix.domainBasisFunctionSet();
    const auto numBaseFunctions = baseSet.size();

    // for constant diffusion "2*discreteFunctionSpace_.order()" is sufficient, for the general case, it is better to
    // use a higher order quadrature:
    const QuadratureType quadrature(sub_grid_entity, 2 * subDiscreteFunctionSpace_.order() + 2);
    const size_t numQuadraturePoints = quadrature.nop();
    for (size_t quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint)
    {
      // local (barycentric) coordinates (with respect to local grid entity)
      const typename QuadratureType::CoordinateType& local_point = quadrature.point(quadraturePoint);
      const DomainType global_point = sub_grid_geometry.global(local_point);

      const double weight = quadrature.weight(quadraturePoint)
                            * sub_grid_geometry.integrationElement(local_point);

      baseSet.jacobianAll(quadrature[quadraturePoint], gradient_phi);
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
  global_matrix.communicate();
} // assemble_matrix


//! stiffness matrix for a linear elliptic diffusion operator
void LocalProblemOperator::assemble_matrix(MsFEMLocalProblemSolver::LocProbFEMMatrixType& global_matrix,
                                           const SubGridList::CoarseNodeVectorType& coarse_node_vector ) const
// x_T is the barycenter of the macro grid element T
{
  global_matrix.reserve(DSFe::diagonalAndNeighborStencil(global_matrix));
  global_matrix.clear();

  // local grid basis functions:
  std::vector< RangeType > phi( subDiscreteFunctionSpace_.mapper().maxNumDofs() );

  // gradient of micro scale base function:
  std::vector< typename BasisFunctionSetType::JacobianRangeType > gradient_phi(
    subDiscreteFunctionSpace_.mapper().maxNumDofs() );

  const EntityIteratorType end = subDiscreteFunctionSpace_.end();
  for (EntityIteratorType it = subDiscreteFunctionSpace_.begin(); it != end; ++it)
  {
    const EntityType& sub_grid_entity = *it;
    const GeometryType& sub_grid_geometry = sub_grid_entity.geometry();
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

    DSFe::LocalMatrixProxy<MsFEMLocalProblemSolver::LocProbFEMMatrixType> local_matrix(global_matrix, sub_grid_entity, sub_grid_entity);

    const BasisFunctionSetType& baseSet = local_matrix.domainBasisFunctionSet();
    const auto numBaseFunctions = baseSet.size();

    std::vector<RangeType> value_phi(numBaseFunctions);
    // for constant diffusion "2*discreteFunctionSpace_.order()" is sufficient, for the general case, it is better to
    // use a higher order quadrature:
    const QuadratureType quadrature(sub_grid_entity, 2 * subDiscreteFunctionSpace_.order() + 2);
    const size_t numQuadraturePoints = quadrature.nop();
    for (size_t quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint)
    {
      // local (barycentric) coordinates (with respect to local grid entity)
      const typename QuadratureType::CoordinateType& local_point = quadrature.point(quadraturePoint);
      const DomainType global_point = sub_grid_geometry.global(local_point);

      const double weight = quadrature.weight(quadraturePoint)
                            * sub_grid_geometry.integrationElement(local_point);

      baseSet.jacobianAll(quadrature[quadraturePoint], gradient_phi);
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
        typename LocalFunctionType::JacobianRangeType diffusion_in_gradient_phi;
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
  global_matrix.communicate();
} // assemble_matrix

void LocalProblemOperator::set_zero_boundary_condition_RHS(const HostDiscreteFunctionSpaceType& host_space,
                                                           LocalProblemOperator::DiscreteFunction& rhs ) const {

  const DiscreteFunctionSpaceType& discreteFunctionSpace
    = rhs.space();

  const GridType& subGrid = discreteFunctionSpace.grid();
  const HostGridPart& hostGridPart = host_space.gridPart();
  
  // set Dirichlet Boundary to zero
  const EntityIteratorType endit = discreteFunctionSpace.end();
  for (EntityIteratorType it = discreteFunctionSpace.begin(); it != endit; ++it)
  {
    const EntityType& subgrid_entity = *it;

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
      else if ( iit->boundaryId() != 1 )
      { continue; }

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

  typedef typename DiscreteFunctionSpaceType::GridPartType        GridPartType;
  typedef typename DiscreteFunctionSpaceType::GridType            GridType;

  const DiscreteFunctionSpaceType& discreteFunctionSpace
    = rhs.space();

  const EntityIteratorType endit = discreteFunctionSpace.end();
  for (EntityIteratorType it = discreteFunctionSpace.begin(); it != endit; ++it)
  {
    // entity
    const EntityType& entity = *it;

    // create quadrature for given geometry type
    const Fem::CachingQuadrature< GridPartType, 0 > quadrature(entity, 2 * discreteFunctionSpace.order() + 2);

    // get geoemetry of entity
    const GeometryType& geo = entity.geometry();

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

  typedef typename DiscreteFunctionSpaceType::BasisFunctionSetType BaseFunctionSet;


  const DiscreteFunctionSpaceType& discreteFunctionSpace = local_problem_RHS.space();

  // set entries to zero:
  local_problem_RHS.clear();

  // gradient of micro scale base function:
  std::vector< JacobianRangeType > gradient_phi( discreteFunctionSpace.mapper().maxNumDofs() );

  const EntityIteratorType end = discreteFunctionSpace.end();
  for (EntityIteratorType it = discreteFunctionSpace.begin(); it != end; ++it)
  {
    const EntityType& local_grid_entity = *it;
    const GeometryType& geometry = local_grid_entity.geometry();

    LocalFunctionType elementOfRHS = local_problem_RHS.localFunction(local_grid_entity);

    const BaseFunctionSet& baseSet = elementOfRHS.basisFunctionSet();
    const auto numBaseFunctions = baseSet.size();

    const QuadratureType quadrature(local_grid_entity, 2 * discreteFunctionSpace.order() + 2);
    const size_t numQuadraturePoints = quadrature.nop();
    for (size_t quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint)
    {
      const typename QuadratureType::CoordinateType& local_point = quadrature.point(quadraturePoint);

      // remember, we are concerned with: - \int_{U(T)} (A^eps)(x) e · ∇ \phi(x)

      // global point in the subgrid
      const DomainType global_point = geometry.global(local_point);
      const double weight = quadrature.weight(quadraturePoint) * geometry.integrationElement(local_point);

      // A^eps(x) e
      // diffusion operator evaluated in 'x' multiplied with e
      JacobianRangeType diffusion_in_e;
      diffusion_operator_.diffusiveFlux(global_point, e, diffusion_in_e);

      baseSet.jacobianAll(quadrature[quadraturePoint], gradient_phi);
      for (unsigned int i = 0; i < numBaseFunctions; ++i)
      {
        elementOfRHS[i] -= weight * (diffusion_in_e[0] * gradient_phi[i][0]);
      }
    }
  }
} // assemble_local_RHS

/** Assemble right hand side vectors for all local problems on one coarse cell.
*
* @param[in] coarseEntity The coarse cell.
* @param[in] specifier A MacroMicroGridSpecifier (needed for access to the coarse base function set).
* @param[out] allLocalRHS A vector with pointers to the discrete functions for the right hand sides.
*
* @note The vector allLocalRHS is assumed to have the correct size and contain pointers to all local rhs
* functions. The discrete functions in allLocalRHS will be cleared in this function.
*/
void LocalProblemOperator
::assembleAllLocalRHS(const CoarseEntityType& coarseEntity, const MacroMicroGridSpecifierType& specifier,
                SubDiscreteFunctionVectorType& allLocalRHS) const {
  assert(allLocalRHS.size()>0 && "You need to preallocate the necessary space outside this function!");
  //! @todo correct the error message below (+1 for simplecial, +2 for arbitrary)
  assert((specifier.simplexCoarseGrid() && allLocalRHS.size()==GridType::dimension + 1) ||
         (!(specifier.simplexCoarseGrid()) && allLocalRHS.size() == specifier.fineSpace().mapper().maxNumDofs()+2) &&
          "You need to allocate storage space for the correctors for all unit vector/all coarse basis functions" &&
          " and the dirichlet- and neuman corrector");

  // build unit vectors (needed for cases where rhs is assembled for unit vectors instead of coarse
  // base functions)
  JacobianRangeType unitVectors[dimension];
  for (int i = 0; i < dimension; ++i) {
    for (int j = 0; j < dimension; ++j) {
      if (i == j) {
        unitVectors[i][0][j] = 1.0;
      } else {
        unitVectors[i][0][j] = 0.0;
      }
    }
  }

  // get dirichlet and neumann data
  auto dirichletDataPtr = Dune::Multiscale::Problem::getDirichletData();
  const auto& dirichletData = *dirichletDataPtr;
  auto neumannDataPtr = Dune::Multiscale::Problem::getNeumannData();
  const auto& neumannData = *neumannDataPtr;

  const DiscreteFunctionSpaceType& discreteFunctionSpace = allLocalRHS[0]->space();

  //! @todo we should use the dirichlet constraints here somehow
  SubDiscreteFunctionType dirichletExtension("dirichletExtension", discreteFunctionSpace);
  dirichletExtension.clear();
  HostDiscreteFunction dirichletExtensionCoarse("Dirichlet Extension Coarse", specifier.coarseSpace());
  dirichletExtensionCoarse.clear();
  this->projectDirichletValues(dirichletExtensionCoarse);
  Dune::Stuff::HeterogenousProjection<> projection;
  projection.project(dirichletExtensionCoarse, dirichletExtension);

  // set entries to zero:
  for (auto& rhs : allLocalRHS) rhs->clear();

  // get the base function set of the coarse space for the given coarse entity
  auto& coarseBaseSet = specifier.coarseSpace().basisFunctionSet(coarseEntity);
  std::vector< CoarseBaseFunctionSetType::JacobianRangeType > coarseBaseFuncJacs(coarseBaseSet.size());

  // gradient of micro scale base function:
  std::vector< JacobianRangeType > gradient_phi( discreteFunctionSpace.mapper().maxNumDofs() );
  std::vector< RangeType > phi( discreteFunctionSpace.mapper().maxNumDofs() );

  const int numBoundaryCorrectors = specifier.simplexCoarseGrid() ? 1 : 2;
  const int numInnerCorrectors = allLocalRHS.size() - numBoundaryCorrectors;

  for (auto& localGridCell : discreteFunctionSpace) {
    const GeometryType& geometry = localGridCell.geometry();
    const bool hasBoundaryIntersection = localGridCell.hasBoundaryIntersections();
    LocalFunctionType dirichletLF = dirichletExtension.localFunction(localGridCell);
    JacobianRangeType dirichletJac(0.0);

    for (int coarseBaseFunc=0; coarseBaseFunc<allLocalRHS.size(); ++coarseBaseFunc) {
      LocalFunctionType rhsLocalFunction = allLocalRHS[coarseBaseFunc]->localFunction(localGridCell);

      const BasisFunctionSetType& baseSet = rhsLocalFunction.basisFunctionSet();
      const auto numBaseFunctions = baseSet.size();

      if (coarseBaseFunc<numInnerCorrectors || coarseBaseFunc==numInnerCorrectors+1) {
        const QuadratureType quadrature(localGridCell, 2 * discreteFunctionSpace.order() + 2);
        const size_t numQuadraturePoints = quadrature.nop();
        for (size_t quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint) {
          const auto& local_point = quadrature.point(quadraturePoint);
          // global point in the subgrid
          const DomainType global_point = geometry.global(local_point);

          const double weight = quadrature.weight(quadraturePoint) * geometry.integrationElement(local_point);

          JacobianRangeType diffusion(0.0);
          if (coarseBaseFunc < numInnerCorrectors) {
            if (specifier.simplexCoarseGrid())
              diffusion_operator_.diffusiveFlux(global_point, unitVectors[coarseBaseFunc], diffusion);
            else {
              const DomainType quadInCoarseLocal = coarseEntity.geometry().local(global_point);
              coarseBaseSet.jacobianAll(quadInCoarseLocal, coarseBaseFuncJacs);
              diffusion_operator_.diffusiveFlux(global_point, coarseBaseFuncJacs[coarseBaseFunc], diffusion);
            }
          } else {
            dirichletLF.jacobian(local_point, dirichletJac);
            diffusion_operator_.diffusiveFlux(global_point, dirichletJac, diffusion);

          }
          baseSet.jacobianAll(quadrature[quadraturePoint], gradient_phi);
          for (unsigned int i = 0; i < numBaseFunctions; ++i) {
            rhsLocalFunction[i] -= weight * (diffusion[0] * gradient_phi[i][0]);
          }

        }
      }

      // boundary integrals
      if (coarseBaseFunc==numInnerCorrectors && hasBoundaryIntersection) {
        const auto intEnd = discreteFunctionSpace.gridPart().iend(localGridCell);
        for (auto iIt = discreteFunctionSpace.gridPart().ibegin(localGridCell); iIt!=intEnd; ++iIt) {
          const auto& intersection = *iIt;
          if (Dune::Multiscale::Problem::isNeumannBoundary(intersection)) {
            const int          orderOfIntegrand = (polynomialOrder - 1) + 2 * (polynomialOrder + 1);
            const int          quadOrder        = std::ceil((orderOfIntegrand + 1) / 2);
            FaceQuadratureType faceQuad(discreteFunctionSpace.gridPart(), intersection, quadOrder, FaceQuadratureType::INSIDE);
            RangeType neumannValue(0.0);
            const size_t numQuadPoints = faceQuad.nop();
            // loop over all quadrature points
            for (unsigned int iqP = 0; iqP < numQuadPoints; ++iqP) {
              // get local coordinate of quadrature point
              const auto& xLocal = faceQuad.localPoint(iqP);
              const auto& faceGeometry  = intersection.geometry();

              // the following does not work because subgrid does not implement geometryInInside()
              // const auto& insideGeometry    = intersection.geometryInInside();
              // const typename FaceQuadratureType::CoordinateType& xInInside = insideGeometry.global(xLocal);
              // therefore, we have to do stupid things:
              const typename FaceQuadratureType::CoordinateType& xGlobal = faceGeometry.global(xLocal);
              auto insidePtr = intersection.inside();
              const auto& insideEntity = *insidePtr;
              const auto& xInInside = insideEntity.geometry().local(xGlobal);

              const double factor = faceGeometry.integrationElement(xLocal) * faceQuad.weight(iqP);

              neumannData.evaluate(xGlobal, neumannValue);
              baseSet.evaluateAll(xInInside, phi);
              for (unsigned int i = 0; i < numBaseFunctions; ++i) {
                rhsLocalFunction[i] -= factor * (neumannValue * phi[i]);
              }
            }
          }
        }
      }
    }
  }

  return;
}


void LocalProblemOperator
      ::assemble_local_RHS(const JacobianRangeType &e, // direction 'e'
                           const SubGridList::CoarseNodeVectorType& coarse_node_vector, // for constraints on the space
                           const int& oversampling_strategy,
                           // rhs local msfem problem:
                           DiscreteFunction& local_problem_RHS) const {

  const DiscreteFunctionSpaceType& discreteFunctionSpace = local_problem_RHS.space();

  // set entries to zero:
  local_problem_RHS.clear();

  // gradient of micro scale base function:
  std::vector< JacobianRangeType > gradient_phi( discreteFunctionSpace.mapper().maxNumDofs() );

  const EntityIteratorType end = discreteFunctionSpace.end();
  for (EntityIteratorType it = discreteFunctionSpace.begin(); it != end; ++it)
  {
    const EntityType& local_grid_entity = *it;
    const GeometryType& geometry = local_grid_entity.geometry();
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
    LocalFunctionType elementOfRHS = local_problem_RHS.localFunction(local_grid_entity);

    const BasisFunctionSetType& baseSet = elementOfRHS.basisFunctionSet();
    const auto numBaseFunctions = baseSet.size();

    const QuadratureType quadrature(local_grid_entity, 2 * discreteFunctionSpace.order() + 2);
    const size_t numQuadraturePoints = quadrature.nop();
    for (size_t quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint)
    {
      const typename QuadratureType::CoordinateType& local_point = quadrature.point(quadraturePoint);

      // remember, we are concerned with: - \int_{U(T)} (A^eps)(x) e · ∇ \phi(x)
      // global point in the subgrid
      const DomainType global_point = geometry.global(local_point);
      const double weight = quadrature.weight(quadraturePoint) * geometry.integrationElement(local_point);

      // A^eps(x) e
      // diffusion operator evaluated in 'x' multiplied with e
      JacobianRangeType diffusion_in_e;
      diffusion_operator_.diffusiveFlux(global_point, e, diffusion_in_e);
      baseSet.jacobianAll(quadrature[quadraturePoint], gradient_phi);
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
      ::assemble_local_RHS_Dirichlet_corrector(
             const HostDiscreteFunction& dirichlet_extension,
             const SubGridList::CoarseNodeVectorType& coarse_node_vector, // for constraints on the space
             const int& oversampling_strategy,
             // rhs local msfem problem:
             DiscreteFunction& local_problem_RHS) const {

  typedef typename DiscreteFunction::DiscreteFunctionSpaceType DiscreteFunctionSpace;
  typedef typename DiscreteFunction::LocalFunctionType         LocalFunction;

  typedef typename DiscreteFunctionSpace::IteratorType        Iterator;
  typedef typename Iterator::Entity                           Entity;
  typedef typename Entity::Geometry                           Geometry;

  typedef typename DiscreteFunctionSpace::GridPartType GridPart;
  typedef Fem::CachingQuadrature< GridPart, 0 >             Quadrature;

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

    HostEntityPointer host_entity_pointer = subGrid.getHostEntity< 0 >( local_grid_entity );
    const HostEntity& host_entity = *host_entity_pointer;

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
    const HostLocalFunction loc_dirichlet_extension = dirichlet_extension.localFunction( host_entity );

    const BasisFunctionSetType& baseSet = elementOfRHS.basisFunctionSet();
    const auto numBaseFunctions = baseSet.size();

    const Quadrature quadrature(local_grid_entity, 2 * discreteFunctionSpace.order() + 2);
    const size_t numQuadraturePoints = quadrature.nop();
    for (size_t quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint)
    {
      const typename Quadrature::CoordinateType& local_point = quadrature.point(quadraturePoint);

      // remember, we are concerned with: - \int_T A( \nabla dirichlet_extension ) · ∇ \phi(x)
      // global point in the subgrid
      const DomainType global_point = geometry.global(local_point);
      const double weight = quadrature.weight(quadraturePoint) * geometry.integrationElement(local_point);

      JacobianRangeType gradient_dirichlet_extension;
      loc_dirichlet_extension.jacobian( quadrature[quadraturePoint] , gradient_dirichlet_extension);

      // A( \nabla dirichlet_extension )
      // diffusion operator evaluated in 'x' multiplied with gradient_dirichlet_extension
      JacobianRangeType diffusive_flux;

      diffusion_operator_.diffusiveFlux(global_point, gradient_dirichlet_extension , diffusive_flux );
      baseSet.jacobianAll(quadrature[quadraturePoint], gradient_phi);
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
          elementOfRHS[i] -= weight * (diffusive_flux[0] * gradient_phi[i][0]);
        else
          elementOfRHS[i] = 0.0;
      }

    }
  }

} // assemble_local_RHS_Dirichlet_corrector


void LocalProblemOperator
      ::assemble_local_RHS_Neumann_corrector(
             const NeumannBoundaryType& neumann_bc,
             const HostDiscreteFunctionSpaceType& host_space,
             const SubGridList::CoarseNodeVectorType& coarse_node_vector, // for constraints on the space
             const int& oversampling_strategy,
             // rhs local msfem problem:
             DiscreteFunction& local_problem_RHS) const {

  typedef typename DiscreteFunction::DiscreteFunctionSpaceType DiscreteFunctionSpace;
  typedef typename DiscreteFunction::LocalFunctionType         LocalFunction;

  typedef typename DiscreteFunctionSpace::IteratorType        Iterator;
  typedef typename Iterator::Entity                           Entity;
  typedef typename Entity::Geometry                           Geometry;

  typedef typename DiscreteFunctionSpace::GridPartType GridPart;
  typedef Fem::CachingQuadrature< GridPart, 0 >             Quadrature;

  const DiscreteFunctionSpace& discreteFunctionSpace = local_problem_RHS.space();
  const GridType& subGrid = discreteFunctionSpace.grid();

  // set entries to zero:
  local_problem_RHS.clear();

  // gradient of micro scale base function:
  std::vector<RangeType> phi( discreteFunctionSpace.mapper().maxNumDofs() );

  const Iterator end = discreteFunctionSpace.end();
  for (Iterator it = discreteFunctionSpace.begin(); it != end; ++it)
  {
    const Entity& local_grid_entity = *it;
    const Geometry& geometry = local_grid_entity.geometry();
    assert(local_grid_entity.partitionType() == InteriorEntity);

    HostEntityPointer host_entity_pointer = subGrid.getHostEntity< 0 >( local_grid_entity );
    const HostEntity& host_entity = *host_entity_pointer;

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
         { sub_grid_entity_corner_is_relevant.push_back(c); }
        }
      }
    }

    LocalFunction elementOfRHS = local_problem_RHS.localFunction(local_grid_entity);
    const BasisFunctionSetType& baseSet = elementOfRHS.basisFunctionSet();

    const auto& lagrangePointSet = host_space.lagrangePointSet( host_entity );

    std::vector< std::vector<RangeType> > phi_values(sub_grid_entity_corner_is_relevant.size());
    for(auto j : DSC::valueRange(sub_grid_entity_corner_is_relevant.size())) {
        baseSet.evaluateAll(geometry.local(geometry.corner(sub_grid_entity_corner_is_relevant[j])), phi_values[j]);
    }

    for (const auto& intersection
         : Dune::Stuff::Common::intersectionRange(host_space.gridPart(), host_entity ))
      {
        if ( !intersection.boundary() )
          continue;
        // boundaryId 1 = Dirichlet face; boundaryId 2 = Neumann face;
        if ( intersection.boundary() && (intersection.boundaryId() != 2) )
          continue;

        const auto face = intersection.indexInInside();

        const HostFaceQuadrature faceQuadrature( host_space.gridPart(),
                                                 intersection, 2 * host_space.order() + 2, HostFaceQuadrature::INSIDE );
        const int numFaceQuadraturePoints = faceQuadrature.nop();

        enum { faceCodim = 1 };
        for (int faceQuadraturePoint = 0; faceQuadraturePoint < numFaceQuadraturePoints; ++faceQuadraturePoint)
        {
          baseSet.evaluateAll( faceQuadrature[faceQuadraturePoint], phi );

          const auto local_point_entity = faceQuadrature.point( faceQuadraturePoint ); 
          const auto global_point = host_entity.geometry().global( local_point_entity ); 
          const auto local_point_face = intersection.geometry().local( global_point );

          RangeType neumann_value( 0.0 );
          neumann_bc.evaluate( global_point, neumann_value );

          const double face_weight = intersection.geometry().integrationElement( local_point_face )
                          * faceQuadrature.weight( faceQuadraturePoint );

          auto faceIterator = lagrangePointSet.beginSubEntity< faceCodim >( face );
          const auto faceEndIterator = lagrangePointSet.endSubEntity< faceCodim >( face );

          for ( ; faceIterator != faceEndIterator; ++faceIterator)
          {
             bool zero_entry = false;
             for ( size_t sgec = 0; sgec < sub_grid_entity_corner_is_relevant.size(); ++sgec )
             {
               const auto& value_phi_i = phi_values[sgec][*faceIterator];
               if ( value_phi_i == 1.0 )
               {
                 zero_entry = true;
               }
             }
             
             if (!zero_entry)
                elementOfRHS[ *faceIterator ] -= neumann_value * face_weight * phi[ *faceIterator ];
             else
                elementOfRHS[ *faceIterator ] = 0.0;
          }

        }

      }
  }
} // assemble_local_RHS_Neumann_corrector



void LocalProblemOperator
      ::assemble_local_RHS_lg_problems( const HostDiscreteFunction& coarse_basis_func, double clement_weight,
                                           DiscreteFunction& local_problem_RHS ) const {

  const DiscreteFunctionSpaceType& discreteFunctionSpace = local_problem_RHS.space();

  const GridType& subGrid = discreteFunctionSpace.grid();
  
  // set entries to zero:
  local_problem_RHS.clear();

  // gradient of micro scale base function:
  std::vector< JacobianRangeType > gradient_phi( discreteFunctionSpace.mapper().maxNumDofs() );

  const EntityIteratorType end = discreteFunctionSpace.end();
  for (EntityIteratorType it = discreteFunctionSpace.begin(); it != end; ++it)
  {
    const EntityType& local_grid_entity = *it;
    const GeometryType& geometry = local_grid_entity.geometry();
    assert(local_grid_entity.partitionType() == InteriorEntity);

    LocalFunctionType elementOfRHS = local_problem_RHS.localFunction(local_grid_entity);

    const BasisFunctionSetType& baseSet = elementOfRHS.basisFunctionSet();
    const auto numBaseFunctions = baseSet.size();
        
    HostEntityPointer host_entity_pointer = subGrid.getHostEntity< 0 >( local_grid_entity );
    const HostEntity& host_entity = *host_entity_pointer;

    HostLocalFunction local_coarse_basis_func = coarse_basis_func.localFunction( host_entity );

    const QuadratureType quadrature(local_grid_entity, 2 * discreteFunctionSpace.order() + 2);
    const size_t numQuadraturePoints = quadrature.nop();
    for (size_t quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint)
    {
      const typename QuadratureType::CoordinateType& local_point = quadrature.point(quadraturePoint);

      const double weight = clement_weight * quadrature.weight(quadraturePoint) * geometry.integrationElement(local_point);

      std::vector<RangeType> fine_phi_x( discreteFunctionSpace.mapper().maxNumDofs() );
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
                                            std::vector< std::unique_ptr< LocalProblemOperator::DiscreteFunction > >& local_problem_RHS ) const {
  const DiscreteFunctionSpaceType& discreteFunctionSpace = local_problem_RHS[0]->space();
  
  const GridType& subGrid = discreteFunctionSpace.grid();

  for (auto& rhs : local_problem_RHS)
  {
    rhs->clear();
  }

  for (const EntityType& local_grid_entity : discreteFunctionSpace)
  {
    const GeometryType& geometry = local_grid_entity.geometry();
    assert(local_grid_entity.partitionType() == InteriorEntity);
    
    const BasisFunctionSetType& baseSet = (local_problem_RHS[0]->localFunction(local_grid_entity)).basisFunctionSet();
    const auto numBaseFunctions = baseSet.size();

    HostEntityPointer host_entity_pointer = subGrid.getHostEntity< 0 >( local_grid_entity );
    const HostEntity& host_entity = *host_entity_pointer;

    const QuadratureType quadrature(local_grid_entity, 2 * discreteFunctionSpace.order() + 2);
    const size_t numQuadraturePoints = quadrature.nop();
    for (size_t quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint)
    {
      const typename QuadratureType::CoordinateType& local_point = quadrature.point(quadraturePoint);

      const double weight = quadrature.weight(quadraturePoint) * geometry.integrationElement(local_point);

      std::vector<RangeType> fine_phi_x( discreteFunctionSpace.mapper().maxNumDofs() );
      baseSet.evaluateAll( quadrature[quadraturePoint], fine_phi_x);

      for (std::size_t j = 0; j < local_problem_RHS.size(); ++j)
      {
         int interior_basis_func_id = ids_basis_functions_in_subgrid[j];

         HostLocalFunction local_coarse_basis_func = coarse_basis_func_list[interior_basis_func_id]->localFunction( host_entity );
         LocalFunctionType elementOfRHS = local_problem_RHS[j]->localFunction( local_grid_entity );

         RangeType value_coarse_basis_func;
         local_coarse_basis_func.evaluate( quadrature[quadraturePoint] , value_coarse_basis_func);

         for (unsigned int i = 0; i < numBaseFunctions; ++i)
            elementOfRHS[i] += clement_weights[interior_basis_func_id] * weight * value_coarse_basis_func * fine_phi_x[i];
      }
    }
  }

} // assemble_local_RHS_pre_processing_all

/** Set the dirichlet values to a given discrete function on the sub mesh
*
* @param[in, out] function The function in which the values will be set.
*/
void LocalProblemOperator::projectDirichletValues(HostDiscreteFunction& function) const {
/*  // make sure, we are on a hexahedral element
  assert(function.space().gridPart().grid().leafIndexSet().geomTypes(0).size()==1 &&
         function.space().gridPart().grid().leafIndexSet().geomTypes(0)[0].isCube() &&
         "This method only works for hexahedral elements at the moment!");*/

  const auto& gridPart = function.space().gridPart();
  auto dirichletDataPtr = Multiscale::Problem::getDirichletData();
  const auto& dirichletData = *dirichletDataPtr;
//  Fem::GridFunctionAdapter<Multiscale::Problem::DirichletDataBase,
//  typename SubDiscreteFunctionType::DiscreteFunctionSpaceType::GridPartType> gf("dirichlet", dirichletData , gridPart);
  for (const auto& localCell : function.space()) {
    if (localCell.hasBoundaryIntersections())
      for (const auto& intersection : DSC::intersectionRange(gridPart, localCell)) {
        if (Multiscale::Problem::isDirichletBoundary(intersection)) {
          auto funcLocal = function.localFunction(localCell);
          const auto& lagrangePointSet = function.space().lagrangePointSet(localCell);
          const auto faceNumber = intersection.indexInInside();
          for (auto lp : DSC::lagrangePointSetRange<1>(function.space(), localCell, faceNumber)) {
            auto lagrangePoint = lagrangePointSet.point(lp);
            DomainType lagrangePointGlobal = localCell.geometry().global(lagrangePoint);
            RangeType dirichletVal(0.0);
            dirichletData.evaluate(lagrangePointGlobal, dirichletVal);
            funcLocal[lp] = dirichletVal;
          }
        }
      }
  }

  return;
}


} //namespace MsFEM {
} //namespace Multiscale {
} //namespace Dune {
