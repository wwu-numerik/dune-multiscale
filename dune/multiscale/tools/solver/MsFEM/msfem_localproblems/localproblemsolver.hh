#ifndef DiscreteEllipticMsFEMLocalProblem_HH
#define DiscreteEllipticMsFEMLocalProblem_HH

#include <vector>

#include <dune/fem/operator/matrix/spmatrix.hh>
#include <dune/common/fmatrix.hh>
#include <dune/stuff/common/filesystem.hh>
#include <dune/stuff/fem/functions/checks.hh>

#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/operator/common/operator.hh>

// LOCPROBLEMSOLVER_VERBOSE: 0 = false, 1 = true
#define LOCPROBLEMSOLVER_VERBOSE false

// dune-subgrid include:
#include <dune/subgrid/subgrid.hh>

// dune-fem includes:
#include <dune/fem/gridpart/common/gridpart.hh>
#include <dune/fem/operator/2order/lagrangematrixsetup.hh>

#include <dune/multiscale/tools/disc_func_writer/discretefunctionwriter.hh>
#include <dune/multiscale/tools/misc/outputparameter.hh>
#include <dune/multiscale/msfem/msfem_traits.hh>

namespace Dune {
/** \brief define output parameters for local problems
 *  appends "local_problems" for path
 **/
struct LocalProblemDataOutputParameters
  : public myDataOutputParameters
{
public:
  explicit LocalProblemDataOutputParameters()
    : myDataOutputParameters(DSC_CONFIG_GET("global.datadir", "data") + "/local_problems/")
  {}
};

// Imp stands for Implementation
template< class SubDiscreteFunctionType, class DiffusionOperatorType >
class LocalProblemOperator
  : public Operator< typename SubDiscreteFunctionType::RangeFieldType,
                     typename SubDiscreteFunctionType::RangeFieldType,
                     SubDiscreteFunctionType,
                     SubDiscreteFunctionType >
{
  typedef LocalProblemOperator< SubDiscreteFunctionType, DiffusionOperatorType > This;

public:
  typedef SubDiscreteFunctionType DiscreteFunction;
  typedef DiffusionOperatorType   DiffusionModel;

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
  LocalProblemOperator(const DiscreteFunctionSpace& subDiscreteFunctionSpace, const DiffusionModel& diffusion_op)
    : subDiscreteFunctionSpace_(subDiscreteFunctionSpace)
      , diffusion_operator_(diffusion_op)
  {}

private:
  LocalProblemOperator(const This&);

public:
  // dummy operator
  virtual void operator()(const DiscreteFunction& u, DiscreteFunction& w) const;

  // assemble stiffness matrix for local problems (oversampling strategy 1)
  template< class MatrixType >
  void assemble_matrix(MatrixType& global_matrix) const;

  // assemble stiffness matrix for local problems (oversampling strategy 2 and 3)
  template< class MatrixType, class CoarseNodeVectorType >
  void assemble_matrix(MatrixType& global_matrix, const CoarseNodeVectorType& coarse_node_vector /*for constraints*/) const;
  
  // the right hand side assembler methods
  void assemble_local_RHS(  // direction 'e'
    const JacobianRangeType& e,
    // rhs local msfem problem:
    DiscreteFunction& local_problem_RHS) const;

  // the right hand side assembler methods
  template< class CoarseNodeVectorType >
  void assemble_local_RHS(  // direction 'e'
    const JacobianRangeType& e,
    const CoarseNodeVectorType& coarse_node_vector, /*for constraints*/
    // rhs local msfem problem:
    DiscreteFunction& local_problem_RHS) const;
    
  void printLocalRHS(const DiscreteFunction& rhs) const;

  double normRHS(const DiscreteFunction& rhs) const;

private:
  const DiscreteFunctionSpace& subDiscreteFunctionSpace_;
  const DiffusionModel& diffusion_operator_;
};

// dummy implementation of "operator()"
// 'w' = effect of the discrete operator on 'u'
template< class DiscreteFunctionImp, class DiffusionImp >
void LocalProblemOperator< DiscreteFunctionImp, DiffusionImp >::operator()(const DiscreteFunctionImp& /*u*/,
                                                                           DiscreteFunctionImp& /*w*/) const {
  DUNE_THROW(Dune::NotImplemented,"the ()-operator of the LocalProblemOperator class is not yet implemented and still a dummy.");
}



//! stiffness matrix for a linear elliptic diffusion operator
// for oversampling strategy 1 (no constraints)
template< class SubDiscreteFunctionImp, class DiffusionImp >
template< class MatrixType >
void LocalProblemOperator< SubDiscreteFunctionImp, DiffusionImp >::assemble_matrix(MatrixType& global_matrix) const
// x_T is the barycenter of the macro grid element T
{
  typedef typename MatrixType::LocalMatrixType LocalMatrix;

  Problem::ModelProblemData model_info;

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

    LocalMatrix local_matrix = global_matrix.localMatrix(sub_grid_entity, sub_grid_entity);

    const BaseFunctionSet& baseSet = local_matrix.domainBaseFunctionSet();
    const unsigned int numBaseFunctions = baseSet.size();

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
      const FieldMatrix< double, dimension, dimension >& inverse_jac
        = sub_grid_geometry.jacobianInverseTransposed(local_point);

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




//! stiffness matrix for a linear elliptic diffusion operator
template< class SubDiscreteFunctionImp, class DiffusionImp >
template< class MatrixType, class CoarseNodeVectorType >
void LocalProblemOperator< SubDiscreteFunctionImp, DiffusionImp >::assemble_matrix(MatrixType& global_matrix, const CoarseNodeVectorType& coarse_node_vector ) const
// x_T is the barycenter of the macro grid element T
{
  typedef typename MatrixType::LocalMatrixType LocalMatrix;

  Problem::ModelProblemData model_info;

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
      for ( int coarse_node_local_id = 0; coarse_node_local_id < coarse_node_vector.size(); ++coarse_node_local_id )
       {
	 // if the subgrid corner 'c' is in the 'relevant coarse node vector' and if 'c' was not yet added to the
	 // vector 'sub_grid_entity_corner_is_relevant' then add it to the vector
         if ( (coarse_node_vector[coarse_node_local_id] == sub_grid_geometry.corner(c)) 
	     && (std::find(sub_grid_entity_corner_is_relevant.begin(), sub_grid_entity_corner_is_relevant.end(), c) == sub_grid_entity_corner_is_relevant.end()) )
	     { sub_grid_entity_corner_is_relevant.push_back(c); }
       }
    }
  
//! delete me:
/*
    for ( int coarse_node_local_id = 0; coarse_node_local_id < subgrid_list_.getCoarseNodeVector( coarse_index ).size(); ++coarse_node_local_id )
       {
         std::cout << coarse_node_local_id+1 << " : " << subgrid_list_.getCoarseNodeVector( coarse_index )[coarse_node_local_id] << std :: endl << std :: endl << std :: endl;
       }
*/
    
    LocalMatrix local_matrix = global_matrix.localMatrix(sub_grid_entity, sub_grid_entity);

    const BaseFunctionSet& baseSet = local_matrix.domainBaseFunctionSet();
    const unsigned int numBaseFunctions = baseSet.size();

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
      const FieldMatrix< double, dimension, dimension >& inverse_jac
        = sub_grid_geometry.jacobianInverseTransposed(local_point);

      for (unsigned int i = 0; i < numBaseFunctions; ++i)
      {
        // jacobian of the base functions, with respect to the reference element
        typename BaseFunctionSet::JacobianRangeType gradient_phi_ref_element;
        baseSet.jacobian(i, quadrature[quadraturePoint], gradient_phi_ref_element);

        // multiply it with transpose of jacobian inverse to obtain the jacobian with respect to the real entity
        inverse_jac.mv(gradient_phi_ref_element[0], gradient_phi[i][0]);

        baseSet.evaluate(i, quadrature[quadraturePoint], phi[i]);

        for ( int sgec = 0; sgec < sub_grid_entity_corner_is_relevant.size(); ++sgec )
        {
          RangeType value_phi_i(0.0);
          baseSet.evaluate(i, sub_grid_geometry.local(sub_grid_geometry.corner(sgec)), value_phi_i);
          if ( value_phi_i == 1.0 )
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

template< class DiscreteFunctionImp, class DiffusionImp >
void LocalProblemOperator< DiscreteFunctionImp, DiffusionImp >::printLocalRHS(const DiscreteFunctionImp& rhs) const {
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
      DSC_LOG_DEBUG << "Number of Dof: " << i << " ; " << rhs.name() << " : " << elementOfRHS[i] << std::endl;
    }
  }
}  // end method

template< class DiscreteFunctionImp, class DiffusionImp >
double LocalProblemOperator< DiscreteFunctionImp, DiffusionImp >::normRHS(const DiscreteFunctionImp& rhs) const {
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

// assemble the right hand side of a local problem (reconstruction problem on entity)
// ----------------------------------------------

// assemble method for the case of a linear diffusion operator

// we compute the following entries for each fine-scale base function phi_h_i:
// - \int_{T_0} (A^eps ○ F)(x) ∇ \Phi_H(x_T) · ∇ \phi_h_i(x)
template< class DiscreteFunctionImp, class DiffusionImp >
// template< class MatrixType >
void LocalProblemOperator< DiscreteFunctionImp, DiffusionImp >
      ::assemble_local_RHS(const JacobianRangeType &e, // direction 'e'
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

    LocalFunction elementOfRHS = local_problem_RHS.localFunction(local_grid_entity);

    const BaseFunctionSet& baseSet = elementOfRHS.baseFunctionSet();
    const unsigned int numBaseFunctions = baseSet.size();

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
        elementOfRHS[i] -= weight * (diffusion_in_e[0] * gradient_phi[i][0]);
      }
    }
  }
} // assemble_local_RHS



// assemble method for the case of a linear diffusion operator
// in a constraint space, for oversampling strategy 2

// we compute the following entries for each fine-scale base function phi_h_i:
// - \int_{T_0} (A^eps ○ F)(x) ∇ \Phi_H(x_T) · ∇ \phi_h_i(x)
template< class DiscreteFunctionImp, class DiffusionImp >
template< class CoarseNodeVectorType >
void LocalProblemOperator< DiscreteFunctionImp, DiffusionImp >
      ::assemble_local_RHS(const JacobianRangeType &e, // direction 'e'
                           const CoarseNodeVectorType& coarse_node_vector, // for constraints on the space
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

    std::vector< int > sub_grid_entity_corner_is_relevant;
    for ( int c = 0; c < geometry.corners(); ++c )
    {
      for ( int coarse_node_local_id = 0; coarse_node_local_id < coarse_node_vector.size(); ++coarse_node_local_id )
       {
	 // if the subgrid corner 'c' is in the 'relevant coarse node vector' and if 'c' was not yet added to the
	 // vector 'sub_grid_entity_corner_is_relevant' then add it to the vector
         if ( (coarse_node_vector[coarse_node_local_id] == geometry.corner(c)) 
	     && (std::find(sub_grid_entity_corner_is_relevant.begin(), sub_grid_entity_corner_is_relevant.end(), c) == sub_grid_entity_corner_is_relevant.end()) )
	     { sub_grid_entity_corner_is_relevant.push_back(c); }
       }
    }

    LocalFunction elementOfRHS = local_problem_RHS.localFunction(local_grid_entity);

    const BaseFunctionSet& baseSet = elementOfRHS.baseFunctionSet();
    const unsigned int numBaseFunctions = baseSet.size();

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

      for (unsigned int i = 0; i < numBaseFunctions; ++i)
      {
        // jacobian of the base functions, with respect to the reference element
        JacobianRangeType gradient_phi_ref_element;
        baseSet.jacobian(i, quadrature[quadraturePoint], gradient_phi_ref_element);

        // multiply it with transpose of jacobian inverse to obtain the jacobian with respect to the real entity
        inverse_jac.mv(gradient_phi_ref_element[0], gradient_phi[i][0]);

        for ( int sgec = 0; sgec < sub_grid_entity_corner_is_relevant.size(); ++sgec )
        {
          RangeType value_phi_i(0.0);
          baseSet.evaluate(i, geometry.local(geometry.corner(sgec)), value_phi_i);
          if ( value_phi_i == 1.0 )
          {
	    assert( dimension == 2);
            gradient_phi[i][0][0] = 0.0;
            gradient_phi[i][0][1] = 0.0;
	  }
        }
        
      }

      for (unsigned int i = 0; i < numBaseFunctions; ++i)
      {
        elementOfRHS[i] -= weight * (diffusion_in_e[0] * gradient_phi[i][0]);
      }
    }
  }
} // assemble_local_RHS


//! ------------------------------------------------------------------------------------------------
//! ------------------------------------------------------------------------------------------------

//! ------------------------------------------------------------------------------------------------

//! ------------------------------------------------------------------------------------------------
//! --------------------- the essential local msfem problem solver class ---------------------------

template< class HostDiscreteFunctionType,
          class SubGridListType,
          class MacroMicroGridSpecifierType,
          class DiffusionOperatorType >
class MsFEMLocalProblemSolver
{
public:
  //! ---------------- typedefs for the HostDiscreteFunctionSpace -----------------------

  //! type of discrete function space
  typedef typename HostDiscreteFunctionType::DiscreteFunctionSpaceType
  HostDiscreteFunctionSpaceType;

  //! type of (non-discrete )function space
  typedef typename HostDiscreteFunctionSpaceType::FunctionSpaceType FunctionSpaceType;

  //! type of grid partition
  typedef typename HostDiscreteFunctionSpaceType::GridPartType HostGridPartType;

  //! type of grid
  typedef typename HostDiscreteFunctionSpaceType::GridType HostGridType;

  //! type of range vectors
  typedef typename HostDiscreteFunctionSpaceType::RangeType RangeType;

  //! type of value of a gradient of a function
  typedef typename HostDiscreteFunctionSpaceType::JacobianRangeType JacobianRangeType;

  //! type of range vectors
  typedef typename HostDiscreteFunctionSpaceType::DomainType DomainType;

  typedef typename HostGridType::Traits::LeafIndexSet HostGridLeafIndexSet;

  typedef typename HostDiscreteFunctionSpaceType::IteratorType HostGridEntityIteratorType;

  typedef typename HostGridEntityIteratorType::Entity HostEntityType;

  typedef typename HostEntityType::EntityPointer HostEntityPointerType;

  typedef typename HostGridType::template Codim< 0 >::Geometry HostGridEntityGeometry;

  typedef typename HostDiscreteFunctionType::LocalFunctionType HostLocalFunctionType;

  typedef typename HostGridPartType::IntersectionIteratorType HostIntersectionIterator;

  //! ---------------- typedefs for the SubgridDiscreteFunctionSpace -----------------------
  // ( typedefs for the local grid and the corresponding local ('sub') )discrete space )

  //! type of grid
  typedef typename SubGridListType::SubGridType SubGridType;

  //! type of grid part
  typedef LeafGridPart< SubGridType > SubGridPartType;

  //! type of subgrid discrete function space
  typedef LagrangeDiscreteFunctionSpace< FunctionSpaceType, SubGridPartType, 1 >  // 1=POLORDER
  SubDiscreteFunctionSpaceType;

  //! type of subgrid discrete function
  typedef AdaptiveDiscreteFunction< SubDiscreteFunctionSpaceType > SubDiscreteFunctionType;

  typedef typename SubDiscreteFunctionSpaceType::IteratorType SubgridIteratorType;

  typedef typename SubgridIteratorType::Entity SubgridEntityType;

  typedef typename SubgridEntityType::EntityPointer SubgridEntityPointerType;

  typedef typename SubDiscreteFunctionType::LocalFunctionType SubLocalFunctionType;

  typedef typename SubDiscreteFunctionSpaceType::LagrangePointSetType SGLagrangePointSetType;

  //!-----------------------------------------------------------------------------------------

  //! ------------------ Matrix Traits for the local Problems ---------------------

  typedef typename SubDiscreteFunctionSpaceType::LagrangePointSetType SubgridLagrangePointSetType;

  enum { faceCodim = 1 };
  typedef typename SubgridLagrangePointSetType::template Codim< faceCodim >::SubEntityIteratorType
  SubgridFaceDofIteratorType;

  //! polynomial order of base functions
  enum { polynomialOrder = SubDiscreteFunctionSpaceType::polynomialOrder };

  struct LocProbMatrixTraits
  {
    typedef SubDiscreteFunctionSpaceType                          RowSpaceType;
    typedef SubDiscreteFunctionSpaceType                          ColumnSpaceType;
    typedef LagrangeMatrixSetup< false >                          StencilType;
    typedef ParallelScalarProduct< SubDiscreteFunctionSpaceType > ParallelScalarProductType;

    template< class M >
    struct Adapter
    {
      typedef LagrangeParallelMatrixAdapter< M > MatrixAdapterType;
    };
  };

  typedef SparseRowMatrixOperator< SubDiscreteFunctionType, SubDiscreteFunctionType,
                                   LocProbMatrixTraits > LocProbFEMMatrix;

  #ifdef SYMMETRIC_DIFFUSION_MATRIX
  typedef CGInverseOperator< SubDiscreteFunctionType, LocProbFEMMatrix > InverseLocProbFEMMatrix;
  #else
  // OEMGMRESOp //OEMBICGSQOp // OEMBICGSTABOp
  typedef OEMBICGSTABOp< SubDiscreteFunctionType, LocProbFEMMatrix > InverseLocProbFEMMatrix;
  #endif // ifdef SYMMETRIC_DIFFUSION_MATRIX

  // discrete elliptic operator describing the elliptic local msfem problems
  typedef LocalProblemOperator< SubDiscreteFunctionType, DiffusionOperatorType > LocalProblemOperatorType;

private:
  const HostDiscreteFunctionSpaceType& hostDiscreteFunctionSpace_;
  const DiffusionOperatorType& diffusion_;
  const MacroMicroGridSpecifierType& specifier_;
  SubGridListType& subgrid_list_;

public:
  /** \brief constructor - with diffusion operator A^{\epsilon}(x)
   * \param subgrid_list cannot be const because Dune::Fem does not provide Gridparts that can be build on a const grid
   * \param DSC_LOG_INFO does not take ownership
   **/
  MsFEMLocalProblemSolver(const HostDiscreteFunctionSpaceType& hostDiscreteFunctionSpace,
                          const MacroMicroGridSpecifierType& specifier,
                          SubGridListType& subgrid_list,
                          const DiffusionOperatorType& diffusion_operator)
    : hostDiscreteFunctionSpace_(hostDiscreteFunctionSpace)
      , diffusion_(diffusion_operator)
      , specifier_(specifier)
      , subgrid_list_(subgrid_list)
  {}

  template< class Stream >
  void oneLinePrint(Stream& stream, const SubDiscreteFunctionType& func) const {
    typedef typename SubDiscreteFunctionType::ConstDofIteratorType
    DofIteratorType;
    DofIteratorType it = func.dbegin();
    stream << "\n" << func.name() << ": [ ";
    for ( ; it != func.dend(); ++it)
      stream << std::setw(5) << *it << "  ";

    stream << " ] " << std::endl;
  } // oneLinePrint

  //! ----------- method: solve the local MsFEM problem ------------------------------------------

  void solvelocalproblem(JacobianRangeType& e,
                         SubDiscreteFunctionType& local_problem_solution,
                         const int coarse_index = -1 ) const {
    // set solution equal to zero:
    local_problem_solution.clear();

    const SubDiscreteFunctionSpaceType& subDiscreteFunctionSpace = local_problem_solution.space();

    //! the matrix in our linear system of equations
    // in the non-linear case, it is the matrix for each iteration step
    LocProbFEMMatrix locprob_system_matrix("Local Problem System Matrix",
                                           subDiscreteFunctionSpace,
                                           subDiscreteFunctionSpace);

    //! define the discrete (elliptic) local MsFEM problem operator
    // ( effect of the discretized differential operator on a certain discrete function )
    LocalProblemOperatorType local_problem_op(subDiscreteFunctionSpace, diffusion_);

    const SubGridType& subGrid = subDiscreteFunctionSpace.grid();

    typedef typename SubDiscreteFunctionSpaceType::IteratorType SGIteratorType;
    typedef typename SubGridPartType::IntersectionIteratorType  SGIntersectionIteratorType;

    //! right hand side vector of the algebraic local MsFEM problem
    SubDiscreteFunctionType local_problem_rhs("rhs of local MsFEM problem", subDiscreteFunctionSpace);
    local_problem_rhs.clear();

    // NOTE:
    // is the right hand side of the local MsFEM problem equal to zero or almost identical to zero?
    // if yes, the solution of the local MsFEM problem is also identical to zero. The solver is getting a problem with
    // this situation, which is why we do not solve local msfem problems for zero-right-hand-side, since we already know
    // the result.
    
    // assemble the stiffness matrix
    if ( DSC_CONFIG_GET( "msfem.oversampling_strategy", 1 ) == 1 )
      { local_problem_op.assemble_matrix(locprob_system_matrix); }
    else if ( ( DSC_CONFIG_GET( "msfem.oversampling_strategy", 1 ) == 2 ) ||
              ( DSC_CONFIG_GET( "msfem.oversampling_strategy", 1 ) == 3 ) )
      { if ( coarse_index < 0 )
          DUNE_THROW(Dune::InvalidStateException, "Invalid coarse index: coarse_index < 0");
        local_problem_op.assemble_matrix(locprob_system_matrix, subgrid_list_.getCoarseNodeVector( coarse_index ) ); }
    else
      DUNE_THROW( Dune::InvalidStateException, "Oversampling Strategy must be 1, 2 or 3!");

// can be deleted (just to check the coarse node vector)
/*
    for ( int coarse_node_local_id = 0; coarse_node_local_id < subgrid_list_.getCoarseNodeVector( coarse_index ).size(); ++coarse_node_local_id )
       {
         std::cout << coarse_node_local_id+1 << " : " << subgrid_list_.getCoarseNodeVector( coarse_index )[coarse_node_local_id] << std :: endl << std :: endl << std :: endl;
       }
*/
    
    //! boundary treatment:
    typedef typename LocProbFEMMatrix::LocalMatrixType LocalMatrix;

    typedef typename SGLagrangePointSetType::template Codim< faceCodim >::SubEntityIteratorType
    FaceDofIteratorType;

    const HostGridPartType& hostGridPart = hostDiscreteFunctionSpace_.gridPart();
      
    const SubgridIteratorType sg_end = subDiscreteFunctionSpace.end();
    for (SubgridIteratorType sg_it = subDiscreteFunctionSpace.begin(); sg_it != sg_end; ++sg_it)
    {
      const SubgridEntityType& subgrid_entity = *sg_it;

      HostEntityPointerType host_entity_pointer = subGrid.template getHostEntity< 0 >(subgrid_entity);
      const HostEntityType& host_entity = *host_entity_pointer;

      LocalMatrix local_matrix = locprob_system_matrix.localMatrix(subgrid_entity, subgrid_entity);

      const SGLagrangePointSetType& lagrangePointSet = subDiscreteFunctionSpace.lagrangePointSet(subgrid_entity);

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
        const FaceDofIteratorType fdend = lagrangePointSet.template endSubEntity< 1 >(face);
        for (FaceDofIteratorType fdit = lagrangePointSet.template beginSubEntity< 1 >(face); fdit != fdend; ++fdit)
          local_matrix.unitRow(*fdit);
      }
    }


    // assemble right hand side of algebraic local msfem problem
    if ( DSC_CONFIG_GET( "msfem.oversampling_strategy", 1 ) == 1 )
      { local_problem_op.assemble_local_RHS(e, local_problem_rhs); }
    else if ( ( DSC_CONFIG_GET( "msfem.oversampling_strategy", 1 ) == 2 ) ||
              ( DSC_CONFIG_GET( "msfem.oversampling_strategy", 1 ) == 3 ) )
      { if ( coarse_index < 0 )
          DUNE_THROW(Dune::InvalidStateException, "Invalid coarse index: coarse_index < 0");
	local_problem_op.assemble_local_RHS(e, subgrid_list_.getCoarseNodeVector( coarse_index ), local_problem_rhs ); }
    else
      DUNE_THROW(Dune::InvalidStateException, "Oversampling Strategy must be 1, 2 or 3!");
    // oneLinePrint( DSC_LOG_DEBUG, local_problem_rhs );

    // zero boundary condition for 'cell problems':
    // set Dirichlet Boundary to zero
    for (SubgridIteratorType sg_it = subDiscreteFunctionSpace.begin(); sg_it != sg_end; ++sg_it)
    {
      const SubgridEntityType& subgrid_entity = *sg_it;

      HostEntityPointerType host_entity_pointer = subGrid.template getHostEntity< 0 >(subgrid_entity);
      const HostEntityType& host_entity = *host_entity_pointer;

      HostIntersectionIterator iit = hostGridPart.ibegin(host_entity);
      const HostIntersectionIterator endiit = hostGridPart.iend(host_entity);
      for ( ; iit != endiit; ++iit)
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

        SubLocalFunctionType rhsLocal = local_problem_rhs.localFunction(subgrid_entity);
        const SGLagrangePointSetType& lagrangePointSet
          = subDiscreteFunctionSpace.lagrangePointSet(subgrid_entity);

        const int face = (*iit).indexInInside();

        FaceDofIteratorType faceIterator
          = lagrangePointSet.template beginSubEntity< faceCodim >(face);
        const FaceDofIteratorType faceEndIterator
          = lagrangePointSet.template endSubEntity< faceCodim >(face);
        for ( ; faceIterator != faceEndIterator; ++faceIterator)
          rhsLocal[*faceIterator] = 0;
      }
    }

    // After boundary treatment:
    // oneLinePrint( DSC_LOG_DEBUG, local_problem_rhs );

    const double norm_rhs = local_problem_op.normRHS(local_problem_rhs);

    if ( !( local_problem_rhs.dofsValid() ) )
    {
      DUNE_THROW(Dune::InvalidStateException, "Local MsFEM Problem RHS invalid.");
    }

    if (norm_rhs < /*1e-06*/ 1e-30)
    {
      local_problem_solution.clear();
      DSC_LOG_ERROR << "Local MsFEM problem with solution zero." << std::endl;
    } else {
      InverseLocProbFEMMatrix locprob_fem_biCGStab(locprob_system_matrix, 1e-8, 1e-8, 20000, LOCPROBLEMSOLVER_VERBOSE);
      locprob_fem_biCGStab(local_problem_rhs, local_problem_solution);
    }

    if ( !( local_problem_solution.dofsValid() ) )
    {
      DUNE_THROW(Dune::InvalidStateException,"Current solution of the local msfem problem invalid!");
    }

    // oneLinePrint( DSC_LOG_DEBUG, local_problem_solution );
  } // solvelocalproblem

  //! ----------- end method: solve local MsFEM problem ------------------------------------------

  // create a hostgrid function from a subgridfunction
  void subgrid_to_hostrid_function(const SubDiscreteFunctionType& sub_func,
                                   HostDiscreteFunctionType& host_func) {
    host_func.clear();

    const SubDiscreteFunctionSpaceType& subDiscreteFunctionSpace = sub_func.space();
    const SubGridType& subGrid = subDiscreteFunctionSpace.grid();

    SubgridIteratorType sub_endit = subDiscreteFunctionSpace.end();
    for (SubgridIteratorType sub_it = subDiscreteFunctionSpace.begin(); sub_it != sub_endit; ++sub_it)
    {
      const SubgridEntityType& sub_entity = *sub_it;

      HostEntityPointerType host_entity_pointer = subGrid.template getHostEntity< 0 >(*sub_it);
      const HostEntityType& host_entity = *host_entity_pointer;

      SubLocalFunctionType sub_loc_value = sub_func.localFunction(sub_entity);
      HostLocalFunctionType host_loc_value = host_func.localFunction(host_entity);

      const unsigned int numBaseFunctions = sub_loc_value.baseFunctionSet().size();
      for (unsigned int i = 0; i < numBaseFunctions; ++i)
      {
        host_loc_value[i] = sub_loc_value[i];
      }
    }
  } // subgrid_to_hostrid_function

  void output_local_solution(const int coarse_index, const int which,
                             const HostDiscreteFunctionType& host_local_solution) const
  {
    if (!DSC_CONFIG_GET("global.local_solution_vtk_output", false))
      return;
    typedef tuple< const HostDiscreteFunctionType* >      IOTupleType;
    typedef DataOutput< HostGridType, IOTupleType > DataOutputType;

    // general output parameters
    LocalProblemDataOutputParameters outputparam;
    // --------- data output local solution --------------

    // create and initialize output class
    IOTupleType local_solution_series(&host_local_solution);

    const std::string ls_name_s = (boost::format("/local_problem_solution_e%d_%d") % which % coarse_index).str();

    outputparam.set_prefix(ls_name_s);
    DataOutputType localsol_dataoutput(
      hostDiscreteFunctionSpace_.gridPart().grid(), local_solution_series, outputparam);
    localsol_dataoutput.writeData( 1.0 /*dummy*/, (boost::format("local-problem-solution-%d") % which).str() );
  }

  // method for solving and saving the solutions of the local msfem problems
  // for the whole set of macro-entities and for every unit vector e_i

  //! ---- method: solve and save the whole set of local msfem problems -----

  // Use the host-grid entities of Level 'computational_level' as computational domains for the subgrid computations
  void assemble_all(bool /*silent*/ = true /* state information on subgrids */) {
    std::string local_path = DSC_CONFIG_GET("global.datadir", "data") + "/local_problems/";
    Dune::Stuff::Common::testCreateDirectory(local_path);

    enum { dimension = MsfemTraits::GridType::dimension };
    enum { maxnumOfBaseFct = 100 };

    JacobianRangeType e[dimension];
    for (int i = 0; i < dimension; ++i)
      for (int j = 0; j < dimension; ++j)
      {
        if (i == j)
        { e[i][0][j] = 1.0; } else
        { e[i][0][j] = 0.0; }
      }

    // number of coarse grid entities (of codim 0).
    int number_of_coarse_grid_entities = specifier_.getNumOfCoarseEntities();

    DSC_LOG_INFO << "in method 'assemble_all': number_of_coarse_grid_entities = " << number_of_coarse_grid_entities
              << std::endl;
    DSC_PROFILER.startTiming("msfem.localproblemsolver.assemble_all");

    // we want to determine minimum, average and maxiumum time for solving a local msfem problem in the current method
    Dune::Stuff::Common::MinMaxAvg<double> cell_time;

    std::vector<int> coarse_indices;
    const HostDiscreteFunctionSpaceType& coarseSpace = specifier_.coarseSpace();
    const HostGridLeafIndexSet& coarseGridLeafIndexSet = coarseSpace.gridPart().grid().leafIndexSet();
    for (HostGridEntityIteratorType coarse_it = coarseSpace.begin(); coarse_it != coarseSpace.end(); ++coarse_it)
    {
      coarse_indices.push_back(coarseGridLeafIndexSet.index(*coarse_it));
    }


    const auto& comm = Dune::MPIHelper::getCollectiveCommunication();
    int slice = coarse_indices.size() / comm.size();
    for(int gc = comm.rank() * slice; gc < std::min(long(comm.rank() +1)* slice, long(coarse_indices.size())); ++gc)
    {
      const int coarse_index = coarse_indices[gc];

      DSC_LOG_INFO << "-------------------------" << std::endl
                   << "Coarse index " << coarse_index << std::endl;

      const std::string locprob_solution_location =
          (boost::format("local_problems/_localProblemSolutions_%d") % coarse_index).str();

      DiscreteFunctionWriter dfw(locprob_solution_location);

      SubGridType& subGrid = subgrid_list_.getSubGrid(coarse_index);

      SubGridPartType subGridPart(subGrid);

      DSC_LOG_INFO  << std::endl
                    << "Number of the local problem: " << dimension * coarse_index << " (of "
                    << (dimension * number_of_coarse_grid_entities) - 1 << " problems in total)" << std::endl
                    << "   Subgrid " << coarse_index << " contains " << subGrid.size(0) << " elements and "
                    << subGrid.size(2) << " nodes." << std::endl;

      const SubDiscreteFunctionSpaceType subDiscreteFunctionSpace(subGridPart);

      const std::string name_local_solution = (boost::format("Local Problem Solution %d") % coarse_index).str();

      //! only for dimension 2!
      SubDiscreteFunctionType local_problem_solution_0(name_local_solution, subDiscreteFunctionSpace);
      local_problem_solution_0.clear();

      SubDiscreteFunctionType local_problem_solution_1(name_local_solution, subDiscreteFunctionSpace);
      local_problem_solution_1.clear();

      // take time
      DSC_PROFILER.startTiming("none.local_problem_solution");

      // solve the problems
      solvelocalproblem(e[0], local_problem_solution_0, coarse_index);

      cell_time(DSC_PROFILER.stopTiming("none.local_problem_solution") / 1000.f);
      DSC_PROFILER.resetTiming("none.local_problem_solution");

      dfw.append(local_problem_solution_0);

      HostDiscreteFunctionType host_local_solution(name_local_solution, hostDiscreteFunctionSpace_);
      subgrid_to_hostrid_function(local_problem_solution_0, host_local_solution);
      output_local_solution(coarse_index, 0, host_local_solution);

      DSC_LOG_INFO  << std::endl
                    << "Number of the local problem: "
                    << (dimension * coarse_index) + 1 << " (of "
                    << (dimension * number_of_coarse_grid_entities) - 1 << " problems in total)" << std::endl
                    << "   Subgrid " << coarse_index << " contains " << subGrid.size(0) << " elements and "
                    << subGrid.size(2) << " nodes." << std::endl;

      // take time
      DSC_PROFILER.startTiming("none.local_problem_solution");

      // solve the problems
      solvelocalproblem(e[1], local_problem_solution_1, coarse_index);

      // min/max time
      cell_time(DSC_PROFILER.stopTiming("none.local_problem_solution") / 1000.f);
      DSC_PROFILER.resetTiming("none.local_problem_solution");

      dfw.append(local_problem_solution_1);

      subgrid_to_hostrid_function(local_problem_solution_1, host_local_solution);
      output_local_solution(coarse_index, 1, host_local_solution);
    } //for

    const auto total_time = DSC_PROFILER.stopTiming("msfem.localproblemsolver.assemble_all");
    DSC_LOG_INFO << std::endl;
    DSC_LOG_INFO << "In method: assemble_all." << std::endl << std::endl;
    DSC_LOG_INFO << "MsFEM problems solved for " << number_of_coarse_grid_entities << " coarse grid entities."
                  << std::endl;
    DSC_LOG_INFO << dimension * number_of_coarse_grid_entities << " local MsFEM problems solved in total."
                  << std::endl;
    DSC_LOG_INFO << "Minimum time for solving a local problem = " << cell_time.min() << "s." << std::endl;
    DSC_LOG_INFO << "Maximum time for solving a localproblem = " << cell_time.max() << "s." << std::endl;
    DSC_LOG_INFO << "Average time for solving a localproblem = "
                  << cell_time.average()<< "s." << std::endl;
    DSC_LOG_INFO << "Total time for computing and saving the localproblems = "
                    << total_time << "s," << std::endl << std::endl;
  } // assemble_all
}; // end class
} // end namespace Dune

#endif // #ifndef DiscreteEllipticMsFEMLocalProblem_HH
