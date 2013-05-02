// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DiscreteEllipticMsFEMConservativeFluxSolver_HH
#define DiscreteEllipticMsFEMConservativeFluxSolver_HH

#include <vector>
#include <dune/common/fmatrix.hh>

#include <dune/fem/operator/matrix/spmatrix.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/geometry/quadraturerules.hh>

// FLUX_SOLVER_VERBOSE: 0 = false, 1 = true
#define FLUX_SOLVER_VERBOSE false

// VTK output for conservative flux solution
// #define VTK_OUTPUT

// dune-subgrid include:
#include <dune/multiscale/tools/subgrid_io.hh>
#include <dune/subgrid/subgrid.hh>

// dune-fem includes:
#include <dune/fem/gridpart/common/gridpart.hh>
#include <dune/fem/operator/2order/lagrangematrixsetup.hh>
#include <dune/fem/solver/inverseoperators.hh>
#include <dune/stuff/common/math.hh>

#include <dune/multiscale/tools/misc.hh>
#include <dune/multiscale/tools/misc/outputparameter.hh>
#include <dune/multiscale/hmm/cell_problem_numbering.hh>
#include <dune/stuff/common/profiler.hh>
#include <dune/stuff/fem/localmatrix_proxy.hh>

namespace Dune {
namespace Multiscale {
namespace MsFEM {

/** \brief define output parameters for local problems
 *  appends "local_problems" for path
 **/
struct ConFluxProblemDataOutputParameters
  : public OutputParameters
{
public:
  explicit ConFluxProblemDataOutputParameters()
    : OutputParameters(DSC_CONFIG_GET("global.datadir", "data") + "/cf_problems/")
  {}
};

//! \TODO docme
template< class SubGridDiscreteFunctionType, class DiscreteFunctionType, class DiffusionOperatorType,
          class MacroMicroGridSpecifierType >
class ConservativeFluxOperator
  : public Operator< typename SubGridDiscreteFunctionType::RangeFieldType,
                     typename SubGridDiscreteFunctionType::RangeFieldType,
                     SubGridDiscreteFunctionType,
                     SubGridDiscreteFunctionType >
{
  typedef ConservativeFluxOperator< SubGridDiscreteFunctionType, DiscreteFunctionType, DiffusionOperatorType,
                                    MacroMicroGridSpecifierType > This;

private:
  typedef SubGridDiscreteFunctionType SubGridDiscreteFunction;
  typedef DiscreteFunctionType        DiscreteFunction;
  typedef DiffusionOperatorType       DiffusionModel;

  typedef typename SubGridDiscreteFunction::DiscreteFunctionSpaceType SubGridDiscreteFunctionSpace;
  typedef typename DiscreteFunction::DiscreteFunctionSpaceType        DiscreteFunctionSpace;

  typedef typename DiscreteFunctionSpace::GridPartType GridPart;
  typedef typename DiscreteFunctionSpace::GridType     GridType;

  typedef typename SubGridDiscreteFunctionSpace::GridPartType SubGridPart;
  typedef typename SubGridDiscreteFunctionSpace::GridType     SubGridType;

  typedef typename DiscreteFunctionSpace::RangeFieldType RangeFieldType;
  typedef typename DiscreteFunctionSpace::DomainType     DomainType;
  typedef typename DiscreteFunctionSpace::RangeType      RangeType;
  typedef typename DiscreteFunctionSpace::JacobianRangeType
    JacobianRangeType;

  static const int dimension = GridPart::GridType::dimension;
  static const int polynomialOrder = DiscreteFunctionSpace::polynomialOrder;

  typedef typename DiscreteFunction::LocalFunctionType        LocalFunction;
  typedef typename SubGridDiscreteFunction::LocalFunctionType SubGridLocalFunction;

  typedef typename DiscreteFunctionSpace::BaseFunctionSetType                   BaseFunctionSet;
  typedef typename DiscreteFunctionSpace::LagrangePointSetType                  LagrangePointSet;
  typedef typename LagrangePointSet::template Codim< 1 >::SubEntityIteratorType FaceDofIterator;

  typedef typename DiscreteFunctionSpace::IteratorType Iterator;
  typedef typename Iterator::Entity                    Entity;
  typedef typename Entity::EntityPointer               EntityPointer;
  typedef typename Entity::Geometry                    Geometry;

  typedef typename GridPart::IntersectionIteratorType IntersectionIterator;
  typedef typename IntersectionIterator::Intersection Intersection;
  typedef typename Intersection::LocalCoordinate      LocalCoordinate;

  typedef typename SubGridDiscreteFunctionSpace::BaseFunctionSetType                   SubGridBaseFunctionSet;
  typedef typename SubGridDiscreteFunctionSpace::LagrangePointSetType                  SubGridLagrangePointSet;
  typedef typename SubGridLagrangePointSet::template Codim< 1 >::SubEntityIteratorType SubGridFaceDofIterator;

  typedef typename SubGridDiscreteFunctionSpace::IteratorType SubGridIterator;
  typedef typename SubGridIterator::Entity                    SubGridEntity;
  typedef typename SubGridEntity::Geometry                    SubGridGeometry;

  typedef typename SubGridPart::IntersectionIteratorType     SubGridIntersectionIterator;
  typedef typename SubGridIntersectionIterator::Intersection SubGridIntersection;

  typedef CachingQuadrature< GridPart, 0 > Quadrature;

  typedef QuadratureRule< double, 1 >      FaceQuadratureRule;
  typedef CachingQuadrature< GridPart, 1 > FaceQuadrature;

  typedef CachingQuadrature< SubGridPart, 0 > SubGridQuadrature;

  typedef typename GridType::template Codim< 1 >::Geometry FaceGeometryType;

public:
  ConservativeFluxOperator(const SubGridDiscreteFunctionSpace& subDiscreteFunctionSpace,
                           const DiscreteFunctionSpace& discreteFunctionSpace,
                           const DiffusionModel& diffusion_op,
                           const MacroMicroGridSpecifierType& specifier)
    : subDiscreteFunctionSpace_(subDiscreteFunctionSpace)
      , discreteFunctionSpace_(discreteFunctionSpace)
      , diffusion_operator_(diffusion_op)
      , specifier_(specifier)
  {}

private:
  ConservativeFluxOperator(const This&) = delete;

public:
  // dummy operator
  virtual void operator()(const SubGridDiscreteFunction& u, SubGridDiscreteFunction& w) const;

  // assemble stiffness matrix for local problems
  template< class MatrixType >
  void assemble_matrix(const int sub_grid_id, MatrixType& global_matrix) const;

  // the right hand side assembler methods
  void assemble_RHS(  // direction 'e_i'
    JacobianRangeType& e_i,
    // solution of the local corrector problem
    const SubGridDiscreteFunction& local_corrector_e_i,
    const int sub_grid_id,
    // rhs local msfem problem:
    SubGridDiscreteFunction& rhs_flux_problem) const;

  void printLocalRHS(SubGridDiscreteFunction& rhs) const;

  double normRHS(SubGridDiscreteFunction& rhs) const;

private:
  const SubGridDiscreteFunctionSpace& subDiscreteFunctionSpace_;
  const DiscreteFunctionSpace& discreteFunctionSpace_;
  const DiffusionModel& diffusion_operator_;
  const MacroMicroGridSpecifierType& specifier_;
};

//! dummy implementation of "operator()"
//! 'w' = effect of the discrete operator on 'u'
template< class SubGridDiscreteFunctionImp, class DiscreteFunctionImp, class DiffusionImp,
          class MacroMicroGridSpecifierImp >
void ConservativeFluxOperator< SubGridDiscreteFunctionImp, DiscreteFunctionImp, DiffusionImp,
                               MacroMicroGridSpecifierImp >
  ::operator()(const SubGridDiscreteFunctionImp& /*u*/, SubGridDiscreteFunctionImp& /*w*/) const {
  DUNE_THROW(Dune::NotImplemented,"the ()-operator of the ConservativeFluxOperator class is not yet implemented and still a dummy.");
}

//! assemble system matrix
template< class SubGridDiscreteFunctionImp, class DiscreteFunctionImp, class DiffusionImp,
          class MacroMicroGridSpecifierImp >
template< class MatrixType >
void ConservativeFluxOperator< SubGridDiscreteFunctionImp, DiscreteFunctionImp, DiffusionImp,
                               MacroMicroGridSpecifierImp >
  ::assemble_matrix(const int sub_grid_id, MatrixType& global_matrix) const {
  typedef typename MatrixType::LocalMatrixType LocalMatrix;

  global_matrix.reserve();
  global_matrix.clear();

  // local grid basis functions:
  std::vector< RangeType > phi( subDiscreteFunctionSpace_.mapper().maxNumDofs() );

  const SubGridIterator end = subDiscreteFunctionSpace_.end();
  for (SubGridIterator it = subDiscreteFunctionSpace_.begin(); it != end; ++it)
  {
    const SubGridEntity& sub_grid_entity = *it;

    EntityPointer host_entity_pointer = subDiscreteFunctionSpace_.gridPart().grid().template getHostEntity< 0 >(
      sub_grid_entity);

    typedef typename GridType::Traits::LeafIndexSet CoarseGridLeafIndexSet;
    const CoarseGridLeafIndexSet& coarseGridLeafIndexSet = specifier_.coarseSpace().gridPart().grid().leafIndexSet();

    EntityPointer father_of_sub_grid_entity = Stuff::Grid::make_father(coarseGridLeafIndexSet,
                                                                       host_entity_pointer,
                                                                       specifier_.getLevelDifference());
    const int coarse_index = coarseGridLeafIndexSet.index(*father_of_sub_grid_entity);
    assert(sub_grid_entity.partitionType() == InteriorEntity);

    DSFe::LocalMatrixProxy<MatrixType> local_matrix(global_matrix, sub_grid_entity, sub_grid_entity);

    const SubGridBaseFunctionSet& baseSet = local_matrix.domainBaseFunctionSet();
    const auto numBaseFunctions = baseSet.size();

    const IntersectionIterator iend = discreteFunctionSpace_.gridPart().iend(*host_entity_pointer);
    for (IntersectionIterator iit = discreteFunctionSpace_.gridPart().ibegin(*host_entity_pointer); iit != iend; ++iit)
    {
      FaceQuadrature faceQuadrature(discreteFunctionSpace_.gridPart(), *iit,
                                    2 * subDiscreteFunctionSpace_.order() + 1, FaceQuadrature::INSIDE);

      const FaceGeometryType& faceGeometry = iit->geometry();

      bool set_zero = false;
      if (coarse_index != sub_grid_id)
      {
        set_zero = true;
      }

      if ( ( iit->neighbor() ) && (!set_zero) )
      {
        EntityPointer outside_it = iit->outside();

        EntityPointer father_of_neighbor = Stuff::Grid::make_father(coarseGridLeafIndexSet,
                                                                    outside_it, specifier_.getLevelDifference());
        if (Stuff::Grid::entities_identical(*father_of_sub_grid_entity,*father_of_neighbor))
        {
          set_zero = true;
        }
      }

      const size_t numQuadraturePoints = faceQuadrature.nop();

      RangeType check_sum(0.0);

      for (size_t quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint)
      {
        // das liefert immer den gleichen 'local_point' egal fuer welchen Quadraturpunkt.
        const LocalCoordinate local_point = faceGeometry.local( faceQuadrature.point(quadraturePoint) );

        // integration factors
        const double integrationFactor = faceGeometry.integrationElement(local_point);

        // weight
        const double quadratureWeight = faceQuadrature.weight(quadraturePoint);

        check_sum += integrationFactor * quadratureWeight;

        baseSet.evaluateAll(faceQuadrature[quadraturePoint], phi);

        for (unsigned int i = 0; i < numBaseFunctions; ++i)
        {
          for (unsigned int j = 0; j < numBaseFunctions; ++j)
          {
            if (!set_zero)
            {
              local_matrix.add( j, i, integrationFactor * quadratureWeight * (phi[i][0] * phi[j][0]) );
            } else {
              // stabilization (should be close to zero):
              local_matrix.add( j, i, /*0.0*/ 0.00000001 * integrationFactor * quadratureWeight
                                * (phi[i][0] * phi[j][0]) );
            }
          }
        }
      } // done loop over all quadrature points

      if ( check_sum != faceGeometry.volume() )
      {
        DUNE_THROW(Dune::InvalidStateException,"Error in Face Quadrature.");
      }
    }
  }
} // assemble_matrix

template< class SubGridDiscreteFunctionImp, class DiscreteFunctionImp, class DiffusionImp,
          class MacroMicroGridSpecifierImp >
void ConservativeFluxOperator< SubGridDiscreteFunctionImp, DiscreteFunctionImp, DiffusionImp,
                               MacroMicroGridSpecifierImp >
  ::printLocalRHS(SubGridDiscreteFunctionImp& rhs) const {
  typedef typename SubGridDiscreteFunctionImp::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionSpaceType::IteratorType               IteratorType;
  typedef typename DiscreteFunctionImp::LocalFunctionType                LocalFunctionType;

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

template< class SubGridDiscreteFunctionImp, class DiscreteFunctionImp, class DiffusionImp,
          class MacroMicroGridSpecifierImp >
double ConservativeFluxOperator< SubGridDiscreteFunctionImp, DiscreteFunctionImp, DiffusionImp,
                                 MacroMicroGridSpecifierImp >
  ::normRHS(SubGridDiscreteFunctionImp& rhs) const {
  double norm = 0.0;

  typedef typename SubGridDiscreteFunctionImp::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionSpaceType::IteratorType               IteratorType;
  typedef typename IteratorType::Entity                                  EntityType;
  typedef typename SubGridDiscreteFunctionImp::LocalFunctionType         LocalFunctionType;
  typedef typename DiscreteFunctionSpaceType::GridPartType               GridPartType;
  typedef typename DiscreteFunctionSpaceType::GridType                   GridType;
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

// assemble the right hand side of the conservative flux problem
// ----------------------------------------------

template< class SubGridDiscreteFunctionImp, class DiscreteFunctionImp, class DiffusionImp,
          class MacroMicroGridSpecifierImp >
// template< class MatrixType >
void ConservativeFluxOperator< SubGridDiscreteFunctionImp, DiscreteFunctionImp, DiffusionImp,
                               MacroMicroGridSpecifierImp >
  ::assemble_RHS
  ( // direction 'e'
  JacobianRangeType& e_i,
  // solution of the local corrector problem
  const SubGridDiscreteFunction& local_corrector_e_i,
  const int sub_grid_id,
  // rhs flux problem:
  SubGridDiscreteFunction& rhs_flux_problem) const {
  const SubGridDiscreteFunctionSpace& subDiscreteFunctionSpace = rhs_flux_problem.space();

  // set entries to zero:
  rhs_flux_problem.clear();

  // gradient of micro scale base function:
  std::vector< JacobianRangeType > gradient_phi( subDiscreteFunctionSpace.mapper().maxNumDofs() );

  const SubGridIterator end = subDiscreteFunctionSpace.end();
  for (SubGridIterator it = subDiscreteFunctionSpace.begin(); it != end; ++it)
  {
    const SubGridEntity& local_grid_entity = *it;

    EntityPointer host_entity_pointer = subDiscreteFunctionSpace.gridPart().grid().template getHostEntity< 0 >(
      local_grid_entity);

    typedef typename GridType::Traits::LeafIndexSet CoarseGridLeafIndexSet;
    const CoarseGridLeafIndexSet& coarseGridLeafIndexSet = specifier_.coarseSpace().gridPart().grid().leafIndexSet();

    EntityPointer father_of_sub_grid_entity = Stuff::Grid::make_father(coarseGridLeafIndexSet,
                                                                       host_entity_pointer,
                                                                       specifier_.getLevelDifference());
    const int coarse_index = coarseGridLeafIndexSet.index(*father_of_sub_grid_entity);

    if (coarse_index != sub_grid_id)
    {
      continue;
    }

    const SubGridGeometry& geometry = local_grid_entity.geometry();
    assert(local_grid_entity.partitionType() == InteriorEntity);

    SubGridLocalFunction elementOfRHS = rhs_flux_problem.localFunction(local_grid_entity);

    const SubGridBaseFunctionSet& baseSet = elementOfRHS.baseFunctionSet();
    const auto numBaseFunctions = baseSet.size();

    SubGridQuadrature quadrature(local_grid_entity, 2 * subDiscreteFunctionSpace.order() + 2);
    const size_t numQuadraturePoints = quadrature.nop();
    for (size_t quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint)
    {
      const typename SubGridQuadrature::CoordinateType& local_point = quadrature.point(quadraturePoint);

      // remember, we are concerned with: - \int_{U(T)} (A^eps)(x) e · ∇ \phi(x)

      // global point in the subgrid
      const DomainType global_point = geometry.global(local_point);

      const double weight = quadrature.weight(quadraturePoint) * geometry.integrationElement(local_point);

      // transposed of the the inverse jacobian
      const auto& inverse_jac = geometry.jacobianInverseTransposed(local_point);

      // A^eps(x) ( e
      // diffusion operator evaluated in 'x' multiplied with e
      JacobianRangeType diffusion_in_e_i;
      diffusion_operator_.diffusiveFlux(global_point, e_i, diffusion_in_e_i);

      JacobianRangeType total_diffusive_flux;

      JacobianRangeType grad_corrector_e_i;
      SubGridLocalFunction localized_corrector_e_i = local_corrector_e_i.localFunction(local_grid_entity);
      localized_corrector_e_i.jacobian(quadrature[quadraturePoint], grad_corrector_e_i);
      diffusion_operator_.diffusiveFlux(global_point, grad_corrector_e_i, total_diffusive_flux);

      total_diffusive_flux[0] += diffusion_in_e_i[0];

      baseSet.jacobianAll(quadrature[quadraturePoint], inverse_jac, gradient_phi);

      for (unsigned int i = 0; i < numBaseFunctions; ++i)
      {
        elementOfRHS[i] -= weight * (total_diffusive_flux[0] * gradient_phi[i][0]);
      }
    }
  } // end for-loop of 'Iterator'
} // assemble_RHS


//! the essential local msfem problem solver class
template< class SubGridDiscreteFunctionType,
          class HostDiscreteFunctionType,
          class DiffusionOperatorType,
          class MacroMicroGridSpecifierImp >
class ConservativeFluxProblemSolver
{
private:
  typedef MacroMicroGridSpecifierImp MacroMicroGridSpecifierType;

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

  typedef typename HostGridType::Traits::LeafIndexSet HostGridLeafIndexSet;

  typedef typename HostDiscreteFunctionSpaceType::IteratorType HostGridEntityIteratorType;

  typedef typename HostGridEntityIteratorType::Entity HostEntityType;

  typedef typename HostEntityType::EntityPointer HostEntityPointerType;

  typedef typename HostGridType::template Codim< 0 >::Geometry HostGridEntityGeometry;

  typedef typename HostDiscreteFunctionType::LocalFunctionType HostLocalFunctionType;

  typedef typename HostGridPartType::IntersectionIteratorType HostIntersectionIterator;

  typedef typename HostGridType::Traits::LeafIndexSet LeafIndexSetType;

  static const int dimension = HostGridType::dimension;

  //! ---------------- typedefs for the SubGridDiscreteFunctionSpace -----------------------
  // ( typedefs for the local grid and the corresponding local ('sub') )discrete space )

  //! type of discrete function space
  typedef typename SubGridDiscreteFunctionType::DiscreteFunctionSpaceType
  SubGridDiscreteFunctionSpaceType;

  typedef typename SubGridDiscreteFunctionSpaceType::RangeFieldType RangeFieldType;
  typedef typename SubGridDiscreteFunctionSpaceType::DomainType     DomainType;
  typedef typename SubGridDiscreteFunctionSpaceType::RangeType      RangeType;
  typedef typename SubGridDiscreteFunctionSpaceType::JacobianRangeType
  JacobianRangeType;

  //! type of grid partition
  typedef typename SubGridDiscreteFunctionSpaceType::GridPartType SubGridPartType;

  //! type of grid
  typedef typename SubGridDiscreteFunctionSpaceType::GridType SubGridType;

  typedef typename SubGridDiscreteFunctionSpaceType::IteratorType SubGridIteratorType;

  typedef typename SubGridIteratorType::Entity SubGridEntityType;

  typedef typename SubGridEntityType::EntityPointer SubGridEntityPointerType;

  typedef typename SubGridDiscreteFunctionType::LocalFunctionType SubGridLocalFunctionType;

  typedef typename SubGridDiscreteFunctionSpaceType::LagrangePointSetType SubGridLagrangePointSetType;

  //!-----------------------------------------------------------------------------------------

  //! ------------------ Matrix Traits for the local Problems ---------------------

  enum { faceCodim = 1 };
  typedef typename SubGridLagrangePointSetType::template Codim< faceCodim >::SubEntityIteratorType
  SubGridFaceDofIteratorType;

  //! polynomial order of base functions
  enum { polynomialOrder = SubGridDiscreteFunctionSpaceType::polynomialOrder };

  // flux problem matrix traits
  struct FluxProbMatrixTraits
  {
    typedef SubGridDiscreteFunctionSpaceType                          RowSpaceType;
    typedef SubGridDiscreteFunctionSpaceType                          ColumnSpaceType;
    typedef LagrangeMatrixSetup< false >                              StencilType;
    typedef ParallelScalarProduct< SubGridDiscreteFunctionSpaceType > ParallelScalarProductType;

    template< class M >
    struct Adapter
    {
      typedef LagrangeParallelMatrixAdapter< M > MatrixAdapterType;
    };
  };

  typedef SparseRowMatrixOperator< SubGridDiscreteFunctionType, SubGridDiscreteFunctionType,
                                   FluxProbMatrixTraits > FluxProbFEMMatrix;

  // OEMGMRESOp //OEMBICGSQOp // OEMBICGSTABOp /*CGInverseOp*/
  typedef CGInverseOperator< SubGridDiscreteFunctionType, FluxProbFEMMatrix > InverseFluxProbFEMMatrix;

private:
  const DiffusionOperatorType& diffusion_;
  const HostDiscreteFunctionSpaceType& hostDiscreteFunctionSpace_;

  const MacroMicroGridSpecifierType& specifier_;
  mutable boost::format filename_template_;

public:
  //! constructor - with diffusion operator A^{\epsilon}(x)
  ConservativeFluxProblemSolver(const HostDiscreteFunctionSpaceType& hostDiscreteFunctionSpace,
                                const DiffusionOperatorType& diffusion_operator,
                                const MacroMicroGridSpecifierType& specifier)
    : diffusion_(diffusion_operator)
      , hostDiscreteFunctionSpace_(hostDiscreteFunctionSpace)
      , specifier_(specifier)
      , filename_template_(boost::format("/_conservativeFlux_e_%d_sg_%d"))
  {}

  template< class Stream >
  void oneLinePrint(Stream& stream, const SubGridDiscreteFunctionType& func) const {
    typedef typename SubGridDiscreteFunctionType::ConstDofIteratorType
    DofIteratorType;
    DofIteratorType it = func.dbegin();
    stream << "\n" << func.name() << ": [ ";
    for ( ; it != func.dend(); ++it)
      stream << std::setw(5) << *it << "  ";

    stream << " ] " << std::endl;
  } // oneLinePrint

  //! ----------- method: solve the local MsFEM problem ------------------------------------------

  void solve(JacobianRangeType& e_i,  // direction 'e_i'
             const SubGridDiscreteFunctionType& local_corrector_e_i,
             const int sub_grid_id,
             const int direction_index,
             SubGridDiscreteFunctionType& conservative_flux) const {
    // set solution equal to zero:
    conservative_flux.clear();

    const SubGridDiscreteFunctionSpaceType& localDiscreteFunctionSpace = local_corrector_e_i.space();

    //! the matrix in our linear system of equations
    FluxProbFEMMatrix flux_prob_system_matrix("Conservative Flux Problem System Matrix",
                                              localDiscreteFunctionSpace,
                                              localDiscreteFunctionSpace);

    //! define the discrete (elliptic) local MsFEM problem operator
    // ( effect of the discretized differential operator on a certain discrete function )
    // discrete elliptic operator describing the elliptic local msfem problems
    typedef ConservativeFluxOperator< SubGridDiscreteFunctionType,
                                      HostDiscreteFunctionType,
                                      DiffusionOperatorType,
                                      MacroMicroGridSpecifierType > ConservativeFluxOperatorType;
    ConservativeFluxOperatorType cf_problem_operator(localDiscreteFunctionSpace,
                                                     hostDiscreteFunctionSpace_,
                                                     diffusion_,
                                                     specifier_);

    //! right hand side vector of the algebraic local MsFEM problem
    SubGridDiscreteFunctionType rhs("RHS of Conservative Flux Problem", localDiscreteFunctionSpace);
    rhs.clear();

    // assemble the stiffness matrix
    cf_problem_operator.assemble_matrix(sub_grid_id, flux_prob_system_matrix);

    // assemble right hand side of algebraic local msfem problem
    cf_problem_operator.assemble_RHS(e_i, local_corrector_e_i, sub_grid_id, rhs);

    // oneLinePrint( DSC_LOG_DEBUG, rhs );

    const double norm_rhs = cf_problem_operator.normRHS(rhs);

    if ( !( rhs.dofsValid() ) )
    {
      DUNE_THROW(Dune::InvalidStateException,"Local Flux Problem RHS invalid.");
    }

    if (norm_rhs < /*1e-06*/ 1e-30)
    {
      conservative_flux.clear();
      DSC_LOG_INFO << "Local Flux Problem with solution zero." << std::endl;
    } else {
      InverseFluxProbFEMMatrix flux_prob_biCGStab(flux_prob_system_matrix, 1e-8, 1e-8, 20000, FLUX_SOLVER_VERBOSE);
      flux_prob_biCGStab(rhs, conservative_flux);
    }

    if ( !( conservative_flux.dofsValid() ) )
    {
      DUNE_THROW(Dune::InvalidStateException,"Solution of the Local Flux Problem is invalid!");
    }

    #ifdef VTK_OUTPUT
    vtk_output(conservative_flux, sub_grid_id, direction_index);
    #endif
    file_data_output(conservative_flux, sub_grid_id, direction_index);
  } // solve

  //! ----------- end method: solve local MsFEM problem ------------------------------------------

  // create a hostgrid function from a subgridfunction
  void subgrid_to_hostrid_function(const SubGridDiscreteFunctionType& sub_func,
                                   HostDiscreteFunctionType& host_func) const {
    host_func.clear();

    const SubGridDiscreteFunctionSpaceType& subDiscreteFunctionSpace = sub_func.space();
    const SubGridType& subGrid = subDiscreteFunctionSpace.grid();

    SubGridIteratorType sub_endit = subDiscreteFunctionSpace.end();
    for (SubGridIteratorType sub_it = subDiscreteFunctionSpace.begin(); sub_it != sub_endit; ++sub_it)
    {
      const SubGridEntityType& sub_entity = *sub_it;

      HostEntityPointerType host_entity_pointer = subGrid.template getHostEntity< 0 >(*sub_it);
      const HostEntityType& host_entity = *host_entity_pointer;

      SubGridLocalFunctionType sub_loc_value = sub_func.localFunction(sub_entity);
      HostLocalFunctionType host_loc_value = host_func.localFunction(host_entity);

      const auto numBaseFunctions = sub_loc_value.baseFunctionSet().size();
      for (unsigned int i = 0; i < numBaseFunctions; ++i)
      {
        host_loc_value[i] = sub_loc_value[i];
      }
    }
  } // subgrid_to_hostrid_function

  void vtk_output(const SubGridDiscreteFunctionType& subgrid_disc_func,
                  const int sub_grid_index,
                  const int direction_index) const {
    // --------------- writing data output ---------------------
    // typedefs and initialization
    typedef Dune::tuple< HostDiscreteFunctionType* >      IOTupleType;
    typedef DataOutput< HostGridType, IOTupleType > DataOutputType;

    // general output parameters
    ConFluxProblemDataOutputParameters outputparam;
    // -------------------------------------------------------

    HostDiscreteFunctionType host_disc_func("Conservative Flux", hostDiscreteFunctionSpace_);
    subgrid_to_hostrid_function(subgrid_disc_func, host_disc_func);

    // create and initialize output class
    IOTupleType conservative_flux_series(&host_disc_func);
    filename_template_ % direction_index % sub_grid_index;
    outputparam.set_prefix(filename_template_.str());
    DataOutputType cf_dataoutput(hostDiscreteFunctionSpace_.gridPart().grid(), conservative_flux_series, outputparam);

    cf_dataoutput.writeData( 1.0 /*dummy*/, "conservative-flux" );
  } // vtk_output

  void file_data_output(const SubGridDiscreteFunctionType& subgrid_disc_func,
                        const int sub_grid_index,
                        const int direction_index) const {
    const std::string locprob_solution_location
        = std::string("cf_problems/")
          + (filename_template_ % direction_index % sub_grid_index).str();
    DiscreteFunctionWriter(locprob_solution_location).append(subgrid_disc_func);
  } // file_data_output

  template< typename SubGridListType >
  void solve_all(SubGridListType& subgrid_list) const {
    JacobianRangeType e[dimension];

    for (int i = 0; i < dimension; ++i) {
      for (int j = 0; j < dimension; ++j)
      {
        if (i == j)
        { e[i][0][j] = 1.0; } else
        { e[i][0][j] = 0.0; }
      }
    }

    // number of coarse grid entities (of codim 0).
    const int number_of_coarse_grid_entities = specifier_.getNumOfCoarseEntities();

    DSC_PROFILER.startTiming("msfem.conservative_flux_solver.solve_all_subgrids");

    // we want to determine minimum, average and maxiumum time for solving a local msfem problem in the current method
    Dune::Stuff::Common::MinMaxAvg<double> cell_time;


    const HostDiscreteFunctionSpaceType& coarseSpace = specifier_.coarseSpace();
    const LeafIndexSetType& coarseGridLeafIndexSet = coarseSpace.gridPart().grid().leafIndexSet();

    // Coarse Entity Iterator
    const HostGridEntityIteratorType coarse_grid_end = coarseSpace.end();
    for (HostGridEntityIteratorType coarse_grid_it = coarseSpace.begin();
         coarse_grid_it != coarse_grid_end;
         ++coarse_grid_it)
    {
      int global_index_entity = coarseGridLeafIndexSet.index(*coarse_grid_it);

      // the sub grid U(T) that belongs to the coarse_grid_entity T
      SubGridType& sub_grid_U_T = subgrid_list.getSubGrid(global_index_entity);
      SubGridPartType subGridPart(sub_grid_U_T);

      SubGridDiscreteFunctionSpaceType localDiscreteFunctionSpace(subGridPart);

      SubGridDiscreteFunctionType local_problem_solution_e0("Local problem Solution e_0", localDiscreteFunctionSpace);
      local_problem_solution_e0.clear();

      SubGridDiscreteFunctionType local_problem_solution_e1("Local problem Solution e_1", localDiscreteFunctionSpace);
      local_problem_solution_e1.clear();

      // --------- load local solutions -------
      // the file/place, where we saved the solutions of the cell problems
      const std::string local_solution_location = (boost::format("local_problems/_localProblemSolutions_%d")
                                            % global_index_entity).str();

      // reader for the cell problem data file:
      DiscreteFunctionReader discrete_function_reader(local_solution_location);
      discrete_function_reader.read(0, local_problem_solution_e0);
      discrete_function_reader.read(1, local_problem_solution_e1);

      SubGridDiscreteFunctionType conservative_flux_e0("Conservative Flux for e_0", localDiscreteFunctionSpace);
      SubGridDiscreteFunctionType conservative_flux_e1("Conservative Flux for e_1", localDiscreteFunctionSpace);

      DSC_LOG_INFO << "Number of the 'conservative flux problem': " << dimension * global_index_entity << " (of "
                << (dimension * number_of_coarse_grid_entities) - 1 << " problems in total)" << std::endl;

      // take time
      DSC_PROFILER.startTiming("none.local_problem_solution");

      this->solve(e[0], local_problem_solution_e0, global_index_entity, 0, conservative_flux_e0);

      DSC_LOG_INFO << "Number of the 'conservative flux problem': "
                << (dimension * global_index_entity) + 1 << " (of "
                << (dimension * number_of_coarse_grid_entities) - 1 << " problems in total)" << std::endl;

      cell_time(DSC_PROFILER.stopTiming("none.local_problem_solution") / 1000.f);
      DSC_PROFILER.resetTiming("none.local_problem_solution");
      DSC_PROFILER.startTiming("none.local_problem_solution");

      this->solve(e[1], local_problem_solution_e1, global_index_entity, 1, conservative_flux_e1);

      cell_time(DSC_PROFILER.stopTiming("none.local_problem_solution") / 1000.f);
      DSC_PROFILER.resetTiming("none.local_problem_solution");
    }
    const auto total_time = DSC_PROFILER.stopTiming("msfem.conservative_flux_solver.solve_all_subgrids") / 1000.f;
    DSC_LOG_INFO << std::endl;
    DSC_LOG_INFO << "In: 'assemble all conservatice fluxes'." << std::endl << std::endl;
    DSC_LOG_INFO << "Conservative Flux determined for " << number_of_coarse_grid_entities
                  << " coarse grid entities." << std::endl;
    DSC_LOG_INFO << dimension * number_of_coarse_grid_entities
                  << " conservative flux problems solved in total." << std::endl;
    DSC_LOG_INFO << "Minimum time for solving a conservative flux problem = " << cell_time.min() << "s."
                  << std::endl;
    DSC_LOG_INFO << "Maximum time for solving a conservative flux problem = " << cell_time.max() << "s."
                  << std::endl;
    DSC_LOG_INFO << "Average time for solving a conservative flux problem = "
                  << cell_time.average() << "s." << std::endl;
    DSC_LOG_INFO << "Total time for computing and saving the conservative flux problems = "
                  << total_time << "s," << std::endl << std::endl;

  } // solve_all
}; // end class

} //namespace MsFEM {
} //namespace Multiscale {
} //namespace Dune {

#endif // #ifndef DiscreteEllipticMsFEMLocalProblem_HH
