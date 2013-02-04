#ifndef DiscreteEllipticMSFEMOperator_HH
#define DiscreteEllipticMSFEMOperator_HH

#include <dune/common/fmatrix.hh>

#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/multiscale/tools/solver/MsFEM/msfem_localproblems/subgrid-list.hh>
#include <dune/multiscale/tools/solver/MsFEM/msfem_localproblems/localproblemsolver.hh>
#include <dune/multiscale/tools/misc.hh>
#include <dune/fem/operator/2order/lagrangematrixsetup.hh>
#include <dune/stuff/fem/functions/checks.hh>

namespace Dune {

/**
 * \todo warum ist das nen operator?
 */
template< class CoarseDiscreteFunctionImp, class MacroMicroGridSpecifierImp, class FineDiscreteFunctionImp,
          class DiffusionImp >
class DiscreteEllipticMsFEMOperator
  : public Operator< typename CoarseDiscreteFunctionImp::RangeFieldType,
                     typename CoarseDiscreteFunctionImp::RangeFieldType,
                     CoarseDiscreteFunctionImp, CoarseDiscreteFunctionImp >
    , boost::noncopyable
{
  typedef DiscreteEllipticMsFEMOperator< CoarseDiscreteFunctionImp, MacroMicroGridSpecifierImp, FineDiscreteFunctionImp,
                                         DiffusionImp > This;

private:
  typedef CoarseDiscreteFunctionImp  CoarseDiscreteFunction;
  typedef FineDiscreteFunctionImp    FineDiscreteFunction;
  typedef MacroMicroGridSpecifierImp MacroMicroGridSpecifierType;

  typedef DiffusionImp DiffusionModel;

  typedef typename CoarseDiscreteFunction::DiscreteFunctionSpaceType CoarseDiscreteFunctionSpace;
  typedef typename FineDiscreteFunction::DiscreteFunctionSpaceType   FineDiscreteFunctionSpace;

  typedef typename FineDiscreteFunctionSpace::FunctionSpaceType FunctionSpace;

  typedef typename FineDiscreteFunctionSpace::GridPartType FineGridPart;
  typedef typename FineDiscreteFunctionSpace::GridType     FineGrid;

  typedef typename FineDiscreteFunctionSpace::RangeFieldType RangeFieldType;
  typedef typename FineDiscreteFunctionSpace::DomainType     DomainType;
  typedef typename FineDiscreteFunctionSpace::RangeType      RangeType;
  typedef typename FineDiscreteFunctionSpace::JacobianRangeType
  JacobianRangeType;

  typedef SubGrid< FineGrid::dimensionworld, FineGrid >  SubGridType;
  typedef SubGridList< FineDiscreteFunction, SubGridType, MacroMicroGridSpecifierType > SubGridListType;

  typedef MsFEMLocalProblemSolver< FineDiscreteFunction, SubGridListType, MacroMicroGridSpecifierType,
                                   DiffusionModel > MsFEMLocalProblemSolverType;

  static const int dimension = FineGridPart::GridType::dimension;
  static const int polynomialOrder = FineDiscreteFunctionSpace::polynomialOrder;

  typedef typename FineDiscreteFunction::LocalFunctionType FineLocalFunction;
  typedef typename FineDiscreteFunctionSpace::BaseFunctionSetType                   FineBaseFunctionSet;
  typedef typename FineDiscreteFunctionSpace::LagrangePointSetType                  FineLagrangePointSet;
  typedef typename FineLagrangePointSet::template Codim< 1 >::SubEntityIteratorType FineFaceDofIterator;
  typedef typename FineGrid::Traits::LeafIndexSet FineGridLeafIndexSet;
  typedef typename FineDiscreteFunctionSpace::IteratorType FineIterator;
  typedef typename FineIterator::Entity                    FineEntity;
  typedef typename FineEntity::EntityPointer               FineEntityPointer;
  typedef typename FineEntity::Geometry                    FineGeometry;
  typedef typename FineGridPart::IntersectionIteratorType FineIntersectionIterator;
  typedef typename FineIntersectionIterator::Intersection FineIntersection;
  typedef CachingQuadrature< FineGridPart, 0 > FineQuadrature;

  typedef typename CoarseDiscreteFunctionSpace::GridPartType CoarseGridPart;
  typedef typename CoarseDiscreteFunctionSpace::GridType     CoarseGrid;

  typedef typename CoarseDiscreteFunction::LocalFunctionType CoarseLocalFunction;
  typedef typename CoarseDiscreteFunctionSpace::BaseFunctionSetType                   CoarseBaseFunctionSet;
  typedef typename CoarseDiscreteFunctionSpace::LagrangePointSetType                  CoarseLagrangePointSet;
  typedef typename CoarseLagrangePointSet::template Codim< 1 >::SubEntityIteratorType CoarseFaceDofIterator;
  typedef typename CoarseDiscreteFunctionSpace::IteratorType CoarseIterator;
  typedef typename CoarseGrid::Traits::LeafIndexSet          CoarseGridLeafIndexSet;
  typedef typename CoarseIterator::Entity CoarseEntity;
  typedef typename CoarseEntity::Geometry CoarseGeometry;
  typedef typename CoarseGridPart::IntersectionIteratorType CoarseIntersectionIterator;
  typedef typename CoarseIntersectionIterator::Intersection CoarseIntersection;
  typedef CachingQuadrature< CoarseGridPart, 0 > CoarseQuadrature;

  //! ---------------- typedefs for the SubgridDiscreteFunctionSpace -----------------------
  // ( typedefs for the local grid and the corresponding local ('sub') )discrete space )

  //! type of grid part
  typedef LeafGridPart< SubGridType > SubGridPart;

  //! type of subgrid discrete function space
  typedef LagrangeDiscreteFunctionSpace< FunctionSpace, SubGridPart, 1 >  // 1=POLORDER
  LocalDiscreteFunctionSpace;

  //! type of subgrid discrete function
  typedef AdaptiveDiscreteFunction< LocalDiscreteFunctionSpace > LocalDiscreteFunction;
  typedef typename LocalDiscreteFunctionSpace::IteratorType LocalGridIterator;
  typedef typename LocalGridIterator::Entity LocalGridEntity;
  typedef typename LocalGridEntity::EntityPointer LocalGridEntityPointer;
  typedef typename LocalDiscreteFunction::LocalFunctionType LocalGridLocalFunction;
  typedef typename LocalDiscreteFunctionSpace::LagrangePointSetType LGLagrangePointSet;
  typedef typename LocalDiscreteFunctionSpace::BaseFunctionSetType LocalGridBaseFunctionSet;
  typedef typename LocalGridEntity::Geometry LocalGridGeometry;
  typedef CachingQuadrature< SubGridPart, 0 > LocalGridQuadrature;

  //!-----------------------------------------------------------------------------------------

public:
  DiscreteEllipticMsFEMOperator(MacroMicroGridSpecifierType& specifier,
                                const CoarseDiscreteFunctionSpace& coarseDiscreteFunctionSpace,
                                // number of layers per coarse grid entity T:  U(T) is created by enrichting T with
                                // n(T)-layers:
                                SubGridListType& subgrid_list,
                                const DiffusionModel& diffusion_op)
    : specifier_(specifier)
      , coarseDiscreteFunctionSpace_(coarseDiscreteFunctionSpace)
      , subgrid_list_(subgrid_list)
      , diffusion_operator_(diffusion_op) {
    bool silence = false;

    // coarseDiscreteFunctionSpace_ = specifier_.coarseSpace();
    // fineDiscreteFunctionSpace_ = specifier_.fineSpace();
    MsFEMLocalProblemSolverType loc_prob_solver(
      specifier_.fineSpace(), specifier_, subgrid_list_, diffusion_operator_);

    loc_prob_solver.assemble_all(silence);
  }

  template< class MatrixType >
  void assemble_matrix(MatrixType& global_matrix) const;

private:
  //! dummy operator
  virtual void operator()(const CoarseDiscreteFunction& u, CoarseDiscreteFunction& w) const;

  // oneLinePrint( std :: cout , );
  template< class Stream >
  void oneLinePrint(Stream& stream, const LocalDiscreteFunction& func) {
    typedef typename LocalDiscreteFunction::ConstDofIteratorType
    DofIteratorType;
    DofIteratorType it = func.dbegin();
    stream << "\n" << func.name() << ": [ ";
    for ( ; it != func.dend(); ++it)
      stream << std::setw(5) << *it << "  ";

    stream << " ] " << std::endl;
  } // oneLinePrint

  /** create a hostgrid function from a subgridfunction
   * Note: the maximum gride levels for both underlying grids must be the same
   * \todo use function from estimator utils instead
   **/
  void subgrid_to_hostrid_function(const LocalDiscreteFunction& sub_func, FineDiscreteFunction& host_func);

  MacroMicroGridSpecifierType& specifier_;

  const CoarseDiscreteFunctionSpace& coarseDiscreteFunctionSpace_;
  SubGridListType& subgrid_list_;
  const DiffusionModel& diffusion_operator_;
};

// create a hostgrid function from a subgridfunction
template< class CoarseDiscreteFunctionImp, class MacroMicroGridSpecifierImp, class FineDiscreteFunctionImp,
          class DiffusionImp >
void DiscreteEllipticMsFEMOperator< CoarseDiscreteFunctionImp,
                                    MacroMicroGridSpecifierImp,
                                    FineDiscreteFunctionImp,
                                    DiffusionImp >::subgrid_to_hostrid_function(const LocalDiscreteFunction& sub_func,
                                                                                FineDiscreteFunction& host_func) {
  if ( sub_func.space().gridPart().grid().maxLevel() != host_func.space().gridPart().grid().maxLevel() )
  {
    DSC_LOG_ERROR
    << "Error in method 'subgrid_to_hostrid_function': MaxLevel of SubGrid not identical to MaxLevel of FineGrid."
    << std::endl;
  }

  host_func.clear();

  const LocalDiscreteFunctionSpace& subDiscreteFunctionSpace = sub_func.space();
  const SubGridType& subGrid = subDiscreteFunctionSpace.grid();

  LocalGridIterator sub_endit = subDiscreteFunctionSpace.end();
  for (LocalGridIterator sub_it = subDiscreteFunctionSpace.begin(); sub_it != sub_endit; ++sub_it)
  {
    const LocalGridEntity& sub_entity = *sub_it;

    FineEntityPointer host_entity_pointer = subGrid.template getHostEntity< 0 >(*sub_it);
    const FineEntity& host_entity = *host_entity_pointer;

    LocalGridLocalFunction sub_loc_value = sub_func.localFunction(sub_entity);
    FineLocalFunction host_loc_value = host_func.localFunction(host_entity);

    const auto numBaseFunctions = sub_loc_value.baseFunctionSet().numBaseFunctions();
    for (unsigned int i = 0; i < numBaseFunctions; ++i)
    {
      host_loc_value[i] = sub_loc_value[i];
    }
  }
} // subgrid_to_hostrid_function

// dummy implementation of "operator()"
// 'w' = effect of the discrete operator on 'u'
template< class K, class L, class O, class P >
void DiscreteEllipticMsFEMOperator< K,L,O,P >::operator()(const CoarseDiscreteFunction& /*u*/,
                                                               CoarseDiscreteFunction& /*w*/) const {
  DUNE_THROW(Dune::NotImplemented,"the ()-operator of the DiscreteEllipticMsFEMOperator class is not yet implemented and still a dummy.");
}

template< class CoarseDiscreteFunctionImp, class MacroMicroGridSpecifierImp, class FineDiscreteFunctionImp,
          class DiffusionImp >
template< class MatrixType >
void DiscreteEllipticMsFEMOperator< CoarseDiscreteFunctionImp,
                                    MacroMicroGridSpecifierImp,
                                    FineDiscreteFunctionImp,
                                    DiffusionImp >::assemble_matrix(MatrixType& global_matrix) const {
  // the local problem:
  // Let 'T' denote a coarse grid element and
  // let 'U(T)' denote the environment of 'T' that corresponds with the subgrid.

  // if Petrov-Galerkin-MsFEM
  if ( DSC_CONFIG_GET("msfem.petrov_galerkin", 1 ) )
    DSC_LOG_INFO << "Assembling Petrov-Galerkin-MsFEM Matrix." << std::endl;
  else  // if classical (symmetric) MsFEM
    DSC_LOG_INFO << "Assembling MsFEM Matrix." << std::endl;


  typedef typename MatrixType::LocalMatrixType LocalMatrix;

  global_matrix.reserve();
  global_matrix.clear();

  std::vector< typename CoarseBaseFunctionSet::JacobianRangeType > gradient_Phi(
    coarseDiscreteFunctionSpace_.mapper().maxNumDofs() );

  const auto& coarseGridLeafIndexSet = coarseDiscreteFunctionSpace_.gridPart().grid().leafIndexSet();

  for (const CoarseEntity& coarse_grid_entity : coarseDiscreteFunctionSpace_)
  {
    const CoarseGeometry& coarse_grid_geometry = coarse_grid_entity.geometry();
    assert(coarse_grid_entity.partitionType() == InteriorEntity);

    const int global_index_entity = coarseGridLeafIndexSet.index(coarse_grid_entity);

    LocalMatrix local_matrix = global_matrix.localMatrix(coarse_grid_entity, coarse_grid_entity);

    const CoarseBaseFunctionSet& coarse_grid_baseSet = local_matrix.domainBaseFunctionSet();
    const unsigned int numMacroBaseFunctions = coarse_grid_baseSet.size();

    // the sub grid U(T) that belongs to the coarse_grid_entity T
    SubGridType& sub_grid_U_T = subgrid_list_.getSubGrid(global_index_entity);
    SubGridPart subGridPart(sub_grid_U_T);

    LocalDiscreteFunctionSpace localDiscreteFunctionSpace(subGridPart);

    LocalDiscreteFunction local_problem_solution_e0("Local problem Solution e_0", localDiscreteFunctionSpace);
    local_problem_solution_e0.clear();

    LocalDiscreteFunction local_problem_solution_e1("Local problem Solution e_1", localDiscreteFunctionSpace);
    local_problem_solution_e1.clear();

    // --------- load local solutions -------
    // the file/place, where we saved the solutions of the cell problems
    const std::string local_solution_location = (boost::format("local_problems/_localProblemSolutions_%d")
                                                % global_index_entity).str();
    // reader for the cell problem data file:
    DiscreteFunctionReader discrete_function_reader(local_solution_location);
    discrete_function_reader.read(0, local_problem_solution_e0);
    discrete_function_reader.read(1, local_problem_solution_e1);

    // 1 point quadrature!! We only need the gradient of the base function,
    // which is constant on the whole entity.
    const CoarseQuadrature one_point_quadrature(coarse_grid_entity, 0);

    // the barycenter of the macro_grid_entity
    const typename CoarseQuadrature::CoordinateType& local_coarse_point
      = one_point_quadrature.point(0 /*=quadraturePoint*/);

    // transposed of the the inverse jacobian
    const auto& inverse_jac = coarse_grid_geometry.jacobianInverseTransposed(local_coarse_point);
    coarse_grid_baseSet.jacobianAll(one_point_quadrature[0], inverse_jac, gradient_Phi);

    for (unsigned int i = 0; i < numMacroBaseFunctions; ++i)
    {
      for (unsigned int j = 0; j < numMacroBaseFunctions; ++j)
      {
        RangeType local_integral = 0.0;

        // iterator for the micro grid ( grid for the reference element T_0 )
        const LocalGridIterator local_grid_end = localDiscreteFunctionSpace.end();
        for (LocalGridIterator local_grid_it = localDiscreteFunctionSpace.begin();
             local_grid_it != local_grid_end;
             ++local_grid_it)
        {
          const LocalGridEntity& local_grid_entity = *local_grid_it;

          // check if "local_grid_entity" (which is an entity of U(T)) is in T:
          // -------------------------------------------------------------------

          auto father_of_loc_grid_ent =
              Stuff::Grid::make_father(coarseGridLeafIndexSet,
                                   localDiscreteFunctionSpace.grid().template getHostEntity< 0 >(local_grid_entity),
                                   specifier_.getLevelDifference());
          if (!Stuff::Grid::entities_identical(coarse_grid_entity, *father_of_loc_grid_ent))
          {
            continue;
          }
          // -------------------------------------------------------------------

          const LocalGridGeometry& local_grid_geometry = local_grid_entity.geometry();
          assert(local_grid_entity.partitionType() == InteriorEntity);

          // higher order quadrature, since A^{\epsilon} is highly variable
          LocalGridQuadrature local_grid_quadrature(local_grid_entity, 2 * localDiscreteFunctionSpace.order() + 2);
          const size_t numQuadraturePoints = local_grid_quadrature.nop();

          for (size_t localQuadraturePoint = 0; localQuadraturePoint < numQuadraturePoints; ++localQuadraturePoint)
          {
            // local (barycentric) coordinates (with respect to entity)
            const typename LocalGridQuadrature::CoordinateType& local_subgrid_point = local_grid_quadrature.point(
              localQuadraturePoint);

            DomainType global_point_in_U_T = local_grid_geometry.global(local_subgrid_point);

            const double weight_local_quadrature
              = local_grid_quadrature.weight(localQuadraturePoint) * local_grid_geometry.integrationElement(
              local_subgrid_point);

            LocalGridLocalFunction localized_local_problem_solution_e0 = local_problem_solution_e0.localFunction(
              local_grid_entity);
            LocalGridLocalFunction localized_local_problem_solution_e1 = local_problem_solution_e1.localFunction(
              local_grid_entity);

            // grad coorector for e_0 and e_1
            typename LocalGridBaseFunctionSet::JacobianRangeType grad_loc_sol_e0, grad_loc_sol_e1;
            localized_local_problem_solution_e0.jacobian(local_grid_quadrature[localQuadraturePoint], grad_loc_sol_e0);
            localized_local_problem_solution_e1.jacobian(local_grid_quadrature[localQuadraturePoint], grad_loc_sol_e1);

            // ∇ Phi_H + ∇ Q( Phi_H ) = ∇ Phi_H + ∂_x1 Phi_H Q( e_1 ) + ∂_x2 Phi_H Q( e_2 )
            JacobianRangeType direction_of_diffusion(0.0);
            for (int k = 0; k < dimension; ++k)
            {
              direction_of_diffusion[0][k] += gradient_Phi[i][0][0] * grad_loc_sol_e0[0][k];
              direction_of_diffusion[0][k] += gradient_Phi[i][0][1] * grad_loc_sol_e1[0][k];
              direction_of_diffusion[0][k] += gradient_Phi[i][0][k];
            }

            JacobianRangeType diffusive_flux(0.0);
            diffusion_operator_.diffusiveFlux(global_point_in_U_T, direction_of_diffusion, diffusive_flux);

            // if not Petrov-Galerkin MsFEM:
            if ( !DSC_CONFIG_GET("msfem.petrov_galerkin", true ) )
	    {
              JacobianRangeType reconstruction_grad_phi_j(0.0);
              for (int k = 0; k < dimension; ++k)
              {
               reconstruction_grad_phi_j[0][k] += gradient_Phi[j][0][0] * grad_loc_sol_e0[0][k];
               reconstruction_grad_phi_j[0][k] += gradient_Phi[j][0][1] * grad_loc_sol_e1[0][k];
               reconstruction_grad_phi_j[0][k] += gradient_Phi[j][0][k];
              }

              local_integral += weight_local_quadrature * (diffusive_flux[0] * reconstruction_grad_phi_j[0]);
	    }
            else // if Petrov-Galerkin MsFEM
            {
              local_integral += weight_local_quadrature * (diffusive_flux[0] * gradient_Phi[j][0]);
	    }
          }
        }

        // add entries
        local_matrix.add(j, i, local_integral);
      }
    }
  }

  // boundary treatment
  //!TODO use function call
  const CoarseGridPart& coarseGridPart = coarseDiscreteFunctionSpace_.gridPart();
  for (CoarseIterator it = coarseDiscreteFunctionSpace_.begin(); it != coarseDiscreteFunctionSpace_.end(); ++it)
  {
    const CoarseEntity& entity = *it;
    if ( !entity.hasBoundaryIntersections() )
      continue;

    LocalMatrix local_matrix = global_matrix.localMatrix(entity, entity);

    const CoarseLagrangePointSet& lagrangePointSet = coarseDiscreteFunctionSpace_.lagrangePointSet(entity);

    const CoarseIntersectionIterator iend = coarseGridPart.iend(entity);
    for (CoarseIntersectionIterator iit = coarseGridPart.ibegin(entity); iit != iend; ++iit)
    {
      const CoarseIntersection& intersection = *iit;
      if ( !intersection.boundary() )
        continue;

      const int face = intersection.indexInInside();
      const CoarseFaceDofIterator fdend = lagrangePointSet.template endSubEntity< 1 >(face);
      for (CoarseFaceDofIterator fdit = lagrangePointSet.template beginSubEntity< 1 >(face); fdit != fdend; ++fdit)
        local_matrix.unitRow(*fdit);
    }
  }
} // assemble_matrix
}

#endif // #ifndef DiscreteElliptic_HH
