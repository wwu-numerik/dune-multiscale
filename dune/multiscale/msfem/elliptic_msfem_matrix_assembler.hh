// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef MSFEM_ELLIPTIC_DiscreteEllipticMSFEMOperator_HH
#define MSFEM_ELLIPTIC_DiscreteEllipticMSFEMOperator_HH

#include <config.h>
#include <type_traits>

#include <dune/common/fmatrix.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/operator/2order/lagrangematrixsetup.hh>

#include <dune/multiscale/common/dirichletconstraints.hh>
#include <dune/multiscale/msfem/localproblems/subgrid-list.hh>
#include <dune/multiscale/msfem/localproblems/localproblemsolver.hh>
#include <dune/multiscale/msfem/localproblems/localsolutionmanager.hh>
#include <dune/multiscale/msfem/msfem_grid_specifier.hh>
#include <dune/multiscale/msfem/msfem_traits.hh>
#include <dune/multiscale/tools/misc.hh>
#include <dune/multiscale/tools/discretefunctionwriter.hh>
#include <dune/multiscale/hmm/cell_problem_numbering.hh>
#include <dune/multiscale/problems/base.hh>
#include <dune/multiscale/problems/selector.hh>

#include <dune/stuff/fem/functions/checks.hh>
#include <dune/stuff/fem/localmatrix_proxy.hh>
#include <dune/stuff/fem/matrix_object.hh>

namespace Dune {
namespace Multiscale {
namespace MsFEM {
/**
 * \todo docme
 */
class DiscreteEllipticMsFEMOperator
  : boost::noncopyable
{
private:
  typedef CommonTraits::DiscreteFunctionType CoarseDiscreteFunction;
  typedef CommonTraits::DiscreteFunctionType FineDiscreteFunction;
  typedef MacroMicroGridSpecifier MacroMicroGridSpecifierType;

  typedef CommonTraits::DiffusionType                                DiffusionModel;

  typedef typename CoarseDiscreteFunction::DiscreteFunctionSpaceType CoarseDiscreteFunctionSpace;
  typedef typename FineDiscreteFunction::DiscreteFunctionSpaceType   FineDiscreteFunctionSpace;

  typedef typename FineDiscreteFunctionSpace::FunctionSpaceType      FunctionSpace;

  typedef typename FineDiscreteFunctionSpace::GridPartType           FineGridPart;
  typedef typename FineDiscreteFunctionSpace::GridType               FineGrid;

  typedef typename FineDiscreteFunctionSpace::RangeFieldType         RangeFieldType;
  typedef typename FineDiscreteFunctionSpace::DomainType             DomainType;
  typedef typename FineDiscreteFunctionSpace::RangeType              RangeType;
  typedef typename FineDiscreteFunctionSpace::JacobianRangeType      JacobianRangeType;

  typedef MsFEMLocalProblemSolver                                    MsFEMLocalProblemSolverType;

  static const int dimension       = FineGridPart::GridType::dimension;
  static const int polynomialOrder = FineDiscreteFunctionSpace::polynomialOrder;

  typedef typename FineDiscreteFunction::LocalFunctionType                   FineLocalFunction;
  typedef typename FineDiscreteFunctionSpace::BasisFunctionSetType           FineBaseFunctionSet;
  typedef typename FineDiscreteFunctionSpace::LagrangePointSetType           FineLagrangePointSet;
  typedef typename FineLagrangePointSet::Codim< 1 >::SubEntityIteratorType   FineFaceDofIterator;
  typedef typename FineGrid::Traits::LeafIndexSet                            FineGridLeafIndexSet;
  typedef typename FineDiscreteFunctionSpace::IteratorType                   FineIterator;
  typedef typename FineIterator::Entity                                      FineEntity;
  typedef typename FineEntity::EntityPointer                                 FineEntityPointer;
  typedef typename FineEntity::Geometry                                      FineGeometry;
  typedef typename FineGridPart::IntersectionIteratorType                    FineIntersectionIterator;
  typedef typename FineIntersectionIterator::Intersection                    FineIntersection;
  typedef Fem::CachingQuadrature< FineGridPart, 0 >                          FineQuadratureType;

  typedef typename CoarseDiscreteFunctionSpace::GridPartType                 CoarseGridPart;
  typedef typename CoarseDiscreteFunctionSpace::GridType                     CoarseGrid;

  typedef typename CoarseDiscreteFunction::LocalFunctionType                 CoarseLocalFunction;
  typedef typename CoarseDiscreteFunctionSpace::BasisFunctionSetType         CoarseBaseFunctionSet;
  typedef typename CoarseDiscreteFunctionSpace::LagrangePointSetType         CoarseLagrangePointSet;
  typedef typename CoarseLagrangePointSet::Codim< 1 >::SubEntityIteratorType CoarseFaceDofIterator;
  typedef typename CoarseDiscreteFunctionSpace::IteratorType                 CoarseIterator;
  typedef typename CoarseGrid::Traits::LeafIndexSet                          CoarseGridLeafIndexSet;
  typedef typename CoarseIterator::Entity                                    CoarseEntity;
  typedef typename CoarseEntity::Geometry                                    CoarseGeometry;
  typedef typename CoarseGridPart::IntersectionIteratorType                  CoarseIntersectionIterator;
  typedef typename CoarseIntersectionIterator::Intersection                  CoarseIntersection;
  typedef Fem::CachingQuadrature< CoarseGridPart, 0 >                        CoarseQuadrature;

  typedef MsFEMTraits::SubGridListType SubGridListType;
  typedef typename SubGridListType::SubGridQuadratureType SubGridQuadratureType;

public:
  DiscreteEllipticMsFEMOperator(MacroMicroGridSpecifierType& specifier,
                                const CoarseDiscreteFunctionSpace& coarseDiscreteFunctionSpace,
                                // number of layers per coarse grid entity T:  U(T) is created by enrichting T with
                                // n(T)-layers:
                                MsFEMTraits::SubGridListType& subgrid_list,
                                const DiffusionModel& diffusion_op);

  template< class SPMatrixObject >
  void assemble_matrix(SPMatrixObject& global_matrix) const;

private:
  MacroMicroGridSpecifierType&       specifier_;
  const CoarseDiscreteFunctionSpace& coarseDiscreteFunctionSpace_;
  MsFEMTraits::SubGridListType&      subgrid_list_;
  const DiffusionModel&              diffusion_operator_;
  const bool                         petrovGalerkin_;
};

template< class SPMatrixObject >
void DiscreteEllipticMsFEMOperator::assemble_matrix(SPMatrixObject& global_matrix) const
{
  // the local problem:
  // Let 'T' denote a coarse grid element and
  // let 'U(T)' denote the environment of 'T' that corresponds with the subgrid.

  // if Petrov-Galerkin-MsFEM
  if (petrovGalerkin_)
    DSC_LOG_INFO << "Assembling Petrov-Galerkin-MsFEM Matrix." << std::endl;
  else  // if classical (symmetric) MsFEM
    DSC_LOG_INFO << "Assembling MsFEM Matrix." << std::endl;

  global_matrix.reserve(DSFe::diagonalAndNeighborStencil(global_matrix));
  global_matrix.clear();

  const auto& coarseGridLeafIndexSet = coarseDiscreteFunctionSpace_.gridPart().grid().leafIndexSet();

  for (const CoarseEntity& coarse_grid_entity : coarseDiscreteFunctionSpace_) {
    const CoarseGeometry& coarse_grid_geometry = coarse_grid_entity.geometry();
    assert(coarse_grid_entity.partitionType() == InteriorEntity);

    const int global_index_entity = coarseGridLeafIndexSet.index(coarse_grid_entity);

    DSFe::LocalMatrixProxy< SPMatrixObject > local_matrix(global_matrix, coarse_grid_entity, coarse_grid_entity);

    const CoarseBaseFunctionSet& coarse_grid_baseSet   = local_matrix.domainBasisFunctionSet();
    const unsigned int           numMacroBaseFunctions = coarse_grid_baseSet.size();

    // Load local solutions
    Multiscale::MsFEM::LocalSolutionManager localSolutionManager(coarse_grid_entity, subgrid_list_, specifier_);
    localSolutionManager.loadLocalSolutions();
    Multiscale::MsFEM::LocalSolutionManager::LocalSolutionVectorType& localSolutions
      = localSolutionManager.getLocalSolutions();
    assert(localSolutions.size() > 0);
    std::vector< typename CoarseBaseFunctionSet::JacobianRangeType > gradientPhi(numMacroBaseFunctions);

    const int localQuadratureOrder = 2 * localSolutionManager.getLocalDiscreteFunctionSpace().order() + 2;

    // iterator for the micro grid ( grid for the reference element T_0 )
    for (const auto& localGridEntity : localSolutionManager.getLocalDiscreteFunctionSpace()) {
      // check if "localGridEntity" (which is an entity of U(T)) is in T:
      // -------------------------------------------------------------------
      const auto& hostEntity
        = localSolutionManager.getSubGridPart().grid().template getHostEntity< 0 >(localGridEntity);
      // ignore overlay elements
      if (global_index_entity == subgrid_list_.getEnclosingMacroCellIndex(hostEntity)) {
        assert(hostEntity->partitionType() == InteriorEntity);


        const auto& local_grid_geometry = localGridEntity.geometry();

        // higher order quadrature, since A^{\epsilon} is highly variable
        SubGridQuadratureType localQuadrature(localGridEntity, localQuadratureOrder);
        const size_t       numQuadraturePoints = localQuadrature.nop();

        // number of local solutions without the boundary correctors. Those are only needed for the right hand side
        const int numLocalSolutions = localSolutions.size() - localSolutionManager.numBoundaryCorrectors();
        // evaluate the jacobians of all local solutions in all quadrature points
        std::vector< std::vector< JacobianRangeType > >
        allLocalSolutionEvaluations(numLocalSolutions,
                                    std::vector< JacobianRangeType >(localQuadrature.nop(), JacobianRangeType(0.0)));
        for (int lsNum = 0; lsNum < numLocalSolutions; ++lsNum) {
          auto localFunction = localSolutions[lsNum]->localFunction(localGridEntity);
          localFunction.evaluateQuadrature(localQuadrature, allLocalSolutionEvaluations[lsNum]);
        }

        for (size_t localQuadraturePoint = 0; localQuadraturePoint < numQuadraturePoints; ++localQuadraturePoint) {
          // local (barycentric) coordinates (with respect to entity)
          const typename FineQuadratureType::CoordinateType& local_subgrid_point = localQuadrature.point(
            localQuadraturePoint);

          DomainType   global_point_in_U_T = local_grid_geometry.global(local_subgrid_point);
          const double weight_local_quadrature
            = localQuadrature.weight(localQuadraturePoint) * local_grid_geometry.integrationElement(
            local_subgrid_point);

          // evaluate the jacobian of the coarse grid base set
          const DomainType& local_coarse_point = coarse_grid_geometry.local(global_point_in_U_T);
          coarse_grid_baseSet.jacobianAll(local_coarse_point, gradientPhi);


          for (unsigned int i = 0; i < numMacroBaseFunctions; ++i) {
            for (unsigned int j = 0; j < numMacroBaseFunctions; ++j) {
              RangeType local_integral(0.0);

              // Compute the gradients of the i'th and j'th local problem solutions
              JacobianRangeType gradLocProbSoli(0.0), gradLocProbSolj(0.0);
              if (specifier_.simplexCoarseGrid()) {
                assert(allLocalSolutionEvaluations.size()==dimension);
                // ∇ Phi_H + ∇ Q( Phi_H ) = ∇ Phi_H + ∂_x1 Phi_H ∇Q( e_1 ) + ∂_x2 Phi_H ∇Q( e_2 )
                for (int k = 0; k < dimension; ++k) {
                  gradLocProbSoli.axpy(gradientPhi[i][0][k], allLocalSolutionEvaluations[k][localQuadraturePoint]);
                  gradLocProbSolj.axpy(gradientPhi[j][0][k], allLocalSolutionEvaluations[k][localQuadraturePoint]);
                }
              } else {
                assert(allLocalSolutionEvaluations.size()==numMacroBaseFunctions);
                gradLocProbSoli = allLocalSolutionEvaluations[i][localQuadraturePoint];
                gradLocProbSolj = allLocalSolutionEvaluations[j][localQuadraturePoint];
              }

              JacobianRangeType reconstructionGradPhii(gradientPhi[i]);
              reconstructionGradPhii += gradLocProbSoli;
              JacobianRangeType reconstructionGradPhij(gradientPhi[j]);
              reconstructionGradPhij += gradLocProbSolj;
              JacobianRangeType diffusive_flux(0.0);
              diffusion_operator_.diffusiveFlux(global_point_in_U_T, reconstructionGradPhii, diffusive_flux);
              if (petrovGalerkin_)
                local_integral += weight_local_quadrature * (diffusive_flux[0] * gradientPhi[j][0]);
              else
                local_integral += weight_local_quadrature * (diffusive_flux[0] * reconstructionGradPhij[0]);

              // add entries
              local_matrix.add(j, i, local_integral);
            }
          }
        }
      }
    }
  }

  const auto boundary = Problem::getModelData()->boundaryInfo();
  Dune::DirichletConstraints<CoarseDiscreteFunctionSpace> constraints(*boundary, coarseDiscreteFunctionSpace_);
  constraints.applyToOperator(global_matrix);
} // assemble_matrix
} // namespace MsFEM {
} // namespace Multiscale {
} // namespace Dune {

#endif // #ifndef MSFEM_ELLIPTIC_DiscreteElliptic_HH
