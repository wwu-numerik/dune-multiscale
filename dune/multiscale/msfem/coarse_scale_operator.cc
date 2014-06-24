#include <config.h>
#include <dune/stuff/common/parameter/configcontainer.hh>
#include <dune/stuff/fem/functions/integrals.hh>
#include <dune/stuff/common/profiler.hh>
#include <dune/gdt/operators/projections.hh>
#include <dune/gdt/operators/prolongations.hh>
#include <dune/gdt/spaces/constraints.hh>
#include <dune/gdt/functionals/l2.hh>
#include <dune/multiscale/msfem/diffusion_evaluation.hh>
#include <sstream>

#include "dune/multiscale/msfem/msfem_traits.hh"
#include "coarse_scale_operator.hh"

namespace Dune {
namespace Multiscale {
namespace MsFEM {

class MsFEMCodim0Integral;

class MsFEMCodim0IntegralTraits
{
public:
  typedef MsFEMCodim0Integral derived_type;
};

class MsFEMCodim0Integral
  : public GDT::LocalOperator::Codim0Interface< MsFEMCodim0IntegralTraits >
{
public:
  typedef MsFEMCodim0IntegralTraits Traits;

private:
  static const size_t numTmpObjectsRequired_ = 1;

public:
  explicit MsFEMCodim0Integral(const size_t over_integrate = 0)
    : over_integrate_(over_integrate)
  {}


  size_t numTmpObjectsRequired() const
  {
    return numTmpObjectsRequired_;
  }

  template< class E, class D, int d, class R, int rT, int rCT, int rA, int rCA >
  void apply(Multiscale::MsFEM::LocalSolutionManager& localSolutionManager,
             const MsFEMTraits::LocalEntityType& localGridEntity,
             const Stuff::LocalfunctionSetInterface< E, D, d, R, rT, rCT >& testBase,
             const Stuff::LocalfunctionSetInterface< E, D, d, R, rA, rCA >& ansatzBase,
             Dune::DynamicMatrix< R >& ret,
             std::vector< Dune::DynamicMatrix< R > >& tmpLocalMatrices) const
  {
    auto& diffusion_operator = DMP::getDiffusion();
    const auto& coarse_scale_entity = ansatzBase.entity();

    // quadrature
    typedef Dune::QuadratureRules< D, d > VolumeQuadratureRules;
    typedef Dune::QuadratureRule< D, d > VolumeQuadratureType;
    const size_t integrand_order = diffusion_operator->order() + ansatzBase.order() + testBase.order() + over_integrate_;
    assert(integrand_order < std::numeric_limits< int >::max());
    const VolumeQuadratureType& volumeQuadrature = VolumeQuadratureRules::rule(coarse_scale_entity.type(), int(integrand_order));
    // check matrix and tmp storage
    const size_t rows = testBase.size();
    const size_t cols = ansatzBase.size();
    Dune::Stuff::Common::clear(ret);
    assert(ret.rows() >= rows);
    assert(ret.cols() >= cols);
    assert(tmpLocalMatrices.size() >= numTmpObjectsRequired_);


    const auto numQuadraturePoints = volumeQuadrature.size();
    const auto& localSolutions = localSolutionManager.getLocalSolutions();
    // number of local solutions without the boundary correctors. Those are only needed for the right hand side
    const auto numLocalSolutions = localSolutions.size() - localSolutionManager.numBoundaryCorrectors();
    typedef CommonTraits::DiscreteFunctionSpaceType::BaseFunctionSetType::RangeType RangeType;
    typedef CommonTraits::DiscreteFunctionSpaceType::BaseFunctionSetType::JacobianRangeType JacobianRangeType;
    // evaluate the jacobians of all local solutions in all quadrature points
    std::vector<std::vector<JacobianRangeType>> allLocalSolutionEvaluations(
          numLocalSolutions, std::vector<JacobianRangeType>(numQuadraturePoints, RangeType(0.0)));
    for (auto lsNum : DSC::valueRange(numLocalSolutions)) {
      const auto localFunction = localSolutions[lsNum]->local_function(localGridEntity);
//      assert(localSolutionManager.space().indexSet().contains(localGridEntity));
      localFunction->jacobian(volumeQuadrature, allLocalSolutionEvaluations[lsNum]);
    }

    // loop over all quadrature points
    const auto quadPointEndIt = volumeQuadrature.end();
    std::size_t localQuadraturePoint = 0;
    for (auto quadPointIt = volumeQuadrature.begin(); quadPointIt != quadPointEndIt; ++quadPointIt, ++localQuadraturePoint) {
      const Dune::FieldVector< D, d > x = quadPointIt->position();
      const auto coarseBaseJacs_ = testBase.jacobian(x);
      // integration factors
      const double integrationFactor = coarse_scale_entity.geometry().integrationElement(x);
      const double quadratureWeight = quadPointIt->weight();
      const auto global_point_in_U_T = coarse_scale_entity.geometry().global(x);

      // compute integral
      for (size_t ii = 0; ii < rows; ++ii) {
        auto& retRow = ret[ii];
        for (size_t jj = 0; jj < cols; ++jj) {
          RangeType local_integral(0.0);

          // Compute the gradients of the i'th and j'th local problem solutions
          assert(allLocalSolutionEvaluations.size() == rows /*numMacroBaseFunctions*/);
          const auto gradLocProbSoli = allLocalSolutionEvaluations[ii][localQuadraturePoint];
          const auto gradLocProbSolj = allLocalSolutionEvaluations[jj][localQuadraturePoint];

          auto reconstructionGradPhii = coarseBaseJacs_[ii];
          reconstructionGradPhii += gradLocProbSoli;
          auto reconstructionGradPhij = coarseBaseJacs_[jj];
          reconstructionGradPhij += gradLocProbSolj;
          JacobianRangeType diffusive_flux(0.0);
          diffusion_operator->diffusiveFlux(global_point_in_U_T, reconstructionGradPhii, diffusive_flux);
          local_integral += diffusive_flux[0] * reconstructionGradPhij[0];
          retRow[jj] += local_integral * integrationFactor * quadratureWeight;
        }
      } // compute integral
    } // loop over all quadrature points
  } // ... apply(...)

private:
  const size_t over_integrate_;
}; // class Codim0Integral


// forward, to be used in the traits
template <class>
class MsFemCodim0Matrix;

class MsFemCodim0MatrixTraits
{
public:
  typedef MsFEMCodim0Integral LocalOperatorType;
  typedef MsFemCodim0Matrix<LocalOperatorType> derived_type;
}; // class LocalAssemblerCodim0MatrixTraits


template< class /*TYpe necessary for system assmebler usage only*/>
class MsFemCodim0Matrix
{
public:
  typedef MsFemCodim0MatrixTraits Traits;
  typedef typename Traits::LocalOperatorType LocalOperatorType;

  MsFemCodim0Matrix(const LocalOperatorType& op, LocalGridList& localGridList)
    : localOperator_(op)
    , localGridList_(localGridList)
  {}

  const LocalOperatorType& localOperator() const
  {
    return localOperator_;
  }

private:
  static const size_t numTmpObjectsRequired_ = 1;

public:
  std::vector< size_t > numTmpObjectsRequired() const
  {
    return { numTmpObjectsRequired_, localOperator_.numTmpObjectsRequired() };
  }

  /**
   *  \tparam T           Traits of the SpaceInterface implementation, representing the type of testSpace
   *  \tparam A           Traits of the SpaceInterface implementation, representing the type of ansatzSpace
   *  \tparam EntityType  A model of Dune::Entity< 0 >
   *  \tparam M           Traits of the Dune::Stuff::LA::Container::MatrixInterface implementation, representing the type of systemMatrix
   *  \tparam R           RangeFieldType, i.e. double
   */
  template< class EntityType, class M, class R >
  void assembleLocal(const CommonTraits::GdtSpaceType& testSpace,
                     const CommonTraits::GdtSpaceType& ansatzSpace,
                     const EntityType& coarse_grid_entity,
                     Dune::Stuff::LA::MatrixInterface< M >& systemMatrix,
                     std::vector< std::vector< Dune::DynamicMatrix< R > > >& tmpLocalMatricesContainer,
                     std::vector< Dune::DynamicVector< size_t > >& tmpIndicesContainer) const
  {
    // check
    assert(tmpLocalMatricesContainer.size() >= 1);
    assert(tmpLocalMatricesContainer[0].size() >= numTmpObjectsRequired_);
    assert(tmpLocalMatricesContainer[1].size() >= localOperator_.numTmpObjectsRequired());
    assert(tmpIndicesContainer.size() >= 2);
    // get and clear matrix

    Multiscale::MsFEM::LocalSolutionManager localSolutionManager(testSpace,
                                                                 coarse_grid_entity,
                                                                 localGridList_);
    localSolutionManager.load();
    const auto& localSolutions = localSolutionManager.getLocalSolutions();
    assert(localSolutions.size() > 0);

    for (const auto& localGridEntity : DSC::viewRange(*localSolutionManager.space().grid_view())) {
      // ignore overlay elements
      if (!localGridList_.covers(coarse_grid_entity, localGridEntity))
        continue;

      Dune::DynamicMatrix< R >& localMatrix = tmpLocalMatricesContainer[0][0];
      Dune::Stuff::Common::clear(localMatrix);
      auto& tmpOperatorMatrices = tmpLocalMatricesContainer[1];
      // apply local operator (result is in localMatrix)
      localOperator_.apply(localSolutionManager, localGridEntity,
                           testSpace.base_function_set(coarse_grid_entity),
                           ansatzSpace.base_function_set(coarse_grid_entity),
                           localMatrix,
                           tmpOperatorMatrices);
      // write local matrix to global
      Dune::DynamicVector< size_t >& globalRows = tmpIndicesContainer[0];
      Dune::DynamicVector< size_t >& globalCols = tmpIndicesContainer[1];
      const size_t rows = testSpace.mapper().numDofs(coarse_grid_entity);
      const size_t cols = ansatzSpace.mapper().numDofs(coarse_grid_entity);
      assert(globalRows.size() >= rows);
      assert(globalCols.size() >= cols);
      testSpace.mapper().globalIndices(coarse_grid_entity, globalRows);
      ansatzSpace.mapper().globalIndices(coarse_grid_entity, globalCols);
      for (size_t ii = 0; ii < rows; ++ii) {
        const auto& localRow = localMatrix[ii];
        const size_t globalII = globalRows[ii];
        for (size_t jj = 0; jj < cols; ++jj) {
          const size_t globalJJ = globalCols[jj];
          systemMatrix.add_to_entry(globalII, globalJJ, localRow[jj]);
        }
      } // write local matrix to global
    }
  } // ... assembleLocal(...)

private:
  const LocalOperatorType& localOperator_;
  LocalGridList& localGridList_;
}; // class LocalAssemblerCodim0Matrix

class EllipticCGMsFEM;

class EllipticCGMsFEMTraits
{
public:
  typedef EllipticCGMsFEM derived_type;
  typedef CommonTraits::LinearOperatorType MatrixType;
  typedef CommonTraits::GdtSpaceType SourceSpaceType;
  typedef CommonTraits::GdtSpaceType RangeSpaceType;
  typedef CommonTraits::GridViewType GridViewType;
}; // class EllipticCGTraits

class EllipticCGMsFEM
  : public GDT::Operators::MatrixBased< EllipticCGMsFEMTraits >
  , public GDT::SystemAssembler< EllipticCGMsFEMTraits::RangeSpaceType, EllipticCGMsFEMTraits::GridViewType,
                                 EllipticCGMsFEMTraits::SourceSpaceType, MsFemCodim0Matrix>
{

  typedef GDT::Operators::MatrixBased< EllipticCGMsFEMTraits > OperatorBaseType;
  typedef MsFEMCodim0Integral LocalOperatorType;
  typedef MsFemCodim0Matrix< LocalOperatorType > LocalAssemblerType;
  typedef GDT::SystemAssembler< EllipticCGMsFEMTraits::RangeSpaceType, EllipticCGMsFEMTraits::GridViewType,
  EllipticCGMsFEMTraits::SourceSpaceType, MsFemCodim0Matrix> AssemblerBaseType;
public:
  typedef EllipticCGMsFEMTraits Traits;

  typedef typename Traits::MatrixType       MatrixType;
  typedef typename Traits::SourceSpaceType  SourceSpaceType;
  typedef typename Traits::RangeSpaceType   RangeSpaceType;
  typedef typename Traits::GridViewType     GridViewType;

  using OperatorBaseType::pattern;

  static Stuff::LA::SparsityPatternDefault pattern(const RangeSpaceType& range_space,
                                                   const SourceSpaceType& source_space,
                                                   const GridViewType& grid_view)
  {
    return range_space.compute_volume_pattern(grid_view, source_space);
  }


  EllipticCGMsFEM(MatrixType& mtrx,
             const SourceSpaceType& src_spc,
             LocalGridList& localGridList)
    : OperatorBaseType(mtrx, src_spc)
    , AssemblerBaseType(src_spc)
    , local_assembler_(local_operator_, localGridList)
  {
    this->add(local_assembler_, this->matrix());
  }

  virtual ~EllipticCGMsFEM() {}

  virtual void assemble() DS_OVERRIDE DS_FINAL
  {
    AssemblerBaseType::assemble();
  }

private:
  const LocalOperatorType local_operator_;
  const LocalAssemblerType local_assembler_;
}; // class EllipticCG

CoarseScaleOperator::CoarseScaleOperator(const CoarseDiscreteFunctionSpace& coarseDiscreteFunctionSpace,
                                         LocalGridList& subgrid_list, const CommonTraits::DiffusionType& diffusion_op)
  : coarseDiscreteFunctionSpace_(coarseDiscreteFunctionSpace)
  , subgrid_list_(subgrid_list)
  , diffusion_operator_(diffusion_op)
  , petrovGalerkin_(false)
  , global_matrix_(coarseDiscreteFunctionSpace.mapper().size(), coarseDiscreteFunctionSpace.mapper().size(),
                   CommonTraits::EllipticOperatorType::pattern(coarseDiscreteFunctionSpace)) 
{
    DSC_PROFILER.startTiming("msfem.assembleMatrix");

  const bool is_simplex_grid = DSG::is_simplex_grid(coarseDiscreteFunctionSpace_);
  assert(!is_simplex_grid); //special casing isn't ported yet.

  GDT::SystemAssembler<CommonTraits::DiscreteFunctionSpaceType> global_system_assembler_(coarseDiscreteFunctionSpace_);
  EllipticCGMsFEM elliptic_operator(global_matrix_, coarseDiscreteFunctionSpace_, subgrid_list);
  global_system_assembler_.add(elliptic_operator);
  global_system_assembler_.assemble();

#if 0 //alter referenzcode
  for (const auto& coarse_grid_entity : DSC::viewRange(*coarseDiscreteFunctionSpace_.grid_view())) {
    int cacheCounter = 0;
    const auto& coarse_grid_geometry = coarse_grid_entity.geometry();
    assert(coarse_grid_entity.partitionType() == InteriorEntity);

    //      DSFe::LocalMatrixProxy<MatrixType> local_matrix(global_matrix_, coarse_grid_entity, coarse_grid_entity);
    auto local_matrix = nullptr;

    const auto& coarse_grid_baseSet = local_matrix.domainBasisFunctionSet();
    const auto numMacroBaseFunctions = coarse_grid_baseSet.size();

    Multiscale::MsFEM::LocalSolutionManager localSolutionManager(coarseDiscreteFunctionSpace_,
                                                                 coarse_grid_entity,
                                                                 subgrid_list_);
    localSolutionManager.load();
    const auto& localSolutions = localSolutionManager.getLocalSolutions();
    assert(localSolutions.size() > 0);

    for (const auto& localGridEntity : DSC::viewRange(*localSolutionManager.space().grid_view())) {
      // ignore overlay elements
      if (!subgrid_list_.covers(coarse_grid_entity, localGridEntity))
        continue;
      const auto& local_grid_geometry = localGridEntity.geometry();

      // higher order quadrature, since A^{\epsilon} is highly variable
      const auto localQuadrature = DSFe::make_quadrature(localGridEntity, localSolutionManager.space());
      const auto numQuadraturePoints = localQuadrature.nop();

      // number of local solutions without the boundary correctors. Those are only needed for the right hand side
      const auto numLocalSolutions = localSolutions.size() - localSolutionManager.numBoundaryCorrectors();
      // evaluate the jacobians of all local solutions in all quadrature points
      std::vector<std::vector<JacobianRangeType>> allLocalSolutionEvaluations(
            numLocalSolutions, std::vector<JacobianRangeType>(localQuadrature.nop(), JacobianRangeType(0.0)));
      for (auto lsNum : DSC::valueRange(numLocalSolutions)) {
        auto& sll = localSolutions[lsNum];
        assert(sll.get());
        //! \todo re-enable
        //            assert(localSolutionManager.space().indexSet().contains(localGridEntity));
        auto localFunction = sll->local_function(localGridEntity);
        localFunction.evaluateQuadrature(localQuadrature, allLocalSolutionEvaluations[lsNum]);
      }

      for (size_t localQuadraturePoint = 0; localQuadraturePoint < numQuadraturePoints; ++localQuadraturePoint) {
        // local (barycentric) coordinates (with respect to entity)
        const auto& local_subgrid_point = localQuadrature.point(localQuadraturePoint);

        auto global_point_in_U_T = local_grid_geometry.global(local_subgrid_point);
        const double weight_local_quadrature = localQuadrature.weight(localQuadraturePoint) *
                                               local_grid_geometry.integrationElement(local_subgrid_point);

        // evaluate the jacobian of the coarse grid base set
        const auto& local_coarse_point = coarse_grid_geometry.local(global_point_in_U_T);
        coarse_grid_baseSet.jacobianAll(local_coarse_point, gradientPhi);

        for (unsigned int i = 0; i < numMacroBaseFunctions; ++i) {
          for (unsigned int j = 0; j < numMacroBaseFunctions; ++j) {
            CoarseDiscreteFunctionSpace::RangeType local_integral(0.0);

            // Compute the gradients of the i'th and j'th local problem solutions
            JacobianRangeType gradLocProbSoli(0.0), gradLocProbSolj(0.0);
            if (is_simplex_grid) {
              assert(allLocalSolutionEvaluations.size() == CommonTraits::GridType::dimension);
              // ∇ Phi_H + ∇ Q( Phi_H ) = ∇ Phi_H + ∂_x1 Phi_H ∇Q( e_1 ) + ∂_x2 Phi_H ∇Q( e_2 )
              for (int k = 0; k < CommonTraits::GridType::dimension; ++k) {
                gradLocProbSoli.axpy(gradientPhi[i][0][k], allLocalSolutionEvaluations[k][localQuadraturePoint]);
                gradLocProbSolj.axpy(gradientPhi[j][0][k], allLocalSolutionEvaluations[k][localQuadraturePoint]);
              }
            } else {
              assert(allLocalSolutionEvaluations.size() == numMacroBaseFunctions);
              gradLocProbSoli = allLocalSolutionEvaluations[i][localQuadraturePoint];
              gradLocProbSolj = allLocalSolutionEvaluations[j][localQuadraturePoint];
            }

            JacobianRangeType reconstructionGradPhii(coarseBaseJacs_[cacheCounter][i]);
            reconstructionGradPhii += gradLocProbSoli;
            JacobianRangeType reconstructionGradPhij(coarseBaseJacs_[cacheCounter][j]);
            reconstructionGradPhij += gradLocProbSolj;
            JacobianRangeType diffusive_flux(0.0);
            diffusion_operator_.diffusiveFlux(global_point_in_U_T, reconstructionGradPhii, diffusive_flux);
            if (petrovGalerkin_)
              local_integral += weight_local_quadrature * (diffusive_flux[0] * coarseBaseJacs_[cacheCounter][j][0]);
            else
              local_integral += weight_local_quadrature * (diffusive_flux[0] * reconstructionGradPhij[0]);

            // add entries
            local_matrix.add(j, i, local_integral);
          }
        }
      }
    }
  } // for


  // set unit rows for dirichlet dofs
  Dune::Multiscale::getConstraintsCoarse(coarseDiscreteFunctionSpace_).applyToOperator(global_matrix_);
#endif
  DSC_PROFILER.stopTiming("msfem.assembleMatrix");
  DSC_LOG_DEBUG << "Time to assemble and communicate MsFEM matrix: " << DSC_PROFILER.getTiming("msfem.assembleMatrix") << "ms\n";
}

void CoarseScaleOperator::apply_inverse(const CoarseScaleOperator::CoarseDiscreteFunction& rhs,
                                        CoarseScaleOperator::CoarseDiscreteFunction& solution) {
  BOOST_ASSERT_MSG(rhs.dofs_valid(), "Coarse scale RHS DOFs need to be valid!");
  DSC_PROFILER.startTiming("msfem.solveCoarse");
  const typename BackendChooser<CoarseDiscreteFunctionSpace>::InverseOperatorType inverse(
      global_matrix_, rhs.space().communicator());/*, 1e-8, 1e-8, DSC_CONFIG_GET("msfem.solver.iterations", rhs.size()),
      DSC_CONFIG_GET("msfem.solver.verbose", false), "bcgs",
      DSC_CONFIG_GET("msfem.solver.preconditioner_type", std::string("sor")));*/
  inverse.apply(rhs.vector(), solution.vector());
  if (!solution.dofs_valid())
    DUNE_THROW(InvalidStateException, "Degrees of freedom of coarse solution are not valid!");
  DSC_PROFILER.stopTiming("msfem.solveCoarse");
  DSC_LOG_DEBUG << "Time to solve coarse MsFEM problem: " << DSC_PROFILER.getTiming("msfem.solveCoarse") << "ms."
               << std::endl;
} // constructor

} // namespace MsFEM {
} // namespace Multiscale {
} // namespace Dune {
