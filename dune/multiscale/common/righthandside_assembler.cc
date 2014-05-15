#include <config.h>
#include <assert.h>
#include <dune/common/exceptions.hh>
#include <dune/multiscale/problems/selector.hh>
#include <dune/stuff/common/ranges.hh>
#include <dune/stuff/functions/femadapter.hh>
#include <dune/stuff/common/profiler.hh>
#include <dune/fem/misc/threads/domainthreaditerator.hh>
#include <memory>

#include <dune/stuff/fem/functions/integrals.hh>
#include <dune/stuff/functions/femadapter.hh>
#include <dune/multiscale/common/traits.hh>
#include <dune/multiscale/tools/misc.hh>

#include <dune/multiscale/msfem/localproblems/subgrid-list.hh>
#include <dune/multiscale/msfem/localproblems/localsolutionmanager.hh>

#include "righthandside_assembler.hh"

void Dune::Multiscale::MsFEM::RightHandSideAssembler::assemble(
    const CommonTraits::DiscreteFunctionSpaceType& coarse_space,
    const Dune::Multiscale::CommonTraits::SourceType& f, DMM::LocalGridList& subgrid_list,
    Dune::Multiscale::CommonTraits::DiscreteFunctionType& rhsVector) {
  DSC_PROFILER.startTiming("msfem.assembleRHS");
  cached_ = false;
  
  // cache grid variable
  const bool isSimplexGrid = DSG::is_simplex_grid(coarse_space);
  
  static constexpr int dimension = CommonTraits::GridType::dimension;
  const auto& diffusion = *Problem::getDiffusion();
  const auto& neumannData = *Problem::getNeumannData();

  rhsVector.vector() *= 0;
  RangeType f_x;
//#ifdef _OPENMP
//#pragma omp parallel
//#endif
  {
    for (const auto& coarse_grid_entity : DSC::viewRange(*rhsVector.space().grid_view())) {
      int cacheCounter = 0;
      const auto& coarseGeometry = coarse_grid_entity.geometry();
      auto rhsLocalFunction = rhsVector.local_discrete_function(coarse_grid_entity);
      const auto numLocalBaseFunctions = rhsLocalFunction->vector().size();
      const auto& coarseBaseSet = coarse_space.baseFunctionSet(coarse_grid_entity);

      // --------- add corrector contribution of right hand side --------------------------
      // Load local solutions
      MsFEM::LocalSolutionManager localSolutionManager(coarse_space, coarse_grid_entity, subgrid_list);
      localSolutionManager.load();
      auto& localSolutions = localSolutionManager.getLocalSolutions();
      assert(localSolutions.size() > 0);
      MsFEM::MsFEMTraits::LocalGridDiscreteFunctionType dirichletExtension("Dirichlet Extension",
                                                                           localSolutionManager.space());
      dirichletExtension.clear();
      Dune::Multiscale::copyDirichletValues(rhsVector.space(), dirichletExtension);

      // iterator for the micro grid ( grid for the reference element T_0 )
      const auto& subGrid = subgrid_list.getSubGrid(coarse_grid_entity);
      auto view = subGrid.leafGridView();
      for (const auto& localEntity : DSC::viewRange(view)) {
        if (subgrid_list.covers(coarse_grid_entity, localEntity)) {
          // higher order quadrature, since A^{\epsilon} is highly variable
          const auto localQuadrature = DSFe::make_quadrature(localEntity, localSolutionManager.space());
          auto dirichletExtensionLF = dirichletExtension.local_function(localEntity);

          // evaluate all local solutions and their jacobians in all quadrature points
          std::vector<std::vector<RangeType>> allLocalSolutionEvaluations(
              localSolutions.size(), std::vector<RangeType>(localQuadrature.nop(), 0.0));
          std::vector<std::vector<JacobianRangeType>> allLocalSolutionJacobians(
              localSolutions.size(), std::vector<JacobianRangeType>(localQuadrature.nop(), JacobianRangeType(0.0)));
          // add contributions for all inner correctors
          for (auto lsNum : DSC::valueRange(numLocalBaseFunctions)) {
            auto localFunction = localSolutions[lsNum]->local_function(localEntity);
            // this evaluates the local solutions in all quadrature points...
            localFunction.evaluateQuadrature(localQuadrature, allLocalSolutionEvaluations[lsNum]);
            // while this automatically evaluates their jacobians.
            localFunction.evaluateQuadrature(localQuadrature, allLocalSolutionJacobians[lsNum]);

            // assemble intersection-part
            const auto& view = localSolutionManager.grid_view();
            for (const auto& intersection : DSC::intersectionRange(view, localEntity)) {
              if (DMP::is_neumann(intersection)) {
                const int orderOfIntegrand =
                    (CommonTraits::polynomial_order - 1) + 2 * (CommonTraits::polynomial_order + 1);
                const int quadOrder = std::ceil((orderOfIntegrand + 1) / 2);
                const auto faceQuad = DSFe::make_quadrature(intersection, localSolutions[lsNum]->space(), quadOrder);
                RangeType neumannValue(0.0);
                const auto numQuadPoints = faceQuad.nop();
                // loop over all quadrature points

                std::vector<RangeType> phi_x_vec(numLocalBaseFunctions);
                std::vector<RangeType> localSolutionOnFace(numLocalBaseFunctions);
                localFunction.evaluateQuadrature(faceQuad, localSolutionOnFace);

                for (unsigned int iqP = 0; iqP < numQuadPoints; ++iqP) {
                  // get local coordinate of quadrature point
                  const auto& xLocal = faceQuad.localPoint(iqP);
                  const auto& faceGeometry = intersection.geometry();
                  const auto& xGlobal = faceGeometry.global(xLocal);
                  const auto& xInCoarseLocal = coarse_grid_entity.geometry().local(xGlobal);
                  const double factor = faceGeometry.integrationElement(xLocal) * faceQuad.weight(iqP);

                  neumannData.evaluate(xGlobal, neumannValue);
                  coarseBaseSet.evaluateAll(xInCoarseLocal, phi_x_vec);
                  assert((long long)lsNum < (long long)phi_x_vec.size());
                  assert(iqP < localSolutionOnFace.size());
                  rhsLocalFunction[lsNum] -= factor * (neumannValue * (phi_x_vec[lsNum] + localSolutionOnFace[iqP]));
                }
              }
            }
          }

          // element part of boundary conditions
          std::vector<JacobianRangeType> dirichletVals(localQuadrature.nop(), 0.0);
          dirichletExtensionLF.evaluateQuadrature(localQuadrature, dirichletVals);

          // assemble element-part
          const auto& localGeometry = localEntity.geometry();
          for (size_t qP = 0; qP < localQuadrature.nop(); ++qP) {
            // local (barycentric) coordinates (with respect to entity)
            const auto& quadPoint = localQuadrature.point(qP);
            const auto quadPointGlobal = localGeometry.global(quadPoint);

            const double quadWeight = localQuadrature.weight(qP) * localGeometry.integrationElement(quadPoint);

            // evaluate gradient of basis function
            const auto quadPointLocalInCoarse = coarseGeometry.local(quadPointGlobal);


            if (!cached_) {
              coarseBaseEvals_.push_back(std::vector<RangeType>(numLocalBaseFunctions, 0.0));
              coarseBaseJacs_.push_back(std::vector<JacobianRangeType>(numLocalBaseFunctions, 0.0));
              coarseBaseSet.evaluateAll(quadPointLocalInCoarse, coarseBaseEvals_[cacheCounter]);
              coarseBaseSet.jacobianAll(quadPointLocalInCoarse, coarseBaseJacs_[cacheCounter]);
            }

            for (int coarseBF = 0; coarseBF < numLocalBaseFunctions; ++coarseBF) {
              JacobianRangeType diffusive_flux(0.0);

              JacobianRangeType reconstructionGradPhi(coarseBaseJacs_[cacheCounter][coarseBF]);
              RangeType reconstructionPhi(coarseBaseEvals_[cacheCounter][coarseBF]);

              if (isSimplexGrid) {
                assert(localSolutions.size() == dimension + localSolutionManager.numBoundaryCorrectors());
                for (const auto& i : DSC::valueRange(dimension))
                  reconstructionPhi += coarseBaseJacs_[cacheCounter][coarseBF][0][i] * allLocalSolutionEvaluations[i][qP];
                //! @todo add the dirichlet and neumann-correctors!
              } else {
                assert(localSolutions.size() == numLocalBaseFunctions + localSolutionManager.numBoundaryCorrectors());
                // local corrector for coarse base func
                reconstructionPhi += allLocalSolutionEvaluations[coarseBF][qP];
                // element part of boundary conditions
                JacobianRangeType directionOfFlux = dirichletVals[qP];
                // add dirichlet-corrector
                directionOfFlux += allLocalSolutionJacobians[numLocalBaseFunctions + 1][qP];
                // subtract neumann-corrector
                directionOfFlux -= allLocalSolutionJacobians[numLocalBaseFunctions][qP];

                diffusion.diffusiveFlux(quadPointGlobal, directionOfFlux, diffusive_flux);
                reconstructionGradPhi += allLocalSolutionJacobians[coarseBF][qP];
              }
              f.evaluate(quadPointGlobal, f_x);

              rhsLocalFunction[coarseBF] += quadWeight * (f_x * reconstructionPhi);
              rhsLocalFunction[coarseBF] -= quadWeight * (diffusive_flux[0] * reconstructionGradPhi[0]);
            }
            ++cacheCounter;
          }
        }
      }
      cached_ = true;
    } // for
  }   // omp region

  // set dirichlet dofs to zero
  Dune::Multiscale::getConstraintsCoarse(rhsVector.space()).setValue(0.0, rhsVector);
  rhsVector.communicate();
  DSC_PROFILER.stopTiming("msfem.assembleRHS");
  DSC_LOG_DEBUG << "Time to assemble and communicate MsFEM rhs: " << DSC_PROFILER.getTiming("msfem.assembleRHS") << "ms\n";
}
