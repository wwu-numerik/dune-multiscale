// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_RIGHT_HAND_SIDE_ASSEMBLER_HH
#define DUNE_RIGHT_HAND_SIDE_ASSEMBLER_HH

#include <config.h>
#include <dune/fem/quadrature/quadrature.hh>
#include <dune/multiscale/hmm/cell_problem_numbering.hh>
#include <dune/multiscale/hmm/cell_problem_solver.hh>

#include <dune/multiscale/tools/misc.hh>
#include <dune/multiscale/common/traits.hh>
#include <dune/multiscale/msfem/localproblems/subgrid-list.hh>
#include <dune/multiscale/problems/selector.hh>
#include <dune/stuff/fem/functions/checks.hh>
#include <dune/multiscale/msfem/msfem_traits.hh>
#include <dune/multiscale/msfem/localproblems/localsolutionmanager.hh>
#include <dune/multiscale/common/dirichletconstraints.hh>
#include <dune/multiscale/common/traits.hh>

#include <dune/stuff/common/logging.hh>
#include <dune/stuff/functions/interfaces.hh>

namespace Dune {
namespace Multiscale {

//! Assembler for right rand side
//! We assemble the right hand side in a LSE, i.e. f \cdot \Phi_H + G \cdot \nabala \Phi_H
//! we call f the first Source and G the second Source
template <class DiscreteFunctionImp>
class RightHandSideAssembler {
private:
  typedef DiscreteFunctionImp DiscreteFunctionType;
  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

  typedef typename DiscreteFunctionSpaceType::BasisFunctionSetType BasisFunctionSetType;
  typedef typename DiscreteFunctionSpaceType::RangeType RangeType;
  typedef typename DiscreteFunctionSpaceType::DomainFieldType DomainFieldType;
  typedef DomainFieldType TimeType;
  typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
  typedef typename GridPartType::GridType GridType;
  typedef typename DiscreteFunctionSpaceType::JacobianRangeType JacobianRangeType;
  typedef typename DiscreteFunctionSpaceType::DomainType DomainType;
  typedef typename GridType::template Codim<0>::Entity EntityType;
  typedef typename EntityType::Geometry GeometryType;

  typedef MsFEM::LocalSolutionManager LocalSolutionManagerType;

  static const int dimension = GridType::dimension;
  static const int polynomialOrder = DiscreteFunctionSpaceType::polynomialOrder;

  // need a virtual base to work around local classes not being allowed in templated scopes
  struct FunctorBase {
    virtual RangeType operator()(const DomainType& global_quad_point, const JacobianRangeType& gradientPhi) const = 0;
    virtual ~FunctorBase() {}
  };

  template <class FirstSourceType>
  static void assemble_common(const FirstSourceType& f, const FunctorBase& functor, const int polOrd,
                              DiscreteFunctionType& rhsVector) {
    // set rhsVector to zero:
    rhsVector.clear();
    for (const auto& entity : rhsVector.space()) {
      const auto& geometry = entity.geometry();            // Referenz auf Geometrie
      auto elementOfRHS = rhsVector.localFunction(entity); // entity zeigt auf ein bestimmtes Element der
                                                           // entity
      // hier wird sozusagen ein Pointer von localFunction auf discreteFunction erzeugt. Befinden wir uns auf einer
      // bestimmten entity, so berechnet localFunction alle noetigen Werte und speichert sie (da Pointer) in
      // discreteFunction(aktuelleEntity)

      const BasisFunctionSetType baseSet // BaseFunctions leben immer auf Refernzelement!!!
          = rhsVector.space().basisFunctionSet(
              entity); // entity Referenz auf eine bestimmtes Element der entity. In der
                       // ersten Klasse war das Element fest, deshalb konnte man sich
                       // dort Pointer sparen. //loeschen: discreteFunctionSpace
                       // statt
                       // functionSpace

      const auto quadrature = make_quadrature(entity, rhsVector.space(), polOrd);
      const auto numDofs = elementOfRHS.numDofs();
      for (auto quadraturePoint : DSC::valueRange(quadrature.nop())) {
        // the return values:
        RangeType f_x;
        std::vector<RangeType> phi_x(numDofs);
        std::vector<JacobianRangeType> gradientPhi(numDofs);
        // to save: A \nabla PHI_H * \nabla phi_h;
        RangeType res = 0;
        const auto det = geometry.integrationElement(quadrature.point(quadraturePoint));

        f.evaluate(geometry.global(quadrature.point(quadraturePoint)), f_x);
        baseSet.evaluateAll(quadrature[quadraturePoint], phi_x);
        baseSet.jacobianAll(quadrature[quadraturePoint], gradientPhi);

        for (int i = 0; i < numDofs; ++i) {
          res = functor(geometry.global(quadrature.point(quadraturePoint)), gradientPhi[i]);
          elementOfRHS[i] += det * quadrature.weight(quadraturePoint) * (f_x * phi_x[i]);
          elementOfRHS[i] += det * quadrature.weight(quadraturePoint) * (res);
        }
      }
    }
  } // end method

public:
  /** assemble standard right hand side:
   * if there is only one source (f) (there is no second source):
   * discreteFunction is an output parameter (kind of return value)
   **/
  template <int polOrd, class FirstSourceType>
  static void assemble(const FirstSourceType& f, DiscreteFunctionType& rhsVector) {
    struct Functor : public FunctorBase {
      RangeType operator()(const DomainType&, const JacobianRangeType&) const { return RangeType(0.0); }
    } functor;
    assemble_common(f, functor, polOrd, rhsVector);
  } // end method

  /** if there is a first source f and a second source G:
   * discreteFunction is an output parameter (kind of return value)
   **/
  template <int polOrd, class FirstSourceType, class SecondSourceType>
  static void assemble(const FirstSourceType& f, const SecondSourceType& _G, DiscreteFunctionType& rhsVector) {
    struct Functor : public FunctorBase {
      const SecondSourceType& G;
      Functor(const SecondSourceType& __G) : G(__G) {}

      RangeType operator()(const DomainType& global_quad_point, const JacobianRangeType& gradientPhi) const {
        RangeType res(0.0);
        RangeType G_x(0.0);
        // evaluate the second source at the current quadrature point and save its value in 'G_x':
        for (int k = 0; k < dimension; ++k) {
          G_x = 0;
          G.evaluate(k, global_quad_point, G_x);
          res += G_x * gradientPhi[0][k];
        }
        return res;
      }
    } functor(_G);
    assemble_common(f, functor, polOrd, rhsVector);
  } // end method

  /** if there is a first source f, a second source G and a parameter t:
   * discreteFunction is an output parameter (kind of return value)
   **/
  template <int polOrd, class FirstSourceType, class SecondSourceType>
  static void assemble(const FirstSourceType& f, const SecondSourceType& G, const TimeType& t,
                       DiscreteFunctionType& rhsVector) {
    struct Functor : public FunctorBase {
      const SecondSourceType& G;
      const TimeType& t;
      Functor(const SecondSourceType& __G, const TimeType& _t)
        : G(__G)
        , t(_t) {}

      RangeType operator()(const DomainType& global_quad_point, const JacobianRangeType& gradientPhi) const {
        RangeType res(0.0);
        RangeType G_x(0.0);
        // evaluate the second source at the current quadrature point and save its value in 'G_x':
        for (int k = 0; k < dimension; ++k) {
          G_x = 0;
          G.evaluate(k, global_quad_point, t, G_x);
          res += G_x * gradientPhi[0][k];
        }
        return res;
      }
    } functor(G, t);
    assemble_common(f, functor, polOrd, rhsVector);
  } // end method

  /**
   * The rhs-assemble()-methods for linear elliptic problems
   * with non-homogeneous Dirichlet and Neumann boundary conditions:
   **/

  template <int polOrd, class FirstSourceType, class DiffusionOperatorType, class NeumannBCType>
  static void
  assemble(const FirstSourceType& f, const DiffusionOperatorType& A,
           const DiscreteFunctionType& dirichlet_extension, // discrete function describing dirichlet extension
           const NeumannBCType& neumann_bc, DiscreteFunctionType& rhsVector) {
    rhsVector.clear();

    for (const auto& entity : rhsVector.space()) {

      const auto& geometry = entity.geometry();
      auto elementOfRHS = rhsVector.localFunction(entity);
      const auto baseSet = rhsVector.space().basisFunctionSet(entity);

      const int numDofs = elementOfRHS.numDofs();

      std::vector<RangeType> phi_x(numDofs);
      // gradient of base function and gradient of old_u_H
      std::vector<JacobianRangeType> grad_phi_x(numDofs);

      const auto loc_dirichlet_extension = dirichlet_extension.localFunction(entity);
      const auto quadrature = make_quadrature(entity, rhsVector.space(), polOrd);

      const auto& lagrangePointSet = rhsVector.space().lagrangePointSet(entity);

      for (const auto& intersection : DSC::intersectionRange(rhsVector.space().gridPart(), entity)) {
        if (Problem::isNeumannBoundary(intersection)) {
          const auto face = intersection.indexInInside();
          const auto faceQuadrature = make_quadrature(intersection, rhsVector.space(), polOrd);
          const auto numFaceQuadraturePoints = faceQuadrature.nop();

          for (auto faceQuadraturePoint : DSC::valueRange(numFaceQuadraturePoints)) {
            baseSet.evaluateAll(faceQuadrature[faceQuadraturePoint], phi_x);
            baseSet.jacobianAll(faceQuadrature[faceQuadraturePoint], grad_phi_x);

            const auto local_point_entity = faceQuadrature.point(faceQuadraturePoint);
            const auto global_point = geometry.global(local_point_entity);
            const auto local_point_face = intersection.geometry().local(global_point);

            RangeType neumann_value(0.0);
            neumann_bc.evaluate(global_point, neumann_value);

            const double face_weight = intersection.geometry().integrationElement(local_point_face) *
                                       faceQuadrature.weight(faceQuadraturePoint);

            for (const auto& lp : DSC::lagrangePointSetRange(lagrangePointSet, face)) {
              elementOfRHS[lp] += neumann_value * face_weight * phi_x[lp];
            }
          }
        }
      }

      // the return values:
      RangeType f_x;

      JacobianRangeType gradient_dirichlet_extension;
      JacobianRangeType diffusive_flux_in_gradient_dirichlet_extension;

      for (auto quadraturePoint : DSC::valueRange(quadrature.nop())) {
        // local (barycentric) coordinates (with respect to entity)
        const auto& local_point = quadrature.point(quadraturePoint);
        const auto global_point = geometry.global(local_point);

        const double weight = geometry.integrationElement(local_point) * quadrature.weight(quadraturePoint);
        // evaluate the Right Hand Side Function f at the current quadrature point and save its value in 'f_y':
        f.evaluate(global_point, f_x);
        // evaluate the current base function at the current quadrature point and save its value in 'z':
        baseSet.evaluateAll(quadrature[quadraturePoint], phi_x); // i = i'te Basisfunktion;
        // evaluate the gradient of the current base function at the current quadrature point and save its value in
        // 'returnGradient':
        baseSet.jacobianAll(quadrature[quadraturePoint], grad_phi_x);
        // get gradient of dirichlet extension:
        loc_dirichlet_extension.jacobian(quadrature[quadraturePoint], gradient_dirichlet_extension);
        A.diffusiveFlux(global_point, gradient_dirichlet_extension, diffusive_flux_in_gradient_dirichlet_extension);

        for (int i = 0; i < numDofs; ++i) {
          elementOfRHS[i] += weight * (f_x * phi_x[i]);
          elementOfRHS[i] -= weight * (diffusive_flux_in_gradient_dirichlet_extension[0] * grad_phi_x[i][0]);
        }
      }
    }
  } // end method

  /** assemble right hand side (if there is only one source - f):
   *  assemble-method for MsFEM in symmetric (non-Petrov-Galerkin) formulation
   *  rhsVector is the output parameter (kind of return value)
   **/
  template <int polOrd, class FirstSourceType, class MacroMicroGridSpecifierType, class SubGridListType>
  static void assemble_for_MsFEM_symmetric(const FirstSourceType& f, MacroMicroGridSpecifierType& specifier,
                                           SubGridListType& subgrid_list, DiscreteFunctionType& rhsVector) {

    // gather some problem data
    auto diffusionPtr = Problem::getDiffusion();
    const auto& diffusion = *diffusionPtr;
    auto neumannDataPtr = Problem::getNeumannData();
    const auto& neumannData = *neumannDataPtr;

    DiscreteFunctionType dirichletExtension("Dirichlet Extension", specifier.fineSpace());
    dirichletExtension.clear();
    Dune::Multiscale::projectDirichletValues(rhsVector.space(), dirichletExtension);

    // set rhsVector to zero:
    rhsVector.clear();
    const auto& coarseGridLeafIndexSet = specifier.coarseSpace().gridPart().grid().leafIndexSet();
    RangeType f_x;
    for (const auto& coarse_grid_entity : rhsVector.space()) {
      const auto coarseEntityIndex = coarseGridLeafIndexSet.index(coarse_grid_entity);

      const auto& coarseGeometry = coarse_grid_entity.geometry();
      auto rhsLocalFunction = rhsVector.localFunction(coarse_grid_entity);
      const auto numLocalBaseFunctions = rhsLocalFunction.numDofs();

      const BasisFunctionSetType& coarse_grid_baseSet = specifier.coarseSpace().basisFunctionSet(coarse_grid_entity);

      // --------- add standard contribution of right hand side -------------------------
      {
        //!\TODO warum +5 ???
        const auto quadrature = make_quadrature(coarse_grid_entity, rhsVector.space(), polOrd + 5);
        std::vector<RangeType> phi_x_vec(numLocalBaseFunctions);
        const auto numQuadraturePoints = quadrature.nop();
        for (size_t quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint) {
          const double det = coarseGeometry.integrationElement(quadrature.point(quadraturePoint));
          // evaluate the Right Hand Side Function f at the current quadrature point and save its value in 'f_y':
          f.evaluate(coarseGeometry.global(quadrature.point(quadraturePoint)), f_x);
          coarse_grid_baseSet.evaluateAll(quadrature[quadraturePoint], phi_x_vec);
          for (int i = 0; i < numLocalBaseFunctions; ++i) {
            rhsLocalFunction[i] += det * quadrature.weight(quadraturePoint) * (f_x * phi_x_vec[i]);
          }
        }
      }
      // ----------------------------------------------------------------------------------

      // --------- add corrector contribution of right hand side --------------------------
      // Load local solutions
      LocalSolutionManagerType localSolutionManager(coarse_grid_entity, subgrid_list, specifier);
      localSolutionManager.loadLocalSolutions();
      LocalSolutionManagerType::LocalSolutionVectorType& localSolutions = localSolutionManager.getLocalSolutions();
      assert(localSolutions.size() > 0);

      // iterator for the micro grid ( grid for the reference element T_0 )
      const auto& subGrid = subgrid_list.getSubGrid(coarse_grid_entity);
      for (const auto& localEntity : DSC::viewRange(subGrid.leafView())) {
        const auto& hostCell = subGrid.template getHostEntity<0>(localEntity);
        const auto enclosingCoarseCellIndex = subgrid_list.getEnclosingMacroCellIndex(hostCell);
        auto dirichletExtensionLF = dirichletExtension.localFunction(*hostCell);
        if (enclosingCoarseCellIndex == coarseEntityIndex) {
          // higher order quadrature, since A^{\epsilon} is highly variable
          const auto localQuadrature =
              make_quadrature(localEntity, localSolutionManager.getLocalDiscreteFunctionSpace());

          // evaluate all local solutions and their jacobians in all quadrature points
          std::vector<std::vector<RangeType>> allLocalSolutionEvaluations(
              localSolutions.size(), std::vector<RangeType>(localQuadrature.nop(), 0.0));
          std::vector<std::vector<JacobianRangeType>> allLocalSolutionJacobians(
              localSolutions.size(), std::vector<JacobianRangeType>(localQuadrature.nop(), JacobianRangeType(0.0)));
          for (auto lsNum : DSC::valueRange(localSolutions.size())) {
            auto localFunction = localSolutions[lsNum]->localFunction(localEntity);
            // this evaluates the local solutions in all quadrature points...
            localFunction.evaluateQuadrature(localQuadrature, allLocalSolutionEvaluations[lsNum]);
            // while this automatically evaluates their jacobians.
            localFunction.evaluateQuadrature(localQuadrature, allLocalSolutionJacobians[lsNum]);

            const auto& subGridPart = localSolutionManager.getSubGridPart();
            for (const auto& intersection : DSC::intersectionRange(subGridPart.grid().leafView(), localEntity)) {
              if (Problem::isNeumannBoundary(intersection)) {
                const int orderOfIntegrand = (polynomialOrder - 1) + 2 * (polynomialOrder + 1);
                const int quadOrder = std::ceil((orderOfIntegrand + 1) / 2);
                // get type of face quadrature. Is done in this scope because Patricks methods use another type.
                const auto faceQuad = make_quadrature(intersection, localSolutions[lsNum]->space(), quadOrder);
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

                  // the following does not work because subgrid does not implement geometryInInside()
                  // const auto& insideGeometry    = intersection.geometryInInside();
                  // const typename FaceQuadratureType::CoordinateType& xInInside = insideGeometry.global(xLocal);
                  // therefore, we have to do stupid things:
                  const auto& xGlobal = faceGeometry.global(xLocal);
                  const auto& xInCoarseLocal = coarse_grid_entity.geometry().local(xGlobal);
                  const double factor = faceGeometry.integrationElement(xLocal) * faceQuad.weight(iqP);

                  neumannData.evaluate(xGlobal, neumannValue);
                  coarse_grid_baseSet.evaluateAll(xInCoarseLocal, phi_x_vec);
                  for (auto i : DSC::valueRange(numLocalBaseFunctions)) {
                    assert((long long)i < (long long)phi_x_vec.size());
                    assert(iqP < localSolutionOnFace.size());
                    rhsLocalFunction[i] += factor * (neumannValue * (phi_x_vec[i] + localSolutionOnFace[iqP]));
                  }
                }
              }
            }
          }

          const auto& localGeometry = localEntity.geometry();
          RangeType corrector_phi_x;
          for (size_t qP = 0; qP < localQuadrature.nop(); ++qP) {
            // local (barycentric) coordinates (with respect to entity)
            const auto& quadPoint = localQuadrature.point(qP);
            const auto quadPointGlobal = localGeometry.global(quadPoint);

            const double quadWeight = localQuadrature.weight(qP) * localGeometry.integrationElement(quadPoint);

            for (int coarseBF = 0; coarseBF < numLocalBaseFunctions; ++coarseBF) {
              JacobianRangeType diffusive_flux(0.0);

              // evaluate gradient of basis function
              const auto quadPointLocalInCoarse = coarseGeometry.local(quadPointGlobal);
              std::vector<JacobianRangeType> gradient_Phi_vec(numLocalBaseFunctions);
              coarse_grid_baseSet.jacobianAll(quadPointLocalInCoarse, gradient_Phi_vec);

              JacobianRangeType reconstructionGradPhi(gradient_Phi_vec[coarseBF]);

              if (specifier.simplexCoarseGrid()) {
                assert(localSolutions.size() == GridSelector::dimgrid + localSolutionManager.numBoundaryCorrectors());
                DUNE_THROW(NotImplemented, "Boundary values are not implemented for simplex grids yet!");
              } else {
                assert(localSolutions.size() == numLocalBaseFunctions + localSolutionManager.numBoundaryCorrectors());
                // local corrector for coarse base func
                corrector_phi_x = allLocalSolutionEvaluations[coarseBF][qP];
                // element part of boundary conditions
                JacobianRangeType directionOfFlux(0.0);
                //! @attention At this point we assume, that the quadrature points on the subgrid and hostgrid
                //! are the same (dirichletExtensionLF is a localfunction on the hostgrid, quadPoint stems from
                //! a quadrature on the subgrid)!!
                dirichletExtensionLF.jacobian(quadPoint, directionOfFlux);
                // add dirichlet-corrector
                directionOfFlux += allLocalSolutionJacobians[numLocalBaseFunctions + 1][qP];
                // subtract neumann-corrector
                // directionOfFlux -= allLocalSolutionJacobians[numLocalBaseFunctions][qP];

                diffusion.diffusiveFlux(quadPointGlobal, directionOfFlux, diffusive_flux);
                reconstructionGradPhi += allLocalSolutionJacobians[coarseBF][qP];
              }
              f.evaluate(quadPointGlobal, f_x);
              double val = quadWeight * (f_x * corrector_phi_x);
              rhsLocalFunction[coarseBF] += val;
              val = quadWeight * (diffusive_flux[0] * reconstructionGradPhi[0]);
              rhsLocalFunction[coarseBF] -= val;
            }
          }
        }
      }
    }

    // set dirichlet dofs to zero
    Dune::Multiscale::getConstraintsCoarse(rhsVector.space()).setValue(0.0, rhsVector);
  } // end method



};  // end class
} // end namespace Multiscale
} // end namespace Dune

#endif // ifndef DUNE_RIGHT_HAND_SIDE_ASSEMBLER_HH
