// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_RIGHT_HAND_SIDE_ASSEMBLER_HH
#define DUNE_RIGHT_HAND_SIDE_ASSEMBLER_HH

#ifdef HAVE_CMAKE_CONFIG
  #include "cmake_config.h"
#else
  #include "config.h"
#endif // ifdef HAVE_CMAKE_CONFIG

// - Dune includes
#include <dune/fem/quadrature/quadrature.hh>
#include <dune/multiscale/hmm/cell_problem_numbering.hh>
#include <dune/multiscale/hmm/cell_problem_solver.hh>

#include <dune/multiscale/tools/misc.hh>
#include <dune/multiscale/common/traits.hh>
#include <dune/multiscale/msfem/localproblems/subgrid-list.hh>
#include <dune/stuff/fem/functions/checks.hh>
#include <dune/multiscale/msfem/msfem_traits.hh>
#include <dune/multiscale/msfem/localproblems/localsolutionmanager.hh>

#include <dune/stuff/common/logging.hh>
#include <dune/stuff/function/interface.hh>

namespace Dune {




//! Assembler for right rand side
//! We assemble the right hand side in a LSE, i.e. f \cdot \Phi_H + G \cdot \nabala \Phi_H
//! we call f the first Source and G the second Source
template< class DiscreteFunctionImp >
class RightHandSideAssembler
{
private:
  typedef DiscreteFunctionImp DiscreteFunctionType;
  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType
    DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionType::LocalFunctionType
    LocalFunctionType;
  typedef typename DiscreteFunctionSpaceType::BasisFunctionSetType
    BasisFunctionSetType;
  typedef typename DiscreteFunctionSpaceType::RangeType RangeType;
  typedef typename DiscreteFunctionSpaceType::DomainFieldType DomainFieldType;
  typedef DomainFieldType TimeType;
  typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
  typedef typename GridPartType::GridType GridType;
  typedef typename DiscreteFunctionSpaceType::JacobianRangeType
    JacobianRangeType;
  typedef typename DiscreteFunctionSpaceType::DomainType
    DomainType;
  typedef typename GridType::template Codim< 0 >::Entity EntityType;
  typedef typename EntityType::Geometry GeometryType;
  typedef Fem::CachingQuadrature< GridPartType, 0 > Quadrature;
  typedef Fem::CachingQuadrature< GridPartType, 1 > FaceQuadrature;
  
  enum { dimension = GridType::dimension };

public:
  static void printRHS(const DiscreteFunctionType& rhs) {
    for (auto entity : rhs.space())
    {
      const LocalFunctionType elementOfRHS = rhs.localFunction(entity);
      const int numDofs = elementOfRHS.numDofs();
      for (int i = 0; i < numDofs; ++i)
      {
        DSC_LOG_DEBUG << "Number of Dof: " << i << " ; " << rhs.name() << " : " << elementOfRHS[i] << std::endl;
      }
    }
  }  // end method

  //! --------------- The rhs-assemble()-methods for linear elliptic problems -----------------

private:
  // need a virtual base to work around local classes not being allowed in templated scopes
  struct FunctorBase {
    virtual RangeType operator()(const DomainType& global_quad_point, const JacobianRangeType& gradientPhi) const = 0;
    virtual ~FunctorBase(){}
  };

  template< class FirstSourceType >
  static void assemble_common(const FirstSourceType& f,
                       const FunctorBase& functor,
                       const int polOrd,
                       DiscreteFunctionType& rhsVector)
  {
    // set rhsVector to zero:
    rhsVector.clear();
    for (const auto& entity : rhsVector.space())
    {
      const GeometryType& geometry = entity.geometry(); // Referenz auf Geometrie
      LocalFunctionType elementOfRHS = rhsVector.localFunction(entity);   // entity zeigt auf ein bestimmtes Element der
                                                                       // entity
      // hier wird sozusagen ein Pointer von localFunction auf discreteFunction erzeugt. Befinden wir uns auf einer
      // bestimmten entity, so berechnet localFunction alle noetigen Werte und speichert sie (da Pointer) in
      // discreteFunction(aktuelleEntity)

      const BasisFunctionSetType baseSet // BaseFunctions leben immer auf Refernzelement!!!
        = rhsVector.space().basisFunctionSet(entity);     // entity Referenz auf eine bestimmtes Element der entity. In der
                                                          // ersten Klasse war das Element fest, deshalb konnte man sich
                                                          // dort Pointer sparen. //loeschen: discreteFunctionSpace
                                                          // statt
                                                          // functionSpace

      const Fem::CachingQuadrature< GridPartType, 0 > quadrature(entity, polOrd);   // 0 --> codim 0
      const int numQuadraturePoints = quadrature.nop();
      const int numDofs = elementOfRHS.numDofs();
      for (int quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint)
      {
        // the return values:
        RangeType f_x;
        std::vector<RangeType> phi_x(numDofs);
        std::vector<JacobianRangeType> gradientPhi(numDofs);
        // to save: A \nabla PHI_H * \nabla phi_h;
        RangeType res = 0;
        const double det = geometry.integrationElement(quadrature.point(quadraturePoint));

        f.evaluate(geometry.global( quadrature.point(quadraturePoint) ), f_x);
        baseSet.evaluateAll(quadrature[quadraturePoint], phi_x);
        baseSet.jacobianAll(quadrature[quadraturePoint], gradientPhi);

        for (int i = 0; i < numDofs; ++i)
        {
          res = functor(geometry.global(quadrature.point(quadraturePoint)), gradientPhi[i]);
          elementOfRHS[i] += det * quadrature.weight(quadraturePoint) * (f_x * phi_x[i]);
          elementOfRHS[i] += det * quadrature.weight(quadraturePoint) * (res);
        }
      }
    }
  }  // end method

public:
  /** assemble standard right hand side:
   * if there is only one source (f) (there is no second source):
   * discreteFunction is an output parameter (kind of return value)
   **/
  template< int polOrd, class FirstSourceType >
  static void assemble(const FirstSourceType& f,
                DiscreteFunctionType& rhsVector) {
    struct Functor : public FunctorBase {
      RangeType operator()(const DomainType&, const JacobianRangeType& ) const {
        return RangeType(0.0);
      }
    } functor;
    assemble_common(f, functor, polOrd, rhsVector);
  }  // end method

  /** if there is a first source f and a second source G:
   * discreteFunction is an output parameter (kind of return value)
   **/
  template< int polOrd, class FirstSourceType, class SecondSourceType >
  static void assemble(const FirstSourceType& f,
                const SecondSourceType& _G,
                DiscreteFunctionType& rhsVector)
  {
    struct Functor : public FunctorBase {
      const SecondSourceType& G;
      Functor(const SecondSourceType& __G)
        :G(__G) {}

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
  }  // end method

  /** if there is a first source f, a second source G and a parameter t:
   * discreteFunction is an output parameter (kind of return value)
   **/
  template< int polOrd, class FirstSourceType, class SecondSourceType >
  static void assemble(const FirstSourceType& f,
                const SecondSourceType& G,
                const TimeType& t,
                DiscreteFunctionType& rhsVector)
  {
    struct Functor : public FunctorBase {
      const SecondSourceType& G;
      const TimeType& t;
      Functor(const SecondSourceType& __G, const TimeType& _t)
        :G(__G), t(_t) {}

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
  }  // end method


  /**
   * The rhs-assemble()-methods for linear elliptic problems
   * with non-homogeneous Dirichlet and Neumann boundary conditions:
   **/

  template< int polOrd, class FirstSourceType, class DiffusionOperatorType, class NeumannBCType  >
  static void assemble(const FirstSourceType& f,
                       const DiffusionOperatorType& A,
                       const DiscreteFunctionType& dirichlet_extension, //discrete function describing dirichlet extension 
                       const NeumannBCType& neumann_bc,
                             DiscreteFunctionType& rhsVector) {
    rhsVector.clear();

    for (const auto& entity : rhsVector.space())
    {
      
      const auto& geometry = entity.geometry();
      auto elementOfRHS = rhsVector.localFunction(entity);
      const auto baseSet = rhsVector.space().basisFunctionSet(entity);

      const int numDofs = elementOfRHS.numDofs();

      std::vector<RangeType> phi_x(numDofs);
      // gradient of base function and gradient of old_u_H
      std::vector<JacobianRangeType> grad_phi_x(numDofs);
      
      const LocalFunctionType loc_dirichlet_extension = dirichlet_extension.localFunction(entity);
      const Quadrature quadrature(entity, polOrd);
      
      const auto& lagrangePointSet = rhsVector.space().lagrangePointSet( entity );

      for (const auto& intersection
         : Dune::Stuff::Common::intersectionRange(rhsVector.space().gridPart(), entity))
      {
        if ( !intersection.boundary() )
          continue;
        // boundaryId 1 = Dirichlet face; boundaryId 2 = Neumann face;
        if ( intersection.boundary() && (intersection.boundaryId() != 2) )
          continue;

        const auto face = intersection.indexInInside();
      
        const FaceQuadrature faceQuadrature( rhsVector.space().gridPart(),
                                             intersection, polOrd, FaceQuadrature::INSIDE );
        const int numFaceQuadraturePoints = faceQuadrature.nop();

        enum { faceCodim = 1 };
        for (int faceQuadraturePoint = 0; faceQuadraturePoint < numFaceQuadraturePoints; ++faceQuadraturePoint)
        {
          baseSet.evaluateAll( faceQuadrature[faceQuadraturePoint], phi_x );
          baseSet.jacobianAll( faceQuadrature[faceQuadraturePoint], grad_phi_x );

          const auto local_point_entity = faceQuadrature.point( faceQuadraturePoint ); 
          const auto global_point = geometry.global( local_point_entity ); 
          const auto local_point_face = intersection.geometry().local( global_point );

          RangeType neumann_value( 0.0 );
          neumann_bc.evaluate( global_point, neumann_value );

          const double face_weight = intersection.geometry().integrationElement( local_point_face )
                          * faceQuadrature.weight( faceQuadraturePoint );

          auto faceIterator = lagrangePointSet.template beginSubEntity< faceCodim >( face );
          const auto faceEndIterator = lagrangePointSet.template endSubEntity< faceCodim >( face );

          for ( ; faceIterator != faceEndIterator; ++faceIterator)
          {
             elementOfRHS[ *faceIterator ] += neumann_value * face_weight * phi_x[ *faceIterator ];
          }

        }

      }


      const int numQuadraturePoints = quadrature.nop();
      // the return values:
      RangeType f_x;

      JacobianRangeType gradient_dirichlet_extension;
      JacobianRangeType diffusive_flux_in_gradient_dirichlet_extension;
  
      for (int quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint)
      {
        // local (barycentric) coordinates (with respect to entity)
        const auto& local_point = quadrature.point(quadraturePoint);
        const auto global_point = geometry.global(local_point);

        const double weight = geometry.integrationElement(local_point) * quadrature.weight(quadraturePoint);
        // evaluate the Right Hand Side Function f at the current quadrature point and save its value in 'f_y':
        f.evaluate(global_point, f_x);
        // evaluate the current base function at the current quadrature point and save its value in 'z':
        baseSet.evaluateAll(quadrature[quadraturePoint], phi_x);   // i = i'te Basisfunktion;
        // evaluate the gradient of the current base function at the current quadrature point and save its value in
        // 'returnGradient':
        baseSet.jacobianAll(quadrature[quadraturePoint], grad_phi_x);
        // get gradient of dirichlet extension:
        loc_dirichlet_extension.jacobian(quadrature[quadraturePoint], gradient_dirichlet_extension );
        A.diffusiveFlux(global_point, gradient_dirichlet_extension, diffusive_flux_in_gradient_dirichlet_extension);

        for (int i = 0; i < numDofs; ++i)
        {
          elementOfRHS[i] += weight * (f_x * phi_x[i]);
          elementOfRHS[i] -= weight * (diffusive_flux_in_gradient_dirichlet_extension[0] * grad_phi_x[i][0]);
        }

      }
    }
  }  // end method

  
  /** assemble right hand side (if there is only one source - f):
   *  assemble-method for MsFEM in symmetric (non-Petrov-Galerkin) formulation
   *  rhsVector is the output parameter (kind of return value)
   **/
  template< int polOrd, class FirstSourceType, class MacroMicroGridSpecifierType, class SubGridListType >
  static void assemble_for_MsFEM_symmetric(const FirstSourceType& f,
          MacroMicroGridSpecifierType& specifier,
          SubGridListType& subgrid_list,
          DiscreteFunctionType& rhsVector) {
    // set rhsVector to zero:
    rhsVector.clear();
    const auto& coarseGridLeafIndexSet = specifier.coarseSpace().gridPart().grid().leafIndexSet();
    RangeType f_x;
    for (const auto& coarse_grid_entity : rhsVector.space()) {
      const int global_index_entity = coarseGridLeafIndexSet.index(coarse_grid_entity);

      const GeometryType& coarse_grid_geometry = coarse_grid_entity.geometry();
      auto rhsLocalFunction = rhsVector.localFunction(coarse_grid_entity);
      const int numLocalBaseFunctions = rhsLocalFunction.numDofs();

      const BasisFunctionSetType& coarse_grid_baseSet = specifier.coarseSpace().basisFunctionSet(coarse_grid_entity);

      // the sub grid U(T) that belongs to the coarse_grid_entity T
      typedef typename SubGridListType::SubGridPartType SubGridPartType;
      typedef typename SubGridListType::SubGridDiscreteFunctionSpace LocalDiscreteFunctionSpaceType;
      typedef typename SubGridListType::SubGridDiscreteFunction LocalDiscreteFunction;
      typedef typename LocalDiscreteFunction::LocalFunctionType LocalFunctionType;
      typedef Fem::CachingQuadrature< SubGridPartType, 0 > LocalGridQuadrature;

      // --------- add standard contribution of right hand side -------------------------
      {
        const Fem::CachingQuadrature< GridPartType, 0 > quadrature(coarse_grid_entity, polOrd+5);   // 0 --> codim 0
        std::vector<RangeType> phi_x_vec(numLocalBaseFunctions);
        const int numQuadraturePoints = quadrature.nop();
        for (size_t quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint)
        {
          const double det
                  = coarse_grid_geometry.integrationElement( quadrature.point(quadraturePoint) );
          // evaluate the Right Hand Side Function f at the current quadrature point and save its value in 'f_y':
          f.evaluate(coarse_grid_geometry.global( quadrature.point(quadraturePoint) ), f_x);
          coarse_grid_baseSet.evaluateAll(quadrature[quadraturePoint], phi_x_vec);
          for (int i = 0; i < numLocalBaseFunctions; ++i) {
            rhsLocalFunction[i] += det * quadrature.weight(quadraturePoint) * (f_x * phi_x_vec[i]);
          }
        }
      }
      // ----------------------------------------------------------------------------------


      // --------- add corrector contribution of right hand side --------------------------
      {
        // Load local solutions
        Multiscale::MsFEM::LocalSolutionManager localSolutionManager(coarse_grid_entity, subgrid_list, specifier);
        localSolutionManager.loadLocalSolutions();
        Multiscale::MsFEM::LocalSolutionManager::LocalSolutionVectorType& localSolutions
                = localSolutionManager.getLocalSolutions();
        assert(localSolutions.size()>0);

        for (int coarseBF = 0; coarseBF < numLocalBaseFunctions; ++coarseBF) {
          // iterator for the micro grid ( grid for the reference element T_0 )
          const auto& subGrid = subgrid_list.getSubGrid(coarse_grid_entity);
          for (const auto& localEntity : DSC::viewRange(subGrid.leafView())) {
            const auto& hostCell = subGrid.template getHostEntity< 0 >(localEntity);
            const int enclosingCoarseCellIndex = subgrid_list.getEnclosingMacroCellIndex(hostCell);
            if (enclosingCoarseCellIndex==global_index_entity) {
              //! @todo This should be defined in a macro somewhere, I just couldn't find it...
              const int order = localSolutions[0]->space().order();
              LocalGridQuadrature localQuadrature(localEntity, 2 * order + 2);
              std::vector<std::vector<typename LocalFunctionType::RangeType> >
                      allLocalSolutionEvaluations(localSolutions.size(),
                                                  std::vector<RangeType>(localQuadrature.nop(), 0.0));
              for (int lsNum=0; lsNum < localSolutions.size(); ++lsNum) {
                LocalFunctionType localFunction = localSolutions[lsNum]->localFunction(localEntity);
                localFunction.evaluateQuadrature(localQuadrature, allLocalSolutionEvaluations[lsNum]);
              }

              const auto& local_grid_geometry = localEntity.geometry();

              // higher order quadrature, since A^{\epsilon} is highly variable

              RangeType corrector_phi_x;
              for (size_t localQuadraturePoint = 0; localQuadraturePoint < localQuadrature.nop(); ++localQuadraturePoint) {
                // local (barycentric) coordinates (with respect to entity)
                const typename LocalGridQuadrature::CoordinateType& local_subgrid_point = localQuadrature.point(
                        localQuadraturePoint);
                const auto global_point_in_U_T = local_grid_geometry.global(local_subgrid_point);

                const double weight_local_quadrature
                        = localQuadrature.weight(localQuadraturePoint)
                                * local_grid_geometry.integrationElement(local_subgrid_point);

                corrector_phi_x = 0.0;

                if (specifier.simplexCoarseGrid()) {
                  assert(localSolutions.size()==GridSelector::dimgrid);

                  // transposed of the the inverse jacobian
                  const auto quadPointLocalInCoarse = coarse_grid_geometry.local(global_point_in_U_T);
                  std::vector<JacobianRangeType> gradient_Phi_vec(numLocalBaseFunctions);
                  coarse_grid_baseSet.jacobianAll(quadPointLocalInCoarse, gradient_Phi_vec);

                  //assert(localSolutionValues.size()==2);
                  corrector_phi_x += gradient_Phi_vec[coarseBF][0][0] * allLocalSolutionEvaluations[0][localQuadraturePoint];
                  corrector_phi_x += gradient_Phi_vec[coarseBF][0][1] * allLocalSolutionEvaluations[1][localQuadraturePoint];
                } else {
                  assert(localSolutions.size()==numLocalBaseFunctions);
                  // local corrector for coarse base func
                  corrector_phi_x = allLocalSolutionEvaluations[coarseBF][localQuadraturePoint];
                }
                f.evaluate(global_point_in_U_T , f_x);
                rhsLocalFunction[coarseBF] += weight_local_quadrature * (f_x * corrector_phi_x);

              }
            }
          }
        }
      }
    }
  }  // end method


  /**
   * The rhs-assemble()-methods for non-linear elliptic problems:
   * discreteFunction is an output parameter (kind of return value)
   **/
  template< int polOrd, class FirstSourceType, class DiffusionOperatorType  >
  static void assemble_for_Newton_method(const FirstSourceType& f,
                                  const DiffusionOperatorType& A,
                                  const DiscreteFunctionType& old_u_H, // old_u_H from the last iteration step
                                  DiscreteFunctionType& rhsVector) {
    rhsVector.clear();

    for (const auto& entity : rhsVector.space())
    {
      const auto& geometry = entity.geometry();
      auto elementOfRHS = rhsVector.localFunction(entity);
      const auto baseSet = rhsVector.space().basisFunctionSet(entity);

      const LocalFunctionType old_u_H_loc = old_u_H.localFunction(entity);
      const Quadrature quadrature(entity, polOrd);

      const int numDofs = elementOfRHS.numDofs(); // Dofs = Freiheitsgrade (also die Unbekannten)
      const int numQuadraturePoints = quadrature.nop();
      // the return values:
      RangeType f_x;
      std::vector<RangeType> phi_x(numDofs);
      // gradient of base function and gradient of old_u_H
      std::vector<JacobianRangeType> grad_phi_x(numDofs);
      JacobianRangeType grad_old_u_H;
      // Let A denote the diffusion operator, then we save
      // A( \gradient old_u_H )
      JacobianRangeType diffusive_flux_in_grad_old_u_H;
      for (int quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint)
      {
        // local (barycentric) coordinates (with respect to entity)
        const auto& local_point = quadrature.point(quadraturePoint);
        const auto global_point = geometry.global(local_point);
        const double det = geometry.integrationElement(local_point);
        // evaluate the Right Hand Side Function f at the current quadrature point and save its value in 'f_y':
        f.evaluate(global_point, f_x);
        // evaluate the current base function at the current quadrature point and save its value in 'z':
        baseSet.evaluateAll(quadrature[quadraturePoint], phi_x);   // i = i'te Basisfunktion;
        // evaluate the gradient of the current base function at the current quadrature point and save its value in
        // 'returnGradient':
        baseSet.jacobianAll(quadrature[quadraturePoint], grad_phi_x);
        // get gradient of old u_H:
        old_u_H_loc.jacobian(quadrature[quadraturePoint], grad_old_u_H);
        // evaluate diffusion operator in x(=global_point) and grad_old_u_H
        A.diffusiveFlux(global_point, grad_old_u_H, diffusive_flux_in_grad_old_u_H);
        for (int i = 0; i < numDofs; ++i)
        {
          elementOfRHS[i] += det * quadrature.weight(quadraturePoint) * (f_x * phi_x[i]);
          elementOfRHS[i] -= det * quadrature.weight(quadraturePoint)
                             * (diffusive_flux_in_grad_old_u_H[0] * grad_phi_x[i][0]);
        }
      }
    }
  }  // end method


  /**
   * The rhs-assemble()-methods for non-linear elliptic problems
   * if there is a first source f and a lower order term F:
   * discreteFunction is an output parameter (kind of return value)
   **/
  template< int polOrd, class FirstSourceType, class DiffusionOperatorType, class LowerOrderTermType >
  static void assemble_for_Newton_method(const FirstSourceType& f,
                                  const DiffusionOperatorType& A,
                                  const LowerOrderTermType& F,
                                  const DiscreteFunctionType& old_u_H, // old_u_H from the last iteration step
                                  DiscreteFunctionType& rhsVector) {
    rhsVector.clear();

    for (const auto& entity : rhsVector.space())
    {
      const auto& geometry = entity.geometry();
      auto elementOfRHS = rhsVector.localFunction(entity);
      const auto baseSet = rhsVector.space().basisFunctionSet(entity);

      const LocalFunctionType old_u_H_loc = old_u_H.localFunction(entity);
      const Quadrature quadrature(entity, polOrd);

      const int numDofs = elementOfRHS.numDofs(); // Dofs = Freiheitsgrade (also die Unbekannten)
      const int numQuadraturePoints = quadrature.nop();
      // the return values:
      RangeType f_x;
      std::vector<RangeType> phi_x(numDofs);
      // gradient of base function and gradient of old_u_H
      std::vector<JacobianRangeType> grad_phi_x(numDofs);
      JacobianRangeType grad_old_u_H;
      RangeType value_old_u_H;
      // Let A denote the diffusion operator, then we save
      // A( \gradient old_u_H )
      JacobianRangeType diffusive_flux_in_grad_old_u_H;
      for (int quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint)
      {
        // local (barycentric) coordinates (with respect to entity)
        const auto& local_point = quadrature.point(quadraturePoint);
        const auto global_point = geometry.global(local_point);
        const double det = geometry.integrationElement(local_point);
        // evaluate the Right Hand Side Function f at the current quadrature point and save its value in 'f_y':
        f.evaluate(global_point, f_x);
        // evaluate the current base function at the current quadrature point and save its value in 'z':
        baseSet.evaluateAll(quadrature[quadraturePoint], phi_x);   // i = i'te Basisfunktion;
        // evaluate the gradient of the current base function at the current quadrature point and save its value in
        // 'returnGradient':
        baseSet.jacobianAll(quadrature[quadraturePoint], grad_phi_x);
        // get value of old u_H:
        old_u_H_loc.evaluate(quadrature[quadraturePoint], value_old_u_H);
        // get gradient of old u_H:
        old_u_H_loc.jacobian(quadrature[quadraturePoint], grad_old_u_H);
        // evaluate diffusion operator in x(=global_point) and grad_old_u_H
        A.diffusiveFlux(global_point, grad_old_u_H, diffusive_flux_in_grad_old_u_H);

        RangeType F_x;
        F.evaluate( global_point, value_old_u_H, grad_old_u_H, F_x );

        for (int i = 0; i < numDofs; ++i)
        {
          elementOfRHS[i] += det * quadrature.weight(quadraturePoint) * (f_x * phi_x[i]);
          elementOfRHS[i] -= det * quadrature.weight(quadraturePoint)
                             * (diffusive_flux_in_grad_old_u_H[0] * grad_phi_x[i][0]);
          elementOfRHS[i] -= det * quadrature.weight(quadraturePoint) * (F_x * phi_x[i]);
        }
      }
    }
  }  // end method
  
  
  /**
   * The rhs-assemble()-methods for non-linear elliptic problems
   * if there is a first source f, a lower order term F
   * and Dirichlet and Neumann boundary conditions
   * discreteFunction is an output parameter (kind of return value)
   **/
  template< int polOrd, class FirstSourceType, class DiffusionOperatorType, class LowerOrderTermType, class NeumannBCType >
  static void assemble_for_Newton_method(const FirstSourceType& f,
                                  const DiffusionOperatorType& A,
                                  const LowerOrderTermType& F,
                                  const DiscreteFunctionType& old_u_H, // old_u_H from the last iteration step
                                  const DiscreteFunctionType& dirichlet_extension, //discrete function describing dirichlet extension 
                                  const NeumannBCType& neumann_bc,
                                  DiscreteFunctionType& rhsVector) {
    rhsVector.clear();

    for (const auto& entity : rhsVector.space())
    {
      const auto& geometry = entity.geometry();
      auto elementOfRHS = rhsVector.localFunction(entity);
      const auto baseSet = rhsVector.space().basisFunctionSet(entity);

      const int numDofs = elementOfRHS.numDofs(); // Dofs = Freiheitsgrade (also die Unbekannten)

      std::vector<RangeType> phi_x(numDofs);
      // gradient of base function and gradient of old_u_H
      std::vector<JacobianRangeType> grad_phi_x(numDofs);
      
      const LocalFunctionType old_u_H_loc = old_u_H.localFunction(entity);
      const LocalFunctionType loc_dirichlet_extension = dirichlet_extension.localFunction(entity);
      const Quadrature quadrature(entity, polOrd);

      const auto& lagrangePointSet = rhsVector.space().lagrangePointSet( entity );

      for (const auto& intersection
         : Dune::Stuff::Common::intersectionRange(rhsVector.space().gridPart(), entity))
      {
        if ( !intersection.boundary() )
          continue;
        // boundaryId 1 = Dirichlet face; boundaryId 2 = Neumann face;
        if ( intersection.boundary() && (intersection.boundaryId() != 2) )
          continue;

        const auto face = intersection.indexInInside();
      
        const FaceQuadrature faceQuadrature( rhsVector.space().gridPart(),
                                             intersection, polOrd, FaceQuadrature::INSIDE );
        const int numFaceQuadraturePoints = faceQuadrature.nop();

        enum { faceCodim = 1 };
        for (int faceQuadraturePoint = 0; faceQuadraturePoint < numFaceQuadraturePoints; ++faceQuadraturePoint)
        {
          baseSet.evaluateAll( faceQuadrature[faceQuadraturePoint], phi_x );
          baseSet.jacobianAll( faceQuadrature[faceQuadraturePoint], grad_phi_x );

          const auto local_point_entity = faceQuadrature.point( faceQuadraturePoint ); 
          const auto global_point = geometry.global( local_point_entity ); 
          const auto local_point_face = intersection.geometry().local( global_point );

          RangeType neumann_value( 0.0 );
          neumann_bc.evaluate( global_point, neumann_value );

          const double face_weight = intersection.geometry().integrationElement( local_point_face )
                          * faceQuadrature.weight( faceQuadraturePoint );

          auto faceIterator = lagrangePointSet.template beginSubEntity< faceCodim >( face );
          const auto faceEndIterator = lagrangePointSet.template endSubEntity< faceCodim >( face );

          for ( ; faceIterator != faceEndIterator; ++faceIterator)
          {
             elementOfRHS[ *faceIterator ] += neumann_value * face_weight * phi_x[ *faceIterator ];
          }

        }

      }

      const int numQuadraturePoints = quadrature.nop();
      // the return values:
      RangeType f_x;

      JacobianRangeType gradient_dirichlet_extension;
      RangeType value_dirichlet_extension;
      JacobianRangeType grad_old_u_H;
      RangeType value_old_u_H;
      // Let A denote the diffusion operator, then we save
      // A( \gradient old_u_H )
      JacobianRangeType diffusive_flux;
      
      JacobianRangeType direction;

      for (int quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint)
      {
        // local (barycentric) coordinates (with respect to entity)
        const auto& local_point = quadrature.point(quadraturePoint);
        const auto global_point = geometry.global(local_point);
        const double weight = geometry.integrationElement(local_point) * quadrature.weight(quadraturePoint);
        // evaluate the Right Hand Side Function f at the current quadrature point and save its value in 'f_y':
        f.evaluate(global_point, f_x);
        // evaluate the current base function at the current quadrature point and save its value in 'z':
        baseSet.evaluateAll(quadrature[quadraturePoint], phi_x);   // i = i'te Basisfunktion;
        // evaluate the gradient of the current base function at the current quadrature point and save its value in
        // 'returnGradient':
        baseSet.jacobianAll(quadrature[quadraturePoint], grad_phi_x);
        // get value of old u_H:
        old_u_H_loc.evaluate(quadrature[quadraturePoint], value_old_u_H);
        // get gradient of old u_H:
        old_u_H_loc.jacobian(quadrature[quadraturePoint], grad_old_u_H);
        // get value of dirichlet extension:
        loc_dirichlet_extension.evaluate(quadrature[quadraturePoint], value_dirichlet_extension );
        // get gradient of dirichlet extension:
        loc_dirichlet_extension.jacobian(quadrature[quadraturePoint], gradient_dirichlet_extension );
        direction[0] = grad_old_u_H[0] + gradient_dirichlet_extension[0];
        // evaluate diffusion operator in x(=global_point) and grad_old_u_H
        A.diffusiveFlux(global_point, direction, diffusive_flux );

        RangeType F_x;
        F.evaluate( global_point, value_old_u_H + value_dirichlet_extension, direction, F_x );

        for (int i = 0; i < numDofs; ++i)
        {
          elementOfRHS[i] += weight * (f_x * phi_x[i]);
          elementOfRHS[i] -= weight
                             * (diffusive_flux[0] * grad_phi_x[i][0]);
          elementOfRHS[i] -= weight * (F_x * phi_x[i]);
        }
      }
    }
  }  // end method

  
  //! The rhs-assemble()-methods for non-linear elliptic problems, solved with the heterogenous multiscale method
  // ( requires reconstruction of old_u_H and local fine scale averages )
  template< int polOrd, class FirstSourceType, class DiffusionOperatorType, class PeriodicDiscreteFunctionType,
            class CellProblemNumberingManagerType >
  static void assemble_for_HMM_Newton_method(const FirstSourceType& f,
                                      const DiffusionOperatorType& A,
                                      const DiscreteFunctionType& old_u_H, // old_u_H from the last iteration step
                                      // to obtain some information about the periodic discrete function space (space
                                      // for the cell problems)
                                      const CellProblemNumberingManagerType& cp_num_manager,
                                      const PeriodicDiscreteFunctionType& dummy_func,
                                      DiscreteFunctionType& rhsVector)
  {
    typedef typename PeriodicDiscreteFunctionType::LocalFunctionType
      PeriodicLocalFunctionType;

    typedef Multiscale::HMM::CellProblemSolver CellProblemSolverType;
    const std::string cell_solution_location_baseSet = "/cell_problems/_cellSolutions_baseSet";
    const std::string cell_solution_location_discFunc ="/cell_problems/_cellSolutions_discFunc";

    // reader for the cell problem data file:
    DiscreteFunctionReader discrete_function_reader_baseSet(cell_solution_location_baseSet);
    // reader for the cell problem data file:
    DiscreteFunctionReader discrete_function_reader_discFunc(cell_solution_location_discFunc);

    const double delta = DSC_CONFIG_GET("hmm.delta", 1.0f);
    const double epsilon_estimated = DSC_CONFIG_GET("hmm.epsilon_guess", 1.0f);

    const DiscreteFunctionSpaceType& discreteFunctionSpace = rhsVector.space();
    const auto& periodicDiscreteFunctionSpace = dummy_func.space();

    // set rhsVector to zero:
    rhsVector.clear();

    int number_of_entity = 0;

    const auto macro_grid_endit = discreteFunctionSpace.end();
    for (auto macro_grid_it = discreteFunctionSpace.begin(); macro_grid_it != macro_grid_endit; ++macro_grid_it)
    {
      // it* Pointer auf ein Element der Entity
      const auto& macro_grid_geometry = (*macro_grid_it).geometry(); // Referenz auf Geometrie
      auto elementOfRHS = rhsVector.localFunction(*macro_grid_it);

      const BasisFunctionSetType macro_grid_baseSet
        = discreteFunctionSpace.basisFunctionSet(*macro_grid_it);
      const LocalFunctionType old_u_H_loc = old_u_H.localFunction(*macro_grid_it);
      // for \int_{\Omega} f \Phi
      const Quadrature macro_quadrature(*macro_grid_it, polOrd);
      // for - \int_{\Omega} \in_Y A^{\epsilon}( gradient reconstruction ) \nabla \Phi
      const Quadrature one_point_macro_quadrature(*macro_grid_it, 0);
      // the fine scale reconstructions are only available for the barycenter of the macro grid entity (=> only
      // available for the canonical one point quadrature on this element)
      const auto& local_macro_point = one_point_macro_quadrature.point(0 /*=quadraturePoint*/);
      // barycenter of macro grid entity
      const auto macro_entity_barycenter = macro_grid_geometry.global(local_macro_point);

      const double macro_entity_volume = one_point_macro_quadrature.weight(0 /*=quadraturePoint*/)
                                         * macro_grid_geometry.integrationElement(local_macro_point);

      const int numDofs = elementOfRHS.numDofs(); // Dofs = Freiheitsgrade
      // gradient of base function and gradient of old_u_H
      std::vector<JacobianRangeType> grad_Phi_x_vec(numDofs);
      std::vector<RangeType> phi_x(numDofs);
      JacobianRangeType grad_old_u_H_x;

      // get gradient of old u_H:
      old_u_H_loc.jacobian(one_point_macro_quadrature[0], grad_old_u_H_x);

      // Q_h(u_H^{(n-1}))(x_T,y):
      PeriodicDiscreteFunctionType corrector_old_u_H("Corrector of u_H^(n-1)", periodicDiscreteFunctionSpace);
      corrector_old_u_H.clear();

      PeriodicDiscreteFunctionType corrector_Phi_i("Corrector of Phi_i", periodicDiscreteFunctionSpace);
      discrete_function_reader_discFunc.read(number_of_entity, corrector_old_u_H);
      macro_grid_baseSet.jacobianAll(one_point_macro_quadrature[0], grad_Phi_x_vec);

      for (int i = 0; i < numDofs; ++i)
      {
        const auto& grad_Phi_x = grad_Phi_x_vec[i];
        // --------------- the source contribution ( \int_{\Omega} f \Phi ) -------------------------------

        // the return values:
        RangeType f_x;

        const int numMacroQuadraturePoints = macro_quadrature.nop();
        for (int quadraturePoint = 0; quadraturePoint < numMacroQuadraturePoints; ++quadraturePoint)
        {
          // local (barycentric) coordinates (with respect to entity)
          const auto& local_point = macro_quadrature.point(quadraturePoint);
          const auto global_point = macro_grid_geometry.global(local_point);
          const double quad_weight
            = macro_grid_geometry.integrationElement(local_point) * macro_quadrature.weight(quadraturePoint);
          // evaluate the Right Hand Side Function f at the current quadrature point and save its value in 'f_y':
          f.evaluate(global_point, f_x);
          //!TODO order of loops sucks
          macro_grid_baseSet.evaluateAll(macro_quadrature[quadraturePoint], phi_x);
          elementOfRHS[i] += quad_weight * (f_x * phi_x[i]);
        }

        // --------------- end of source contribution -----------------------------------------------------

        // --------------- the contribution of the jacobian of the diffusion operator, evaluated in the old
        // reconstructed macro solution -------------------------------


        corrector_Phi_i.clear();
        if ( !DSC_CONFIG_GET("hmm.petrov_galerkin", true ) ) {
            typename EntityType::EntityPointer entity_ptr(macro_grid_it);
            discrete_function_reader_baseSet.read(cp_num_manager.get_number_of_cell_problem(entity_ptr, i)
                                                  , corrector_Phi_i);
        }

        RangeType fine_scale_contribution = 0.0;
        for (const auto& micro_grid_entity : periodicDiscreteFunctionSpace)
        {
          const GeometryType& micro_grid_geometry = micro_grid_entity.geometry();
          assert(micro_grid_entity.partitionType() == InteriorEntity);

          auto loc_corrector_old_u_H = corrector_old_u_H.localFunction(micro_grid_entity);
          auto loc_corrector_Phi_i = corrector_Phi_i.localFunction(micro_grid_entity);

          // higher order quadrature, since A^{\epsilon} is highly variable
          const Quadrature micro_grid_quadrature(micro_grid_entity, 2 * periodicDiscreteFunctionSpace.order() + 2);
          const size_t numQuadraturePoints = micro_grid_quadrature.nop();

          for (size_t microQuadraturePoint = 0; microQuadraturePoint < numQuadraturePoints; ++microQuadraturePoint)
          {
            // local (barycentric) coordinates (with respect to entity)
            const typename Quadrature::CoordinateType& local_micro_point = micro_grid_quadrature.point(
              microQuadraturePoint);

            const DomainType global_point_in_Y = micro_grid_geometry.global(local_micro_point);

            const double weight_micro_quadrature = micro_grid_quadrature.weight(microQuadraturePoint)
                                                   * micro_grid_geometry.integrationElement(local_micro_point);

            JacobianRangeType grad_corrector_old_u_H;
            loc_corrector_old_u_H.jacobian(micro_grid_quadrature[microQuadraturePoint], grad_corrector_old_u_H);

            JacobianRangeType grad_corrector_Phi_i;
            if ( !DSC_CONFIG_GET("hmm.petrov_galerkin", true ) )
              loc_corrector_Phi_i.jacobian(micro_grid_quadrature[microQuadraturePoint], grad_corrector_Phi_i);

            // x_T + (delta * y)
            DomainType current_point_in_macro_grid;
            for (int k = 0; k < dimension; ++k)
              current_point_in_macro_grid[k] = macro_entity_barycenter[k] + (delta * global_point_in_Y[k]);

            // evaluate jacobian matrix of diffusion operator in 'position_vector' in direction 'direction_vector':

            JacobianRangeType direction_vector;
            for (int k = 0; k < dimension; ++k)
              direction_vector[0][k] = grad_old_u_H_x[0][k] + grad_corrector_old_u_H[0][k];

            JacobianRangeType diffusive_flux;
            A.diffusiveFlux(current_point_in_macro_grid, direction_vector, diffusive_flux);

            double cutting_function = 1.0;
            for (int k = 0; k < dimension; ++k)
            {
              // is the current quadrature point in the relevant cell?
              if ( fabs(global_point_in_Y[k]) > ( 0.5 * (epsilon_estimated / delta) ) )
              { cutting_function *= 0.0; }
            }

            // if test function reconstruction = non-Petrov-Galerkin
            if ( !DSC_CONFIG_GET("hmm.petrov_galerkin", true ) ) {
              JacobianRangeType grad_reconstruction_Phi_i;
              for (int k = 0; k < dimension; ++k)
                grad_reconstruction_Phi_i[0][k] = grad_Phi_x[0][k] + grad_corrector_Phi_i[0][k];

              fine_scale_contribution += cutting_function * weight_micro_quadrature
                                         * (diffusive_flux[0] * grad_reconstruction_Phi_i[0]);
            } else {
              fine_scale_contribution += cutting_function * weight_micro_quadrature * (diffusive_flux[0] * grad_Phi_x[0]);
            }
          }
        }

        elementOfRHS[i] -= pow(delta / epsilon_estimated, dimension) * macro_entity_volume * fine_scale_contribution;

        // --------------- end of diffusion contribution -----------------------------------------------------
      }

      number_of_entity += 1;
    }
  }  // end method
}; // end class
} // end namespace

#endif // ifndef DUNE_RIGHT_HAND_SIDE_ASSEMBLER_HH
