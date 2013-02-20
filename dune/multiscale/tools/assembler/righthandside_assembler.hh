#ifndef DUNE_RIGHT_HAND_SIDE_ASSEMBLER_HH
#define DUNE_RIGHT_HAND_SIDE_ASSEMBLER_HH

// - Dune includes
#include <dune/fem/quadrature/quadrature.hh>
#include <dune/multiscale/tools/disc_func_writer/discretefunctionwriter.hh>
#include <dune/multiscale/tools/solver/HMM/cell_problem_solving/solver.hh>

//#include <dune/multiscale/tools/solver/MsFEM/msfem_localproblems/subgrid-list.hh>
#include <dune/multiscale/tools/misc.hh>
#include <dune/stuff/fem/functions/checks.hh>

namespace Dune {
// Assembler for right rand side
// We assemble the right hand side in a LSE, i.e. f \cdot \Phi_H + G \cdot \nabala \Phi_H
// we call f the first Source and G the second Source
template< class DiscreteFunctionImp >
class RightHandSideAssembler
{
private:
  typedef DiscreteFunctionImp DiscreteFunctionType;
  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType
    DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionType::LocalFunctionType
    LocalFunctionType;
  typedef typename DiscreteFunctionSpaceType::BaseFunctionSetType
    BaseFunctionSetType;
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
  typedef CachingQuadrature< GridPartType, 0 > Quadrature;
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

      const BaseFunctionSetType baseSet // BaseFunctions leben immer auf Refernzelement!!!
        = rhsVector.space().baseFunctionSet(entity);     // entity Referenz auf eine bestimmtes Element der entity. In der
                                                          // ersten Klasse war das Element fest, deshalb konnte man sich
                                                          // dort Pointer sparen. //loeschen: discreteFunctionSpace
                                                          // statt
                                                          // functionSpace

      const CachingQuadrature< GridPartType, 0 > quadrature(entity, polOrd);   // 0 --> codim 0
      const int numQuadraturePoints = quadrature.nop();
      const int numDofs = elementOfRHS.numDofs();
      for (int quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint)
      {
        // the return values:
        RangeType f_x;
        std::vector<RangeType>phi_x;
        std::vector<JacobianRangeType> gradientPhi;
        // to save: A \nabla PHI_H * \nabla phi_h;
        RangeType res = 0;
        const double det = geometry.integrationElement(quadrature.point(quadraturePoint));
        const auto& inv = geometry.jacobianInverseTransposed(quadrature.point(quadraturePoint));
        f.evaluate(geometry.global( quadrature.point(quadraturePoint) ), f_x);
        baseSet.evaluateAll(quadrature[quadraturePoint], phi_x);
        baseSet.jacobianAll(quadrature[quadraturePoint], inv, gradientPhi);

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
                DiscreteFunctionType& rhsVector) {
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
                DiscreteFunctionType& rhsVector)  {
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

  // /############################


  /** assemble right hand side (if there is only one source - f):
   *  assemble-method for MsFEM in symmetric (non-Petrov-Galerkin) formulation 
   *  rhsVector is the output parameter (kind of return value)
   **/
  template< int polOrd, class FirstSourceType, class MacroMicroGridSpecifierType, class SubGridListType >
  static void assemble_for_MsFEM_symmetric(const FirstSourceType& f, MacroMicroGridSpecifierType& specifier, SubGridListType& subgrid_list,
                                           DiscreteFunctionType& rhsVector)
  {
    // set rhsVector to zero:
    rhsVector.clear();
    
    const auto& coarseGridLeafIndexSet = specifier.coarseSpace().gridPart().grid().leafIndexSet();
    for (const auto& coarse_grid_entity : rhsVector.space())
    {

      const int global_index_entity = coarseGridLeafIndexSet.index(coarse_grid_entity);
    
      const GeometryType& coarse_grid_geometry = coarse_grid_entity.geometry();
      auto elementOfRHS = rhsVector.localFunction(coarse_grid_entity);

      const BaseFunctionSetType coarse_grid_baseSet = specifier.coarseSpace().baseFunctionSet(coarse_grid_entity);

      // the sub grid U(T) that belongs to the coarse_grid_entity T
      typedef typename SubGridListType::SubGridPartType SubGridPartType;
      typedef typename SubGridListType::SubGridDiscreteFunctionSpace LocalDiscreteFunctionSpace;
      typedef typename SubGridListType::SubGridDiscreteFunction LocalDiscreteFunction;
      typedef CachingQuadrature< SubGridPartType, 0 > LocalGridQuadrature;
      
      auto& sub_grid_U_T = subgrid_list.getSubGrid(global_index_entity);
      SubGridPartType subGridPart(sub_grid_U_T);      
      LocalDiscreteFunctionSpace localDiscreteFunctionSpace(subGridPart);

      LocalDiscreteFunction local_problem_solution_e0("Local problem Solution e_0", localDiscreteFunctionSpace);
      local_problem_solution_e0.clear();

      LocalDiscreteFunction local_problem_solution_e1("Local problem Solution e_1", localDiscreteFunctionSpace);
      local_problem_solution_e1.clear();
      
      // -- load local solutions --
      // the file/place, where we saved the solutions of the cell problems
      const std::string local_solution_location = (boost::format("local_problems/_localProblemSolutions_%d")
                                                  % global_index_entity).str();
      // reader for the cell problem data file:
      DiscreteFunctionReader discrete_function_reader(local_solution_location);
      discrete_function_reader.read(0, local_problem_solution_e0);
      discrete_function_reader.read(1, local_problem_solution_e1);
      
      // --------- add standard contribution of right hand side -------------------------
      const CachingQuadrature< GridPartType, 0 > quadrature(coarse_grid_entity, polOrd+5);   // 0 --> codim 0
      const int numDofs = elementOfRHS.numDofs();
      std::vector<RangeType> phi_x_vec(numDofs);
      RangeType f_x;
      const int numQuadraturePoints = quadrature.nop();
      for (int quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint)
      {
        const double det
          = coarse_grid_geometry.integrationElement( quadrature.point(quadraturePoint) );
        // evaluate the Right Hand Side Function f at the current quadrature point and save its value in 'f_y':
        f.evaluate(coarse_grid_geometry.global( quadrature.point(quadraturePoint) ), f_x);
        coarse_grid_baseSet.evaluateAll(quadrature[quadraturePoint], phi_x_vec);
        for (int i = 0; i < numDofs; ++i)
        {
          elementOfRHS[i] += det * quadrature.weight(quadraturePoint) * (f_x * phi_x_vec[i]);
        }
      }
      // ----------------------------------------------------------------------------------
      

      // --------- add corrector contribution of right hand side --------------------------
      
      // 1 point quadrature!! We only need the gradient of the base function,
      // which is constant on the whole entity.
      const CachingQuadrature< GridPartType, 0 > one_point_quadrature(coarse_grid_entity, 0);

      // the barycenter of the macro_grid_entity
      const typename CachingQuadrature< GridPartType, 0 >::CoordinateType& local_coarse_point
           = one_point_quadrature.point(0 /*=quadraturePoint*/);

      // transposed of the the inverse jacobian
      const auto& inverse_jac = coarse_grid_geometry.jacobianInverseTransposed(local_coarse_point);
      std::vector<JacobianRangeType> gradient_Phi_vec(numDofs);
      coarse_grid_baseSet.jacobianAll(one_point_quadrature[0], inverse_jac, gradient_Phi_vec);

      for (int i = 0; i < numDofs; ++i)
      {
        // iterator for the micro grid ( grid for the reference element T_0 )
        for (const auto& local_grid_entity : localDiscreteFunctionSpace)
        {
          // check if "local_grid_entity" (which is an entity of U(T)) is in T:
          // -------------------------------------------------------------------

          const auto father_of_loc_grid_ent =
              Stuff::Grid::make_father(coarseGridLeafIndexSet,
                                       localDiscreteFunctionSpace.grid().template getHostEntity< 0 >(local_grid_entity),
                                       specifier.getLevelDifference());
          if (!Stuff::Grid::entities_identical(coarse_grid_entity, *father_of_loc_grid_ent))
          {
            continue;
          }
          // -------------------------------------------------------------------
          
          const auto& local_grid_geometry = local_grid_entity.geometry();
          assert(local_grid_entity.partitionType() == InteriorEntity);

          // higher order quadrature, since A^{\epsilon} is highly variable
          LocalGridQuadrature local_grid_quadrature(local_grid_entity, 2 * localDiscreteFunctionSpace.order() + 2);
          for (size_t localQuadraturePoint = 0; localQuadraturePoint < local_grid_quadrature.nop(); ++localQuadraturePoint)
          {
            RangeType corrector_phi_x;

            // local (barycentric) coordinates (with respect to entity)
            const typename LocalGridQuadrature::CoordinateType& local_subgrid_point = local_grid_quadrature.point(
                                                                                        localQuadraturePoint);
            const auto global_point_in_U_T = local_grid_geometry.global(local_subgrid_point);

            const double weight_local_quadrature
                = local_grid_quadrature.weight(localQuadraturePoint)
                  * local_grid_geometry.integrationElement(local_subgrid_point);

            const auto localized_local_problem_solution_e0 = local_problem_solution_e0.localFunction(
                                                               local_grid_entity);
            const auto localized_local_problem_solution_e1 = local_problem_solution_e1.localFunction(
                                                               local_grid_entity);

            // local corrector for e_0 and e_1
            RangeType loc_sol_e0, loc_sol_e1;
            localized_local_problem_solution_e0.evaluate(local_grid_quadrature[localQuadraturePoint], loc_sol_e0);
            localized_local_problem_solution_e1.evaluate(local_grid_quadrature[localQuadraturePoint], loc_sol_e1);

            corrector_phi_x = 0.0;
            corrector_phi_x += gradient_Phi_vec[i][0][0] * loc_sol_e0;
            corrector_phi_x += gradient_Phi_vec[i][0][1] * loc_sol_e1;

            f.evaluate( global_point_in_U_T , f_x);

            elementOfRHS[i] += weight_local_quadrature * (f_x * corrector_phi_x);
          }
        }
      }
    }
  }  // end method

 
  /**
   * The rhs-assemble()-methods for non-linear elliptic problems
   * if there is a first source f and a second source G:
   * discreteFunction is an output parameter (kind of return value)
   **/
  template< int polOrd, class FirstSourceType, class DiffusionOperatorType >
  static void assemble_for_Newton_method(const FirstSourceType& f,
                                  const DiffusionOperatorType& A,
                                  const DiscreteFunctionType& old_u_H, // old_u_H from the last iteration step
                                  DiscreteFunctionType& rhsVector) {
    rhsVector.clear();

    for (const auto& entity : rhsVector.space())
    {
      const auto& geometry = entity.geometry();
      auto elementOfRHS = rhsVector.localFunction(entity);
      const auto baseSet = rhsVector.space().baseFunctionSet(entity);

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
        const auto& inv = geometry.jacobianInverseTransposed(local_point);
        // evaluate the Right Hand Side Function f at the current quadrature point and save its value in 'f_y':
        f.evaluate(global_point, f_x);
        // evaluate the current base function at the current quadrature point and save its value in 'z':
        baseSet.evaluateAll(quadrature[quadraturePoint], phi_x);   // i = i'te Basisfunktion;
        // evaluate the gradient of the current base function at the current quadrature point and save its value in
        // 'returnGradient':
        baseSet.jacobianAll(quadrature[quadraturePoint], inv, grad_phi_x);
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

    typedef CellProblemSolver< PeriodicDiscreteFunctionType, DiffusionOperatorType > CellProblemSolverType;
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

      const BaseFunctionSetType macro_grid_baseSet
        = discreteFunctionSpace.baseFunctionSet(*macro_grid_it);
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
      JacobianRangeType grad_old_u_H_x;

      const auto& inv = macro_grid_geometry.jacobianInverseTransposed(local_macro_point);
      // get gradient of old u_H:
      old_u_H_loc.jacobian(one_point_macro_quadrature[0], grad_old_u_H_x);

      // Q_h(u_H^{(n-1}))(x_T,y):
      PeriodicDiscreteFunctionType corrector_old_u_H("Corrector of u_H^(n-1)", periodicDiscreteFunctionSpace);
      corrector_old_u_H.clear();

      PeriodicDiscreteFunctionType corrector_Phi_i("Corrector of Phi_i", periodicDiscreteFunctionSpace);
      /* // if the cell problems are not precomputed, we might use:
        CellProblemSolverType cell_problem_solver(periodicDiscreteFunctionSpace, A);
        cell_problem_solver.template solvecellproblem< JacobianRangeType >
          (grad_old_u_H_x, macro_entity_barycenter, corrector_old_u_H);
        if (DSC_CONFIG_GET("TFR", false))
          cell_problem_solver.template solvecellproblem< JacobianRangeType >
            (grad_Phi_x, macro_entity_barycenter, corrector_Phi_i);

      */
      discrete_function_reader_discFunc.read(number_of_entity, corrector_old_u_H);
      macro_grid_baseSet.jacobianAll(one_point_macro_quadrature[0], inv, grad_Phi_x_vec);


      for (int i = 0; i < numDofs; ++i)
      {
        const auto& grad_Phi_x = grad_Phi_x_vec[i];
        // --------------- the source contribution ( \int_{\Omega} f \Phi ) -------------------------------

        // the return values:
        RangeType f_x, phi_x;

        const int numMacroQuadraturePoints = macro_quadrature.nop();
        for (int quadraturePoint = 0; quadraturePoint < numMacroQuadraturePoints; ++quadraturePoint)
        {
          // local (barycentric) coordinates (with respect to entity)
          const typename Quadrature::CoordinateType& local_point = macro_quadrature.point(quadraturePoint);

          const auto global_point = macro_grid_geometry.global(local_point);

          const double quad_weight
            = macro_grid_geometry.integrationElement(local_point) * macro_quadrature.weight(quadraturePoint);

          // evaluate the Right Hand Side Function f at the current quadrature point and save its value in 'f_y':
          f.evaluate(global_point, f_x);

          // evaluate the current base function at the current quadrature point and save its value in 'z':
          macro_grid_baseSet.evaluate(i, macro_quadrature[quadraturePoint], phi_x);   // i = i'te Basisfunktion;

          elementOfRHS[i] += quad_weight * (f_x * phi_x);
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
