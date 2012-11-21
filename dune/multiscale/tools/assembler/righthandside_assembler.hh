#ifndef DUNE_RIGHT_HAND_SIDE_ASSEMBLER_HH
#define DUNE_RIGHT_HAND_SIDE_ASSEMBLER_HH

// - Dune includes
#include <dune/fem/quadrature/quadrature.hh>
#include <dune/multiscale/tools/disc_func_writer/discretefunctionwriter.hh>
#include <dune/multiscale/tools/solver/HMM/cell_problem_solving/solver.hh>


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

  // /############################ The rhs-assemble()-methods for linear elliptic problems
  // #########################################

private:
  //! need a virtual base to work around local classes not being allowed in templated scopes
  struct FunctorBase {
    virtual RangeType operator()(const DomainType& global_quad_point, const JacobianRangeType& gradientPhi) const = 0;
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
      // it* Pointer auf ein Element der Entity
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
      const int numDofs = elementOfRHS.numDofs(); // Dofs = Freiheitsgrade (also die Unbekannten)
      for (int i = 0; i < numDofs; ++i)  // Laufe ueber alle Knoten des entity-elements auf dem wir uns befinden
      {
        // the return values:
        RangeType f_x, phi_x;
        JacobianRangeType gradientPhi;
        // to save: A \nabla PHI_H * \nabla phi_h;
        RangeType res = 0;

        const int numQuadraturePoints = quadrature.nop();
        for (int quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint)
        {
          const double det
            = geometry.integrationElement( quadrature.point(quadraturePoint) );
          // evaluate the Right Hand Side Function f at the current quadrature point and save its value in 'f_y':
          f.evaluate(geometry.global( quadrature.point(quadraturePoint) ), f_x);
          // evaluate the current base function at the current quadrature point and save its value in 'z':
          baseSet.evaluate(i, quadrature[quadraturePoint], phi_x);   // i = i'te Basisfunktion;
          // evaluate the gradient of the current base function at the current quadrature point and save its value in
          // 'returnGradient':
          baseSet.jacobian(i, quadrature[quadraturePoint], gradientPhi);
          // Bis jetzt nur der Gradient auf dem Referenzelement!!!!!!! Es muss noch transformiert werden, um den
          // Gradienten auf dem echten Element zu bekommen! Das passiert folgendermassen:
          const FieldMatrix< double, dimension, dimension >& inv
            = geometry.jacobianInverseTransposed( quadrature.point(quadraturePoint) );
          // multiply with transpose of jacobian inverse
          gradientPhi[0] = FMatrixHelp::mult(inv, gradientPhi[0]);

          res = functor(geometry.global(quadrature.point(quadraturePoint)), gradientPhi);

          elementOfRHS[i] += det * quadrature.weight(quadraturePoint) * (f_x * phi_x);
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
  // #########################################

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
      const GeometryType& geometry = entity.geometry();
      LocalFunctionType elementOfRHS = rhsVector.localFunction(entity);
      const BaseFunctionSetType baseSet // BaseFunctions leben immer auf Refernzelement!!!
        = rhsVector.space().baseFunctionSet(entity);

      const LocalFunctionType old_u_H_loc = old_u_H.localFunction(entity);
      const Quadrature quadrature(entity, polOrd);

      const int numDofs = elementOfRHS.numDofs(); // Dofs = Freiheitsgrade (also die Unbekannten)
      for (int i = 0; i < numDofs; ++i)  // Laufe ueber alle Knoten des entity-elements auf dem wir uns befinden
      {
        // the return values:
        RangeType f_x, phi_x;
        // gradient of base function and gradient of old_u_H
        JacobianRangeType grad_phi_x, grad_old_u_H;
        // Let A denote the diffusion operator, then we save
        // A( \gradient old_u_H )
        JacobianRangeType diffusive_flux_in_grad_old_u_H;

        const int numQuadraturePoints = quadrature.nop();
        for (int quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint)
        {
          // local (barycentric) coordinates (with respect to entity)
          const auto& local_point = quadrature.point(quadraturePoint);
          const DomainType global_point = geometry.global(local_point);

          const double det = geometry.integrationElement(local_point);
          // evaluate the Right Hand Side Function f at the current quadrature point and save its value in 'f_y':
          f.evaluate(global_point, f_x);
          // evaluate the current base function at the current quadrature point and save its value in 'z':
          baseSet.evaluate(i, quadrature[quadraturePoint], phi_x);   // i = i'te Basisfunktion;
          // evaluate the gradient of the current base function at the current quadrature point and save its value in
          // 'returnGradient':
          baseSet.jacobian(i, quadrature[quadraturePoint], grad_phi_x);
          // Bis jetzt nur der Gradient auf dem Referenzelement!!!!!!! Es muss noch transformiert werden, um den
          // Gradienten auf dem echten Element zu bekommen! Das passiert folgendermassen:
          const FieldMatrix< double, dimension, dimension >& inv
            = geometry.jacobianInverseTransposed(local_point);
          // multiply with transpose of jacobian inverse
          grad_phi_x[0] = FMatrixHelp::mult(inv, grad_phi_x[0]);
          // get gradient of old u_H:
          old_u_H_loc.jacobian(quadrature[quadraturePoint], grad_old_u_H);
          // evaluate diffusion operator in x(=global_point) and grad_old_u_H
          A.diffusiveFlux(global_point, grad_old_u_H, diffusive_flux_in_grad_old_u_H);
          elementOfRHS[i] += det * quadrature.weight(quadraturePoint) * (f_x * phi_x);
          elementOfRHS[i] -= det * quadrature.weight(quadraturePoint)
                             * (diffusive_flux_in_grad_old_u_H[0] * grad_phi_x[0]);
        }
      }
    }
  }  // end method

  /** /############################ The rhs-assemble()-method for linear elliptic problems, solved with MsFEM in
   * non-Petriv-Galerkin-Formulation #########################################
   * assemble standard right hand side:
   * if there is only one source (f) (there is no second source):
   * discreteFunction is an output parameter (kind of return value)
   **/
  template< int polOrd, class FirstSourceType, class LocalProblemNumberingManagerType >
  static void assemble_msfem(// get number of local problem to determine the reconstruction
                      const LocalProblemNumberingManagerType& lp_num_manager,
                      const FirstSourceType& f,
                      DiscreteFunctionType& rhsVector) {
    // get the local discrete function space
    const DiscreteFunctionSpaceType& localDiscreteFunctionSpace = lp_num_manager.get_local_discrete_function_space();
    // reader for the local problem data file:
    DiscreteFunctionReader discrete_function_reader(lp_num_manager.get_location());

    // set discreteFunction to zero:
    rhsVector.clear();
    for (const auto& entity : rhsVector.space())
    {
      const GeometryType& geometry = (entity).geometry();
      LocalFunctionType elementOfRHS = rhsVector.localFunction(entity);

      const BaseFunctionSetType& baseSet // BaseFunctions leben immer auf Refernzelement!!!
        = rhsVector.space().baseFunctionSet(entity);
      const CachingQuadrature< GridPartType, 0 > quadrature(entity, polOrd);

      // transformation F : T_0 -> T
      // describe the mapping F(x) = Ax + b with F(T_0)=T for an entity T and the reference element T_0:
      // arguments: entity T, point in T_0, point in T.

      // Let (a_0,a_1,a_2) deonte the corners of the 2-simplex T, then the matrix A in the affine transformation
      // F(x) = Ax + a_0, F : T_0 -> T is given by
      // A_11 = a_1( 1 ) - a_0( 1 )     A_12 = a_2( 1 ) - a_0( 1 )
      // A_21 = a_1( 2 ) - a_0( 2 )     A_22 = a_2( 2 ) - a_0( 2 )

      // corners of the reference element:
      typedef typename CachingQuadrature< GridPartType, 0 >::CoordinateType CoordinateType;

      const CoordinateType ref_corner_0 = { 0.0, 0.0 };
      const CoordinateType ref_corner_1 = { 0.0, 1.0 };
      const CoordinateType ref_corner_2 = { 1.0, 0.0 };

      // corner of the global element:
      const DomainType corner_0_of_T = geometry.global(ref_corner_0);
      const DomainType corner_1_of_T = geometry.global(ref_corner_1);
      const DomainType corner_2_of_T = geometry.global(ref_corner_2);

      // value of the matrix A (in F(x) = Ax + a_0)
      double val_A[dimension][dimension];
      val_A[0][0] = corner_1_of_T[0] - corner_0_of_T[0];
      val_A[0][1] = corner_2_of_T[0] - corner_0_of_T[0];
      val_A[1][0] = corner_1_of_T[1] - corner_0_of_T[1];
      val_A[1][1] = corner_2_of_T[1] - corner_0_of_T[1];

      // define 'c := (a_1(1) - a_0(1))·(a_2(2) - a_0(2)) - (a_1(2) - a_0(2))·(a_2(1) - a_0(1))
      const double c = 1.0 / ( (val_A[0][0] * val_A[1][1]) - (val_A[0][1] * val_A[1][0]) );

      // |det(A)|:
      const double abs_det_A = fabs(1.0 / c);

      const int numDofs = elementOfRHS.numDofs(); // Dofs = Freiheitsgrade (also die Unbekannten)

      for (int i = 0; i < numDofs; ++i)  // Laufe ueber alle Knoten des entity-elements auf dem wir uns befinden
      {
        // the return values:
        RangeType f_x, phi_x;
        const int numQuadraturePoints = quadrature.nop();
        for (int quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint)
        {
          const double det
            = geometry.integrationElement( quadrature.point(quadraturePoint) );

          // evaluate the Right Hand Side Function f at the current quadrature point and save its value in 'y':
          f.evaluate(geometry.global( quadrature.point(quadraturePoint) ), f_x);
          // evaluate the current base function at the current quadrature point and save its value in 'z':
          baseSet.evaluate(i, quadrature[quadraturePoint], phi_x);   // i = i'te Basisfunktion;
          elementOfRHS[i] += det * quadrature.weight(quadraturePoint) * (f_x * phi_x);
        }

        // If we do not use the Petrov-Galerkin formulation,
        // we need to add the contribution of the corrector
        // (we splitt:  \int_T f ( \phi + Q(\phi) ) = \int_T f \phi + \int_{T_0} (f ○ F) (Q(\phi)○F) |det A|
        // where we already added \int_T f \phi in the previous step )
        if ( !DSC_CONFIG_GET("PGF", false) )
        {
          // get number of cell problem from entity and number of base function
          const int cell_problem_id = lp_num_manager.get_number_of_local_problem(entity, i);

          DiscreteFunctionType corrector_phi("Corrector Function of Phi", localDiscreteFunctionSpace);
          corrector_phi.clear();

          discrete_function_reader.read(cell_problem_id, corrector_phi);

          RangeType f_x_transformed, Q_phi_transformed;
          // iterator for the micro grid ( grid for the reference element T_0 )
          for (const auto& micro_grid_entity : localDiscreteFunctionSpace)
          {
            const GeometryType& micro_grid_geometry = micro_grid_entity.geometry();

            // ( Q^eps(\Phi) ○ F ):
            typename DiscreteFunctionType::LocalFunctionType localized_corrector = corrector_phi.localFunction(
              micro_grid_entity);

            // higher order quadrature
            const Quadrature micro_grid_quadrature(micro_grid_entity, 2 * localDiscreteFunctionSpace.order() + 2);
            const size_t locNumQuadraturePoints = micro_grid_quadrature.nop();

            for (size_t microQuadraturePoint = 0; microQuadraturePoint < locNumQuadraturePoints; ++microQuadraturePoint)
            {
              // local (barycentric) coordinates (with respect to entity)
              const typename Quadrature::CoordinateType& local_micro_point = micro_grid_quadrature.point(
                microQuadraturePoint);

              const DomainType global_point_in_T_0 = micro_grid_geometry.global(local_micro_point);

              const double weight_micro_quadrature = micro_grid_quadrature.weight(microQuadraturePoint)
                                               * micro_grid_geometry.integrationElement(local_micro_point)
                                                     * abs_det_A;

              // Q^eps(\phi) ○ F :
              localized_corrector.evaluate(micro_grid_quadrature[microQuadraturePoint], Q_phi_transformed);

              // 'F(x)', i.e. F ( global point in the reference element T_0 )
              // (the transformation of the global point in T_0 to its position in T)
              DomainType global_point_transformed(0.0);

              for (int k = 0; k < dimension; ++k)
                for (int l = 0; l < dimension; ++l)
                  global_point_transformed[k] += (val_A[k][l] * global_point_in_T_0[l]);

              global_point_transformed += corner_0_of_T;

              // F(x) = Ax + a_0, F : T_0 -> T is given by
              // A_11 = a_1(1) - a_0(1)     A_12 = a_2(1) - a_0(1)
              // A_21 = a_1(2) - a_0(2)     A_22 = a_2(2) - a_0(2)
              f.evaluate(global_point_transformed, f_x_transformed);

              // note that |det A| is already contained in 'weight_micro_quadrature'
              elementOfRHS[i] += weight_micro_quadrature * f_x_transformed * Q_phi_transformed;
            }
          }
        } //#endif // ifndef PGF
      }
    }
  }  // end method

  /** /############################ The rhs-assemble()-methods for non-linear elliptic problems, solved with the
   * heterogenous multiscale method  #########################################
   * requires reconstruction of old_u_H and local fine scale averages
   **/
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
    typedef typename PeriodicDiscreteFunctionType::DiscreteFunctionSpaceType
      PeriodicDiscreteFunctionSpaceType;

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

    const DiscreteFunctionSpaceType& discreteFunctionSpace
      = rhsVector.space();
    const PeriodicDiscreteFunctionSpaceType& periodicDiscreteFunctionSpace
      = dummy_func.space();

    // set rhsVector to zero:
    rhsVector.clear();

    int number_of_entity = 0;

    const auto macro_grid_endit = discreteFunctionSpace.end();
    for (auto macro_grid_it = discreteFunctionSpace.begin(); macro_grid_it != macro_grid_endit; ++macro_grid_it)
    {
      // it* Pointer auf ein Element der Entity
      const GeometryType& macro_grid_geometry = (*macro_grid_it).geometry(); // Referenz auf Geometrie

      LocalFunctionType elementOfRHS = rhsVector.localFunction(*macro_grid_it);

      const BaseFunctionSetType macro_grid_baseSet
        = discreteFunctionSpace.baseFunctionSet(*macro_grid_it);

      const LocalFunctionType old_u_H_loc = old_u_H.localFunction(*macro_grid_it);

      // for \int_{\Omega} f \Phi
      const Quadrature macro_quadrature(*macro_grid_it, polOrd);

      // for - \int_{\Omega} \in_Y A^{\epsilon}( gradient reconstruction ) \nabla \Phi
      const Quadrature one_point_macro_quadrature(*macro_grid_it, 0);
      // the fine scale reconstructions are only available for the barycenter of the macro grid entity (=> only
      // available for the canonical one point quadrature on this element)

      const int numDofs = elementOfRHS.numDofs(); // Dofs = Freiheitsgrade
      for (int i = 0; i < numDofs; ++i)
      {
        // --------------- the source contribution ( \int_{\Omega} f \Phi ) -------------------------------

        // the return values:
        RangeType f_x, phi_x;

        const int numMacroQuadraturePoints = macro_quadrature.nop();
        for (int quadraturePoint = 0; quadraturePoint < numMacroQuadraturePoints; ++quadraturePoint)
        {
          // local (barycentric) coordinates (with respect to entity)
          const typename Quadrature::CoordinateType& local_point = macro_quadrature.point(quadraturePoint);

          const DomainType global_point = macro_grid_geometry.global(local_point);

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

        const typename Quadrature::CoordinateType& local_macro_point = one_point_macro_quadrature.point(0 /*=quadraturePoint*/);

        // barycenter of macro grid entity
        const DomainType macro_entity_barycenter = macro_grid_geometry.global(local_macro_point);

        const double macro_entity_volume = one_point_macro_quadrature.weight(0 /*=quadraturePoint*/)
                                           * macro_grid_geometry.integrationElement(local_macro_point);

        // gradient of base function and gradient of old_u_H
        JacobianRangeType grad_Phi_x, grad_old_u_H_x;

        // evaluate the gradient of the current base function at the current quadrature point and save its value in
        // 'returnGradient':
        macro_grid_baseSet.jacobian(i, one_point_macro_quadrature[0], grad_Phi_x);
        // Bis jetzt nur der Gradient auf dem Referenzelement!!!!!!! Es muss noch transformiert werden, um den
        // Gradienten auf dem echten Element zu bekommen! Das passiert folgendermassen:

        const FieldMatrix< double, dimension, dimension >& inv
          = macro_grid_geometry.jacobianInverseTransposed(local_macro_point);
        // multiply with transpose of jacobian inverse
        grad_Phi_x[0] = FMatrixHelp::mult(inv, grad_Phi_x[0]);

        // get gradient of old u_H:
        old_u_H_loc.jacobian(one_point_macro_quadrature[0], grad_old_u_H_x);

        // Q_h(u_H^{(n-1}))(x_T,y):
        PeriodicDiscreteFunctionType corrector_old_u_H("Corrector of u_H^(n-1)", periodicDiscreteFunctionSpace);
        corrector_old_u_H.clear();

        PeriodicDiscreteFunctionType corrector_Phi_i("Corrector of Phi_i", periodicDiscreteFunctionSpace);
        corrector_Phi_i.clear();

        if (DSC_CONFIG_GET("AD_HOC_COMPUTATION", false)) {
          CellProblemSolverType cell_problem_solver(periodicDiscreteFunctionSpace, A);
          cell_problem_solver.template solvecellproblem< JacobianRangeType >
            (grad_old_u_H_x, macro_entity_barycenter, corrector_old_u_H);
          if (DSC_CONFIG_GET("TFR", false))
            cell_problem_solver.template solvecellproblem< JacobianRangeType >
              (grad_Phi_x, macro_entity_barycenter, corrector_Phi_i);

        } else {
          discrete_function_reader_discFunc.read(number_of_entity, corrector_old_u_H);
          if (DSC_CONFIG_GET("TFR", false)) {
            typename EntityType::EntityPointer entity_ptr(macro_grid_it);
            discrete_function_reader_baseSet.read(cp_num_manager.get_number_of_cell_problem(entity_ptr, i)
                                                  , corrector_Phi_i);
          }
        }

        RangeType fine_scale_contribution = 0.0;

        const auto micro_grid_end = periodicDiscreteFunctionSpace.end();
        for (auto micro_grid_it = periodicDiscreteFunctionSpace.begin();
             micro_grid_it != micro_grid_end;
             ++micro_grid_it)
        {
          const EntityType& micro_grid_entity = *micro_grid_it;
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
            if (DSC_CONFIG_GET("TFR", false))
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

            // if test function reconstruction
            if (DSC_CONFIG_GET("TFR", false)) {
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
