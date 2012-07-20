#ifndef DiscreteEllipticMSFEMOperator_HH
#define DiscreteEllipticMSFEMOperator_HH

#include <dune/common/fmatrix.hh>

#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/multiscale/tools/solver/MsFEM/msfem_localproblems/localproblemsolver.hh>

#include <dune/fem/operator/2order/lagrangematrixsetup.hh>

namespace Dune {
// Imp stands for Implementation
template< class DiscreteFunctionImp, class DiffusionImp, class LocalProblemNumberingManagerImp >
class DiscreteEllipticMsFEMOperator
  : public Operator< typename DiscreteFunctionImp::RangeFieldType, typename DiscreteFunctionImp::RangeFieldType,
                     DiscreteFunctionImp, DiscreteFunctionImp >
{
  typedef DiscreteEllipticMsFEMOperator< DiscreteFunctionImp, DiffusionImp, LocalProblemNumberingManagerImp > This;

public:
  typedef DiscreteFunctionImp             DiscreteFunction;
  typedef DiffusionImp                    DiffusionModel;
  typedef LocalProblemNumberingManagerImp LocalProblemNumberingManager;

  typedef typename DiscreteFunction::DiscreteFunctionSpaceType DiscreteFunctionSpace;

  typedef typename DiscreteFunctionSpace::GridPartType GridPart;
  typedef typename DiscreteFunctionSpace::GridType     GridType;

  typedef typename DiscreteFunctionSpace::RangeFieldType RangeFieldType;
  typedef typename DiscreteFunctionSpace::DomainType     DomainType;
  typedef typename DiscreteFunctionSpace::RangeType      RangeType;
  typedef typename DiscreteFunctionSpace::JacobianRangeType
  JacobianRangeType;

  #ifdef AD_HOC_COMPUTATION
  typedef MsFEMLocalProblemSolver< DiscreteFunction, DiffusionModel > MsFEMLocalProblemSolverType;
  #endif

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
  DiscreteEllipticMsFEMOperator(const DiscreteFunctionSpace& discreteFunctionSpace,
                                const DiscreteFunctionSpace& localDiscreteFunctionSpace,
                                DiffusionModel& diffusion_op,
                                LocalProblemNumberingManager& cp_num_manager)
    : discreteFunctionSpace_(discreteFunctionSpace)
      , localDiscreteFunctionSpace_(localDiscreteFunctionSpace)
      , diffusion_operator_(diffusion_op)
      , filename_(NULL)
  {}

  DiscreteEllipticMsFEMOperator(const DiscreteFunctionSpace& discreteFunctionSpace,
                                const DiscreteFunctionSpace& localDiscreteFunctionSpace,
                                DiffusionModel& diffusion_op,
                                LocalProblemNumberingManager& cp_num_manager,
                                std::string& filename)
    : discreteFunctionSpace_(discreteFunctionSpace)
      , localDiscreteFunctionSpace_(localDiscreteFunctionSpace)
      , diffusion_operator_(diffusion_op)
      , lp_num_manager_(cp_num_manager)
      , filename_(&filename)
  {}

private:
  DiscreteEllipticMsFEMOperator(const This&);

public:
  // dummy operator
  virtual void operator()(const DiscreteFunction& u, DiscreteFunction& w) const;

  template< class MatrixType >
  void assemble_matrix(MatrixType& global_matrix) const;

private:
  const DiscreteFunctionSpace& discreteFunctionSpace_;
  const DiscreteFunctionSpace& localDiscreteFunctionSpace_;
  DiffusionModel& diffusion_operator_;

  // name of data file, e.g. required if we want to use the saved solutions of the cell problems
  std::string* filename_;
  LocalProblemNumberingManager& lp_num_manager_;
};

// dummy implementation of "operator()"
// 'w' = effect of the discrete operator on 'u'
template< class DiscreteFunctionImp, class DiffusionImp, class LocalProblemNumberingManagerImp >
void DiscreteEllipticMsFEMOperator< DiscreteFunctionImp, DiffusionImp, LocalProblemNumberingManagerImp >::operator()(
  const DiscreteFunction& u,
  DiscreteFunction& w) const {
  DUNE_THROW(Dune::NotImplemented,"the ()-operator of the DiscreteEllipticMsFEMOperator class is not yet implemented and still a dummy.");
}

template< class DiscreteFunctionImp, class DiffusionImp, class LocalProblemNumberingManagerImp >
template< class MatrixType >
void DiscreteEllipticMsFEMOperator< DiscreteFunctionImp, DiffusionImp,
                                    LocalProblemNumberingManagerImp >::assemble_matrix(MatrixType& global_matrix) const
{
  // the local problem:
  // Let 'T' denote a coarse grid element and
  // let 'T_0' denote the reference element in 2D
  // ( the rectangular triangle with corners (0,0), (1,0) and (0,1) )
  // let 'Fx = Ax + b' denote the corresponding affine transformation with ' F(T) = T_0 '
  // Now, for fixed 'T' and for fixed coarse grid basis function Phi_i,
  // we solved for the solution (Q^eps(\Phi_i) ○ F) \in H^1(T_0) (and with zero boundary) of
  // \int_{T_0} (A^eps ○ F)(x) (A^{-1})^T ∇( Q^eps(\Phi_i) ○ F )(x) · ∇ \phi(x)
  // + \int_{T_0} (A^eps ○ F)(x) ∇ \Phi_H(x_T) · ∇ \phi(x) = 0

  // Now, we need to compute the stiffness matrix for the MsFEM, which is formally given by:
  // \sum_{T triangle} \int_T A^eps(x) ( ∇\Phi_i(x_T) + ∇Q^eps(\Phi_i)(x)) ·  ( ∇\Phi_j(x_T) +
  // ∇Q^eps(\Phi_j)(x))

  // After applying the transformation formula, we get (for the midpart):
  // \int_T A^eps(x) ( ∇\Phi_i(x_T) + ∇Q^eps(\Phi_i)(x)) ·  ( ∇\Phi_j(x_T) + ∇Q^eps(\Phi_j)(x)) =
  // |det(A)| \int_{T_0} (A^eps ○ F)(x) ( ∇\Phi_i(x_T) + (A^{-1})^T ∇( Q^eps(\Phi_i) ○ F )(x)) · (
  // ∇\Phi_j(x_T) + (A^{-1})^T ∇( Q^eps(\Phi_j) ○ F )(x))

  // Here, (A^{-1})^T denotes the transposed of the inverse of the matrix 'A' of F(x)=Ax+b
  // Note that ( Q^eps(\Phi_i) ○ F ) is precomputed for every i!

  // if test function reconstruction
  #ifdef PGF
  std::cout << "Assembling Petrov-Galerkin-MsFEM Matrix." << std::endl;
  #else
  std::cout << "Assembling MsFEM Matrix." << std::endl;
  #endif // ifdef PGF

  std::string cell_solution_location;

  // if we know the file, where we saved the solutions of the cell problems, use it
  if (filename_)
  {
    // place, where we saved the solutions of the cell problems
    cell_solution_location = "data/MsFEM/" + (*filename_) + "/local_problems/_localProblemSolutions_baseSet";
  } else {
    #ifndef AD_HOC_COMPUTATION
    DUNE_THROW(Dune::InvalidStateException,"ERROR! No 'filename_' in class 'DiscreteEllipticMsFEMOperator', but no AD_HOC_COMPUTATION initialized. Therefore the location of the saved cell problems is not available. Please define AD_HOC_COMPUTATION (ad hoc computation of the cell problems) or pass a corresponding 'filename_'-variable!");
    #endif // ifndef AD_HOC_COMPUTATION
  }

  bool reader_is_open = false;
  // reader for the cell problem data file:
  DiscreteFunctionReader discrete_function_reader( (cell_solution_location).c_str() );
  reader_is_open = discrete_function_reader.open();

  typedef typename MatrixType::LocalMatrixType LocalMatrix;

  global_matrix.reserve();
  global_matrix.clear();

  std::vector< typename BaseFunctionSet::JacobianRangeType > gradient_Phi( discreteFunctionSpace_.mapper().maxNumDofs() );

  const Iterator macro_grid_end = discreteFunctionSpace_.end();
  for (Iterator macro_grid_it = discreteFunctionSpace_.begin(); macro_grid_it != macro_grid_end; ++macro_grid_it)
  {
    // we compute (and sum up) the following:
    // |det(A)| \int_{T_0} (A^eps ○ F)(x) ( ∇\Phi_i(x_T) + (A^{-1})^T ∇( Q^eps(\Phi_i) ○ F )(x)) · (
    // ∇\Phi_j(x_T) + (A^{-1})^T ∇( Q^eps(\Phi_j) ○ F )(x))

    // the coarse grid element T:
    const Entity& macro_grid_entity = *macro_grid_it;
    const Geometry& macro_grid_geometry = macro_grid_entity.geometry();
    assert(macro_grid_entity.partitionType() == InteriorEntity);

    // transformation F : T_0 -> T
    // describe the mapping F(x) = Ax + b with F(T_0)=T for an entity T and the reference element T_0:
    // arguments: entity T, point in T_0, point in T.

    // Let (a_0,a_1,a_2) deonte the corners of the 2-simplex T, then the matrix A in the affine transformation
    // F(x) = Ax + a_0, F : T_0 -> T is given by
    // A_11 = a_1( 1 ) - a_0( 1 )     A_12 = a_2( 1 ) - a_0( 1 )
    // A_21 = a_1( 2 ) - a_0( 2 )     A_22 = a_2( 2 ) - a_0( 2 )

    // corners of the reference element:
    typename Quadrature::CoordinateType ref_corner_0, ref_corner_1, ref_corner_2;

    ref_corner_0[0] = 0.0;
    ref_corner_0[1] = 0.0;

    ref_corner_1[0] = 1.0;
    ref_corner_1[1] = 0.0;

    ref_corner_2[0] = 0.0;
    ref_corner_2[1] = 1.0;

    // corner of the global element:
    DomainType corner_0_of_T = macro_grid_geometry.global(ref_corner_0);
    DomainType corner_1_of_T = macro_grid_geometry.global(ref_corner_1);
    DomainType corner_2_of_T = macro_grid_geometry.global(ref_corner_2);

    // value of the matrix A (in F(x) = Ax + a_0)
    double val_A[dimension][dimension];
    val_A[0][0] = corner_1_of_T[0] - corner_0_of_T[0];
    val_A[0][1] = corner_2_of_T[0] - corner_0_of_T[0];
    val_A[1][0] = corner_1_of_T[1] - corner_0_of_T[1];
    val_A[1][1] = corner_2_of_T[1] - corner_0_of_T[1];

    // define 'c := (a_1(1) - a_0(1))·(a_2(2) - a_0(2)) - (a_1(2) - a_0(2))·(a_2(1) - a_0(1))
    double c = 1.0 / ( (val_A[0][0] * val_A[1][1]) - (val_A[0][1] * val_A[1][0]) );
    // then the inverse A^{-1} is given by:
    // A^{-1}_11 = (1/c) (a_2(2) - a_0(2))     A^{-1}_12 = (1/c) (a_0(1) - a_2(1))
    // A^{-1}_21 = (1/c) (a_0(2) - a_1(2))     A^{-1}_22 = (1/c) (a_1(1) - a_0(1))

    // |det(A)|:
    double abs_det_A = fabs(1.0 / c);

    // (A^{-1})^T:
    double val_A_inverse_transposed[dimension][dimension];
    val_A_inverse_transposed[0][0] = c * val_A[1][1];
    val_A_inverse_transposed[1][0] = c * (-1.0) * val_A[0][1];
    val_A_inverse_transposed[0][1] = c * (-1.0) * val_A[1][0];
    val_A_inverse_transposed[1][1] = c * val_A[0][0];

    LocalMatrix local_matrix = global_matrix.localMatrix(macro_grid_entity, macro_grid_entity);

    const BaseFunctionSet& macro_grid_baseSet = local_matrix.domainBaseFunctionSet();
    const unsigned int numMacroBaseFunctions = macro_grid_baseSet.numBaseFunctions();

    // 1 point quadrature!! That is how we compute and save the cell problems.
    // If you want to use a higher order quadrature, you also need to change the computation of the cell problems!
    Quadrature one_point_quadrature(macro_grid_entity, 0);

    // the barycenter of the macro_grid_entity
    const typename Quadrature::CoordinateType& local_macro_point = one_point_quadrature.point(0 /*=quadraturePoint*/);
    DomainType macro_entity_barycenter = macro_grid_geometry.global(local_macro_point);

    const double macro_entity_volume = one_point_quadrature.weight(0 /*=quadraturePoint*/)
                                       * macro_grid_geometry.integrationElement(local_macro_point);

    // transposed of the the inverse jacobian
    const FieldMatrix< double, dimension, dimension >& inverse_jac
      = macro_grid_geometry.jacobianInverseTransposed(local_macro_point);

    int cell_problem_id[numMacroBaseFunctions];

    DiscreteFunction* corrector_Phi[discreteFunctionSpace_.mapper().maxNumDofs()];

    for (unsigned int i = 0; i < numMacroBaseFunctions; ++i)
    {
      // get number of cell problem from entity and number of base function
      cell_problem_id[i] = lp_num_manager_.get_number_of_local_problem(macro_grid_it, i);

      // jacobian of the base functions, with respect to the reference element
      typename BaseFunctionSet::JacobianRangeType gradient_Phi_ref_element;
      macro_grid_baseSet.jacobian(i, one_point_quadrature[0], gradient_Phi_ref_element);

      // multiply it with transpose of jacobian inverse to obtain the jacobian with respect to the real entity
      inverse_jac.mv(gradient_Phi_ref_element[0], gradient_Phi[i][0]);

      // ( Q^eps(\Phi_i) ○ F ):
      corrector_Phi[i] = new DiscreteFunction("Corrector Function of Phi", localDiscreteFunctionSpace_);
      corrector_Phi[i]->clear();
      #ifdef AD_HOC_COMPUTATION
      MsFEMLocalProblemSolverType cell_problem_solver(localDiscreteFunctionSpace_, diffusion_operator_);
      cell_problem_solver.template solvelocalproblem< JacobianRangeType >
        ( gradient_Phi[i], macro_entity_barycenter, *(corrector_Phi[i]) );
      #else // ifdef AD_HOC_COMPUTATION
      if (reader_is_open)
      {
        discrete_function_reader.read( cell_problem_id[i], *(corrector_Phi[i]) );
      }
      #endif // ifdef AD_HOC_COMPUTATION
    }

    for (unsigned int i = 0; i < numMacroBaseFunctions; ++i)
    {
      for (unsigned int j = 0; j < numMacroBaseFunctions; ++j)
      {
        RangeType local_integral = 0.0;
        // iterator for the micro grid ( grid for the reference element T_0 )
        const Iterator micro_grid_end = localDiscreteFunctionSpace_.end();
        for (Iterator micro_grid_it = localDiscreteFunctionSpace_.begin();
             micro_grid_it != micro_grid_end;
             ++micro_grid_it)
        {
          // remember:
          // |det(A)| \int_{T_0} (A^eps ○ F)(x) ( ∇\Phi_i(x_T) + (A^{-1})^T ∇( Q^eps(\Phi_i) ○ F )(x)) · (
          // ∇\Phi_j(x_T) + (A^{-1})^T ∇( Q^eps(\Phi_j) ○ F )(x))
          const Entity& micro_grid_entity = *micro_grid_it;
          const Geometry& micro_grid_geometry = micro_grid_entity.geometry();
          assert(micro_grid_entity.partitionType() == InteriorEntity);

          // ( Q^eps(\Phi_i) ○ F ):
          typename DiscreteFunction::LocalFunctionType localized_corrector_i = corrector_Phi[i]->localFunction(
            micro_grid_entity);
          // ( Q^eps(\Phi_j) ○ F ):
          typename DiscreteFunction::LocalFunctionType localized_corrector_j = corrector_Phi[j]->localFunction(
            micro_grid_entity);

          // higher order quadrature, since A^{\epsilon} is highly variable
          Quadrature micro_grid_quadrature(micro_grid_entity, 2 * localDiscreteFunctionSpace_.order() + 2);
          const size_t numQuadraturePoints = micro_grid_quadrature.nop();

          for (size_t microQuadraturePoint = 0; microQuadraturePoint < numQuadraturePoints; ++microQuadraturePoint)
          {
            // local (barycentric) coordinates (with respect to entity)
            const typename Quadrature::CoordinateType& local_micro_point = micro_grid_quadrature.point(
              microQuadraturePoint);

            DomainType global_point_in_T_0 = micro_grid_geometry.global(local_micro_point);

            const double weight_micro_quadrature = micro_grid_quadrature.weight(microQuadraturePoint)
                                                   * micro_grid_geometry.integrationElement(local_micro_point);

            // ∇( Q^eps(\Phi_i) ○ F )   and   ∇( Q^eps(\Phi_j) ○ F )
            JacobianRangeType grad_corrector_i, grad_corrector_j;
            localized_corrector_i.jacobian(micro_grid_quadrature[microQuadraturePoint], grad_corrector_i);
            localized_corrector_j.jacobian(micro_grid_quadrature[microQuadraturePoint], grad_corrector_j);

            // global point in the reference element T_0
            DomainType global_point = micro_grid_geometry.global(local_micro_point);

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

            // direction of the diffusion: ( ∇\Phi_i(x_T) + (A^{-1})^T ∇( Q^eps(\Phi_i) ○ F )(x))
            JacobianRangeType direction_of_diffusion(0.0);
            for (int k = 0; k < dimension; ++k)
            {
              for (int l = 0; l < dimension; ++l)
              {
                direction_of_diffusion[0][k] += val_A_inverse_transposed[k][l] * grad_corrector_i[0][l];
              }

              direction_of_diffusion[0][k] += gradient_Phi[i][0][k];
            }

            JacobianRangeType diffusion_in_gradient_Phi_reconstructed;
            diffusion_operator_.diffusiveFlux(global_point_transformed,
                                              direction_of_diffusion, diffusion_in_gradient_Phi_reconstructed);

            // if test function reconstruction
            #ifndef PGF

            // ( ∇\Phi_j(x_T) + (A^{-1})^T ∇( Q^eps(\Phi_j) ○ F )(x)):
            JacobianRangeType grad_reconstruction_Phi_j(0.0);
            for (int k = 0; k < dimension; ++k)
            {
              for (int l = 0; l < dimension; ++l)
              {
                grad_reconstruction_Phi_j[0][k] += val_A_inverse_transposed[k][l] * grad_corrector_j[0][l];
              }
              grad_reconstruction_Phi_j[0][k] += gradient_Phi[j][0][k];
            }

            local_integral += weight_micro_quadrature
                              * (diffusion_in_gradient_Phi_reconstructed[0] * grad_reconstruction_Phi_j[0]);
            #else // ifndef PGF
            local_integral += weight_micro_quadrature
                              * (diffusion_in_gradient_Phi_reconstructed[0] * gradient_Phi[j][0]);
            #endif // ifndef PGF
          }
        }
        // add |det(A)|*\int_{T_0} ...
        local_matrix.add(j, i, abs_det_A * local_integral);
      }
    }

    // delete?
    // delete[] corrector_Phi;
  }

  // discrete_function_reader.close();

  // boundary treatment
  const GridPart& gridPart = discreteFunctionSpace_.gridPart();
  for (Iterator it = discreteFunctionSpace_.begin(); it != macro_grid_end; ++it)
  {
    const Entity& entity = *it;
    if ( !entity.hasBoundaryIntersections() )
      continue;

    LocalMatrix local_matrix = global_matrix.localMatrix(entity, entity);

    const LagrangePointSet& lagrangePointSet = discreteFunctionSpace_.lagrangePointSet(entity);

    const IntersectionIterator iend = gridPart.iend(entity);
    for (IntersectionIterator iit = gridPart.ibegin(entity); iit != iend; ++iit)
    {
      const Intersection& intersection = *iit;
      if ( !intersection.boundary() )
        continue;

      const int face = intersection.indexInInside();
      const FaceDofIterator fdend = lagrangePointSet.template endSubEntity< 1 >(face);
      for (FaceDofIterator fdit = lagrangePointSet.template beginSubEntity< 1 >(face); fdit != fdend; ++fdit)
        local_matrix.unitRow(*fdit);
    }
  }
} // assemble_matrix

// ! ------------------------------------------------------------------------------------------------
// ! ------------------------------------------------------------------------------------------------
}

#endif // #ifndef DiscreteElliptic_HH
