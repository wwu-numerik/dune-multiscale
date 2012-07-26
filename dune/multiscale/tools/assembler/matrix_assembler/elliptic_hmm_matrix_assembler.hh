#ifndef DiscreteEllipticHMMOperator_HH
#define DiscreteEllipticHMMOperator_HH

#include <boost/noncopyable.hpp>

#include <dune/common/fmatrix.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/operator/2order/lagrangematrixsetup.hh>
#include <dune/multiscale/tools/solver/HMM/cell_problem_solving/discreteoperator.hh>


namespace Dune {
// Imp stands for Implementation
template< class DiscreteFunctionImp, class PeriodicDiscreteFunctionImp, class DiffusionImp,
          class CellProblemNumberingManagerImp >
class DiscreteEllipticHMMOperator
  : public Operator< typename DiscreteFunctionImp::RangeFieldType, typename DiscreteFunctionImp::RangeFieldType,
                     DiscreteFunctionImp, DiscreteFunctionImp >
    , boost::noncopyable
{
  typedef DiscreteEllipticHMMOperator< DiscreteFunctionImp, PeriodicDiscreteFunctionImp, DiffusionImp,
                                       CellProblemNumberingManagerImp > This;

public:
  typedef DiscreteFunctionImp            DiscreteFunction;
  typedef PeriodicDiscreteFunctionImp    PeriodicDiscreteFunction;
  typedef DiffusionImp                   DiffusionModel;
  typedef CellProblemNumberingManagerImp CellProblemNumberingManager;

  typedef typename DiscreteFunction::DiscreteFunctionSpaceType         DiscreteFunctionSpace;
  typedef typename PeriodicDiscreteFunction::DiscreteFunctionSpaceType PeriodicDiscreteFunctionSpace;

  typedef typename DiscreteFunctionSpace::GridPartType GridPart;
  typedef typename DiscreteFunctionSpace::GridType     GridType;

  typedef typename PeriodicDiscreteFunctionSpace::GridPartType PeriodicGridPart;
  typedef typename PeriodicDiscreteFunctionSpace::GridType     PeriodicGridType;

  typedef typename DiscreteFunctionSpace::RangeFieldType RangeFieldType;
  typedef typename DiscreteFunctionSpace::DomainType     DomainType;
  typedef typename DiscreteFunctionSpace::RangeType      RangeType;
  typedef typename DiscreteFunctionSpace::JacobianRangeType
  JacobianRangeType;

  #ifdef AD_HOC_COMPUTATION
  typedef CellProblemSolver< PeriodicDiscreteFunction, DiffusionModel > CellProblemSolverType;
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
  DiscreteEllipticHMMOperator(const DiscreteFunctionSpace& discreteFunctionSpace,
                              const PeriodicDiscreteFunctionSpace& periodicDiscreteFunctionSpace,
                              DiffusionModel& diffusion_op,
                              CellProblemNumberingManager& cp_num_manager)
    : discreteFunctionSpace_(discreteFunctionSpace)
      , periodicDiscreteFunctionSpace_(periodicDiscreteFunctionSpace)
      , diffusion_operator_(diffusion_op)
      , cp_num_manager_(cp_num_manager)
      , filename_(NULL)
  {}

  DiscreteEllipticHMMOperator(const DiscreteFunctionSpace& discreteFunctionSpace,
                              const PeriodicDiscreteFunctionSpace& periodicDiscreteFunctionSpace,
                              DiffusionModel& diffusion_op,
                              CellProblemNumberingManager& cp_num_manager,
                              const std::string& filename)
    : discreteFunctionSpace_(discreteFunctionSpace)
      , periodicDiscreteFunctionSpace_(periodicDiscreteFunctionSpace)
      , diffusion_operator_(diffusion_op)
      , cp_num_manager_(cp_num_manager)
      , filename_(&filename)
  {}

public:
  // dummy operator
  virtual void operator()(const DiscreteFunction& u, DiscreteFunction& w) const;

  template< class MatrixType >
  void assemble_matrix(MatrixType& global_matrix) const;

  template< class MatrixType >
  void assemble_jacobian_matrix(DiscreteFunction& old_macro_function, MatrixType& global_matrix) const;

private:
  const DiscreteFunctionSpace& discreteFunctionSpace_;
  const PeriodicDiscreteFunctionSpace& periodicDiscreteFunctionSpace_;
  DiffusionModel& diffusion_operator_;
  CellProblemNumberingManager& cp_num_manager_;

  // name of data file, e.g. required if we want to use the saved solutions of the cell problems
  const std::string* filename_;
};

// dummy implementation of "operator()"
// 'w' = effect of the discrete operator on 'u'
template< class DiscreteFunctionImp, class PeriodicDiscreteFunctionImp, class DiffusionImp,
          class CellProblemNumberingManagerImp >
void DiscreteEllipticHMMOperator< DiscreteFunctionImp, PeriodicDiscreteFunctionImp, DiffusionImp,
                                  CellProblemNumberingManagerImp >::operator()(const DiscreteFunction& /*u*/,
                                                                               DiscreteFunction& /*w*/) const {
  DUNE_THROW(Dune::NotImplemented,"the ()-operator of the DiscreteEllipticHMMOperator class is not yet implemented and still a dummy.");
}

template< class DiscreteFunctionImp, class PeriodicDiscreteFunctionImp, class DiffusionImp,
          class CellProblemNumberingManagerImp >
template< class MatrixType >
void DiscreteEllipticHMMOperator< DiscreteFunctionImp, PeriodicDiscreteFunctionImp, DiffusionImp,
                                  CellProblemNumberingManagerImp >::assemble_matrix(MatrixType& global_matrix) const {
  // if test function reconstruction
  #ifdef TFR
  DSC_LOG_INFO << "Assembling TFR-HMM Matrix." << std::endl;
  #else
  DSC_LOG_INFO << "Assembling HMM Matrix." << std::endl;
  #endif // ifdef TFR

  std::string cell_solution_location;

  // if we know the file, where we saved the solutions of the cell problems, use it
  if (filename_)
  {
    // place, where we saved the solutions of the cell problems
    cell_solution_location = "data/" + (*filename_) + "/cell_problems/_cellSolutions_baseSet";
  } else {
    #ifndef AD_HOC_COMPUTATION
    DUNE_THROW(Dune::InvalidStateException,"ERROR! No 'filename_' in class 'DiscreteEllipticHMMOperator', but no AD_HOC_COMPUTATION initialized. Therefore the location of the saved cell problems is not available. Please define AD_HOC_COMPUTATION (ad hoc computation of the cell problems) or pass a corresponding 'filename_'-variable!");
    #endif // ifndef AD_HOC_COMPUTATION
  }

  Problem::ModelProblemData model_info;
  const double delta = model_info.getDelta();
  const double epsilon_estimated = model_info.getEpsilonEstimated();


  // reader for the cell problem data file:
  DiscreteFunctionReader discrete_function_reader( (cell_solution_location).c_str() );
  const bool reader_is_open = discrete_function_reader.open();

  typedef typename MatrixType::LocalMatrixType LocalMatrix;

  global_matrix.reserve();
  global_matrix.clear();

  std::vector< typename BaseFunctionSet::JacobianRangeType > gradient_Phi( discreteFunctionSpace_.mapper().maxNumDofs() );

  const Iterator macro_grid_end = discreteFunctionSpace_.end();
  for (Iterator macro_grid_it = discreteFunctionSpace_.begin(); macro_grid_it != macro_grid_end; ++macro_grid_it)
  {
    const Entity& macro_grid_entity = *macro_grid_it;
    const Geometry& macro_grid_geometry = macro_grid_entity.geometry();
    assert(macro_grid_entity.partitionType() == InteriorEntity);

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

    PeriodicDiscreteFunction* corrector_Phi[discreteFunctionSpace_.mapper().maxNumDofs()];

    for (unsigned int i = 0; i < numMacroBaseFunctions; ++i)
    {
      // get number of cell problem from entity and number of base function
      typename Entity::EntityPointer macro_entity_pointer(*macro_grid_it);
      cell_problem_id[i] = cp_num_manager_.get_number_of_cell_problem(macro_entity_pointer, i);

      // jacobian of the base functions, with respect to the reference element
      typename BaseFunctionSet::JacobianRangeType gradient_Phi_ref_element;
      macro_grid_baseSet.jacobian(i, one_point_quadrature[0], gradient_Phi_ref_element);

      // multiply it with transpose of jacobian inverse to obtain the jacobian with respect to the real entity
      inverse_jac.mv(gradient_Phi_ref_element[0], gradient_Phi[i][0]);

      corrector_Phi[i] = new PeriodicDiscreteFunction("Corrector Function of Phi", periodicDiscreteFunctionSpace_);
      corrector_Phi[i]->clear();
      #ifdef AD_HOC_COMPUTATION
      CellProblemSolverType cell_problem_solver(periodicDiscreteFunctionSpace_, diffusion_operator_);
      cell_problem_solver.template solvecellproblem< JacobianRangeType >
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
        RangeType fine_scale_average = 0.0;

        // nur checken ob der momentane Quadraturpunkt in der Zelle delta/epsilon_estimated*Y ist (0 bzw. 1 =>
        // Abschneidefunktion!)

        const Iterator micro_grid_end = periodicDiscreteFunctionSpace_.end();
        for (Iterator micro_grid_it = periodicDiscreteFunctionSpace_.begin();
             micro_grid_it != micro_grid_end;
             ++micro_grid_it)
        {
          const Entity& micro_grid_entity = *micro_grid_it;
          const Geometry& micro_grid_geometry = micro_grid_entity.geometry();
          assert(micro_grid_entity.partitionType() == InteriorEntity);

          typename PeriodicDiscreteFunction::LocalFunctionType localized_corrector_i = corrector_Phi[i]->localFunction(
            micro_grid_entity);
          typename PeriodicDiscreteFunction::LocalFunctionType localized_corrector_j = corrector_Phi[j]->localFunction(
            micro_grid_entity);

          // higher order quadrature, since A^{\epsilon} is highly variable
          Quadrature micro_grid_quadrature(micro_grid_entity, 2 * periodicDiscreteFunctionSpace_.order() + 2);
          const size_t numQuadraturePoints = micro_grid_quadrature.nop();

          for (size_t microQuadraturePoint = 0; microQuadraturePoint < numQuadraturePoints; ++microQuadraturePoint)
          {
            // local (barycentric) coordinates (with respect to entity)
            const typename Quadrature::CoordinateType& local_micro_point = micro_grid_quadrature.point(
              microQuadraturePoint);

            DomainType global_point_in_Y = micro_grid_geometry.global(local_micro_point);

            const double weight_micro_quadrature = micro_grid_quadrature.weight(microQuadraturePoint)
                                                   * micro_grid_geometry.integrationElement(local_micro_point);

            JacobianRangeType grad_corrector_i, grad_corrector_j;
            localized_corrector_i.jacobian(micro_grid_quadrature[microQuadraturePoint], grad_corrector_i);
            localized_corrector_j.jacobian(micro_grid_quadrature[microQuadraturePoint], grad_corrector_j);

            // x_T + (delta * y)
            DomainType current_point_in_macro_grid;
            for (int k = 0; k < dimension; ++k)
              current_point_in_macro_grid[k] = macro_entity_barycenter[k] + (delta * global_point_in_Y[k]);

            JacobianRangeType direction_of_diffusion;
            for (int k = 0; k < dimension; ++k)
              direction_of_diffusion[0][k] = gradient_Phi[i][0][k] + grad_corrector_i[0][k];

            JacobianRangeType diffusion_in_gradient_Phi_reconstructed;
            diffusion_operator_.diffusiveFlux(current_point_in_macro_grid,
                                              direction_of_diffusion, diffusion_in_gradient_Phi_reconstructed);

            double cutting_function = 1.0;
            for (int k = 0; k < dimension; ++k)
            {
              // is the current quadrature point in the relevant cell?
              if ( fabs(global_point_in_Y[k]) > ( 0.5 * (epsilon_estimated / delta) ) )
              {
                cutting_function *= 0.0;
              }
            }

            // if test function reconstruction
            #ifdef TFR
            JacobianRangeType grad_reconstruction_Phi_j;
            for (int k = 0; k < dimension; ++k)
              grad_reconstruction_Phi_j[0][k] = gradient_Phi[j][0][k] + grad_corrector_j[0][k];

            fine_scale_average += cutting_function * weight_micro_quadrature
                                  * (diffusion_in_gradient_Phi_reconstructed[0] * grad_reconstruction_Phi_j[0]);
            #else // ifdef TFR
            fine_scale_average += cutting_function * weight_micro_quadrature
                                  * (diffusion_in_gradient_Phi_reconstructed[0] * gradient_Phi[j][0]);
            #endif // ifdef TFR
          }
        }

        // add |T| * (delta/epsilon)^N \int_Y ...
        local_matrix.add(j, i,
                         pow(delta / epsilon_estimated, dimension) * macro_entity_volume * fine_scale_average);
      }
    }

    // delete?
    // delete[] corrector_Phi;
  }

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

// assemble stiffness matrix for HMM with Newton Method
template< class DiscreteFunctionImp, class PeriodicDiscreteFunctionImp, class DiffusionImp,
          class CellProblemNumberingManagerImp >
template< class MatrixType >
void DiscreteEllipticHMMOperator< DiscreteFunctionImp, PeriodicDiscreteFunctionImp, DiffusionImp,
                                  CellProblemNumberingManagerImp >::assemble_jacobian_matrix(DiscreteFunction& old_u_H /*u_H^(n-1)*/,
                                                                                             MatrixType& global_matrix)
const {
  // if test function reconstruction
  #ifdef TFR
  DSC_LOG_INFO << "Assembling TFR-HMM Matrix for Newton Iteration." << std::endl;
  #else
  DSC_LOG_INFO << "Assembling HMM Matrix for Newton Iteration." << std::endl;
  #endif // ifdef TFR

  std::string cell_solution_location_baseSet;
  std::string cell_solution_location_discFunc;
  std::string jac_cor_cell_solution_location_baseSet_discFunc;

  // if we know the file, where we saved the solutions of the cell problems, use it
  if (filename_)
  {
    // place, where we saved the solutions of the cell problems
    cell_solution_location_baseSet = "data/HMM/" + (*filename_) + "/cell_problems/_cellSolutions_baseSet";
    cell_solution_location_discFunc = "data/HMM/" + (*filename_) + "/cell_problems/_cellSolutions_discFunc";
    jac_cor_cell_solution_location_baseSet_discFunc = "data/HMM/" + (*filename_)
                                                      + "/cell_problems/_JacCorCellSolutions_baseSet_discFunc";
  } else {
    #ifndef AD_HOC_COMPUTATION
    DUNE_THROW(Dune::InvalidStateException,"ERROR! No 'filename_' in class 'DiscreteEllipticHMMOperator', but no AD_HOC_COMPUTATION initialized. Therefore the location of the saved cell problems is not available. Please define AD_HOC_COMPUTATION (ad hoc computation of the cell problems) or pass a corresponding 'filename_'-variable!");
    #endif // ifndef AD_HOC_COMPUTATION
  }

  Problem::ModelProblemData model_info;
  const double delta = model_info.getDelta();
  const double epsilon_estimated = model_info.getEpsilonEstimated();

  // reader for the cell problem data file:
  DiscreteFunctionReader discrete_function_reader_baseSet( (cell_solution_location_baseSet).c_str() );
  discrete_function_reader_baseSet.open();

  // reader for the cell problem data file:
  DiscreteFunctionReader discrete_function_reader_discFunc( (cell_solution_location_discFunc).c_str() );
  discrete_function_reader_discFunc.open();

  // reader for the cell problem data file:
  DiscreteFunctionReader discrete_function_reader_jac_cor( (jac_cor_cell_solution_location_baseSet_discFunc).c_str() );
//  const bool reader_is_open = discrete_function_reader_jac_cor.open();

  typedef typename MatrixType::LocalMatrixType LocalMatrix;

  typedef typename DiscreteFunction::LocalFunctionType
  LocalFunction;

  global_matrix.reserve();
  global_matrix.clear();

  std::vector< typename BaseFunctionSet::JacobianRangeType > gradient_Phi( discreteFunctionSpace_.mapper().maxNumDofs() );
  std::vector< typename BaseFunctionSet::JacobianRangeType > gradient_Phi_new( discreteFunctionSpace_.mapper().maxNumDofs() );

  int number_of_macro_entity = 0;

  const Iterator macro_grid_end = discreteFunctionSpace_.end();
  for (Iterator macro_grid_it = discreteFunctionSpace_.begin(); macro_grid_it != macro_grid_end; ++macro_grid_it)
  {
    const Entity& macro_grid_entity = *macro_grid_it;
    const Geometry& macro_grid_geometry = macro_grid_entity.geometry();
    assert(macro_grid_entity.partitionType() == InteriorEntity);

    LocalMatrix local_matrix = global_matrix.localMatrix(macro_grid_entity, macro_grid_entity);
    LocalFunction local_old_u_H = old_u_H.localFunction(macro_grid_entity);

    const BaseFunctionSet& macro_grid_baseSet = local_matrix.domainBaseFunctionSet();
    const unsigned int numMacroBaseFunctions = macro_grid_baseSet.size();

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

    std::vector<int> cell_problem_id(numMacroBaseFunctions, -1);

    // \nabla_x u_H^{(n-1})(x_T)
    typename BaseFunctionSet::JacobianRangeType grad_old_u_H;
    local_old_u_H.jacobian(one_point_quadrature[0], grad_old_u_H);
    // here: no multiplication with jacobian inverse transposed required!

    // Q_h(u_H^{(n-1}))(x_T,y):
    PeriodicDiscreteFunction corrector_old_u_H("Corrector of u_H^(n-1)", periodicDiscreteFunctionSpace_);
    corrector_old_u_H.clear();

    #ifdef AD_HOC_COMPUTATION
    CellProblemSolverType cell_problem_solver(periodicDiscreteFunctionSpace_, diffusion_operator_);
    cell_problem_solver.template solvecellproblem< JacobianRangeType >
      (grad_old_u_H, macro_entity_barycenter, corrector_old_u_H);
    #else // ifdef AD_HOC_COMPUTATION
    discrete_function_reader_discFunc.read(number_of_macro_entity, corrector_old_u_H);

    #endif // ifdef AD_HOC_COMPUTATION

    #ifdef TFR
    PeriodicDiscreteFunction* corrector_Phi[discreteFunctionSpace_.mapper().maxNumDofs()];
    #endif

    // gradients of macrocopic base functions:
    for (unsigned int i = 0; i < numMacroBaseFunctions; ++i)
    {
      // get number of cell problem from entity and number of base function
      typename Entity::EntityPointer macro_entity_pointer(*macro_grid_it);
      cell_problem_id[i] = cp_num_manager_.get_number_of_cell_problem(macro_entity_pointer, i);

      // jacobian of the base functions, with respect to the reference element
      typename BaseFunctionSet::JacobianRangeType gradient_Phi_ref_element;
      macro_grid_baseSet.jacobian(i, one_point_quadrature[0], gradient_Phi_ref_element);

      // multiply it with transpose of jacobian inverse to obtain the jacobian with respect to the real entity
      inverse_jac.mv(gradient_Phi_ref_element[0], gradient_Phi[i][0]);

      #ifdef TFR
      corrector_Phi[i] = new PeriodicDiscreteFunction("Corrector Function of Phi_j", periodicDiscreteFunctionSpace_);
      corrector_Phi[i]->clear();
      #ifdef AD_HOC_COMPUTATION
      cell_problem_solver.template solvecellproblem< JacobianRangeType >
        ( gradient_Phi[i], macro_entity_barycenter, *(corrector_Phi[i]) );
      #else // ifdef AD_HOC_COMPUTATION
      discrete_function_reader_baseSet.read( cell_problem_id[i], *(corrector_Phi[i]) );
      #endif // ifdef AD_HOC_COMPUTATION
      #endif // ifdef TFR
    }
    // the multiplication with jacobian inverse is delegated
    macro_grid_baseSet.jacobianAll(one_point_quadrature[0], inverse_jac, gradient_Phi_new);
    assert( gradient_Phi == gradient_Phi_new );

    for (unsigned int i = 0; i < numMacroBaseFunctions; ++i)
    {
      // D_Q(\Phi_i,u_H^{n-1})
      // the jacobian of the corrector operator applied to u_H^{(n-1)} in direction of gradient \Phi_i
      PeriodicDiscreteFunction jacobian_corrector_old_u_H_Phi_i("Jacobian Corrector Function of u_H^(n-1) and Phi_i",
                                                                periodicDiscreteFunctionSpace_);
      jacobian_corrector_old_u_H_Phi_i.clear();

      #ifdef AD_HOC_COMPUTATION
      cell_problem_solver.template solve_jacobiancorrector_cellproblem< JacobianRangeType >
        (gradient_Phi[i],
        grad_old_u_H,
        corrector_old_u_H,
        macro_entity_barycenter,
        jacobian_corrector_old_u_H_Phi_i);
      #else // ifdef AD_HOC_COMPUTATION
      discrete_function_reader_jac_cor.read(cell_problem_id[i], jacobian_corrector_old_u_H_Phi_i);
      #endif // ifdef AD_HOC_COMPUTATION

      for (unsigned int j = 0; j < numMacroBaseFunctions; ++j)
      {
        RangeType fine_scale_average = 0.0;

        const Iterator micro_grid_end = periodicDiscreteFunctionSpace_.end();
        for (Iterator micro_grid_it = periodicDiscreteFunctionSpace_.begin();
             micro_grid_it != micro_grid_end;
             ++micro_grid_it)
        {
          const Entity& micro_grid_entity = *micro_grid_it;
          const Geometry& micro_grid_geometry = micro_grid_entity.geometry();
          assert(micro_grid_entity.partitionType() == InteriorEntity);

          typename PeriodicDiscreteFunction::LocalFunctionType loc_corrector_old_u_H = corrector_old_u_H.localFunction(
            micro_grid_entity);

          typename PeriodicDiscreteFunction::LocalFunctionType loc_D_Q_old_u_H_Phi_i
            = jacobian_corrector_old_u_H_Phi_i.localFunction(micro_grid_entity);

          #ifdef TFR
          typename PeriodicDiscreteFunction::LocalFunctionType loc_corrector_Phi_j = corrector_Phi[j]->localFunction(
            micro_grid_entity);
          #endif // ifdef TFR

          // higher order quadrature, since A^{\epsilon} is highly variable
          Quadrature micro_grid_quadrature(micro_grid_entity, 2 * periodicDiscreteFunctionSpace_.order() + 2);
          const size_t numQuadraturePoints = micro_grid_quadrature.nop();

          for (size_t microQuadraturePoint = 0; microQuadraturePoint < numQuadraturePoints; ++microQuadraturePoint)
          {
            // local (barycentric) coordinates (with respect to entity)
            const typename Quadrature::CoordinateType& local_micro_point = micro_grid_quadrature.point(
              microQuadraturePoint);

            DomainType global_point_in_Y = micro_grid_geometry.global(local_micro_point);

            const double weight_micro_quadrature = micro_grid_quadrature.weight(microQuadraturePoint)
                                                   * micro_grid_geometry.integrationElement(local_micro_point);

            JacobianRangeType grad_corrector_old_u_H, grad_D_Q_old_u_H_Phi_i;
            loc_corrector_old_u_H.jacobian(micro_grid_quadrature[microQuadraturePoint], grad_corrector_old_u_H);
            loc_D_Q_old_u_H_Phi_i.jacobian(micro_grid_quadrature[microQuadraturePoint], grad_D_Q_old_u_H_Phi_i);

            #ifdef TFR
            JacobianRangeType grad_corrector_Phi_j;
            loc_corrector_Phi_j.jacobian(micro_grid_quadrature[microQuadraturePoint], grad_corrector_Phi_j);
            #endif // ifdef TFR

            // x_T + (delta * y)
            DomainType current_point_in_macro_grid;
            for (int k = 0; k < dimension; ++k)
              current_point_in_macro_grid[k] = macro_entity_barycenter[k] + (delta * global_point_in_Y[k]);

            // evaluate jacobian matrix of diffusion operator in 'position_vector' in direction 'direction_vector':

            JacobianRangeType position_vector;
            for (int k = 0; k < dimension; ++k)
              position_vector[0][k] = grad_old_u_H[0][k] + grad_corrector_old_u_H[0][k];

            JacobianRangeType direction_vector;
            for (int k = 0; k < dimension; ++k)
              direction_vector[0][k] = gradient_Phi[i][0][k] + grad_D_Q_old_u_H_Phi_i[0][k];

            typename LocalFunction::JacobianRangeType jac_diffusion_flux;
            diffusion_operator_.jacobianDiffusiveFlux(current_point_in_macro_grid,
                                                      position_vector,
                                                      direction_vector,
                                                      jac_diffusion_flux);

            double cutting_function = 1.0;
            for (int k = 0; k < dimension; ++k)
            {
              // is the current quadrature point in the relevant cell?
              if ( fabs(global_point_in_Y[k]) > ( 0.5 * (epsilon_estimated / delta) ) )
              {
                cutting_function *= 0.0;
              }
            }

            // if test function reconstruction
            #ifdef TFR
            JacobianRangeType grad_reconstruction_Phi_j;
            for (int k = 0; k < dimension; ++k)
              grad_reconstruction_Phi_j[0][k] = gradient_Phi[j][0][k] + grad_corrector_Phi_j[0][k];

            fine_scale_average += cutting_function * weight_micro_quadrature
                                  * (jac_diffusion_flux[0] * grad_reconstruction_Phi_j[0]);
            #else // ifdef TFR
            fine_scale_average += cutting_function * weight_micro_quadrature
                                  * (jac_diffusion_flux[0] * gradient_Phi[j][0]);
            #endif // ifdef TFR
          }
        }

        // add |T| * (delta/epsilon)^N \int_Y ...
        local_matrix.add(j, i,
                         pow(delta / epsilon_estimated, dimension) * macro_entity_volume * fine_scale_average);
      }
    }

    // delete?
    // #ifdef TFR
    // delete[] corrector_Phi;
    // #endif

    number_of_macro_entity += 1;
  }

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
} // assemble_jacobian_matrix

// ! ------------------------------------------------------------------------------------------------
// ! ------------------------------------------------------------------------------------------------
}

#endif // #ifndef DiscreteElliptic_HH
