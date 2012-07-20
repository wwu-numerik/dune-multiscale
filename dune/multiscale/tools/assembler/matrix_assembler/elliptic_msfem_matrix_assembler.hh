#ifndef DiscreteEllipticMSFEMOperator_HH
#define DiscreteEllipticMSFEMOperator_HH

#include <dune/common/fmatrix.hh>

#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/multiscale/tools/solver/MsFEM/msfem_localproblems/subgrid-list.hh>
#include <dune/multiscale/tools/solver/MsFEM/msfem_localproblems/localproblemsolver.hh>

#include <dune/fem/operator/2order/lagrangematrixsetup.hh>

// / done

namespace Dune {
// Imp stands for Implementation
template< class CoarseDiscreteFunctionImp, class MacroMicroGridSpecifierImp, class FineDiscreteFunctionImp,
          class DiffusionImp >
class DiscreteEllipticMsFEMOperator
  : public Operator< typename CoarseDiscreteFunctionImp::RangeFieldType,
                     typename CoarseDiscreteFunctionImp::RangeFieldType,
                     CoarseDiscreteFunctionImp, CoarseDiscreteFunctionImp >
{
  typedef DiscreteEllipticMsFEMOperator< CoarseDiscreteFunctionImp, MacroMicroGridSpecifierImp, FineDiscreteFunctionImp,
                                         DiffusionImp > This;

public:
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

  typedef SubGrid< WORLDDIM, FineGrid >                                                 SubGridType;
  typedef SubGridList< FineDiscreteFunction, SubGridType, MacroMicroGridSpecifierType > SubGridListType;

  typedef MsFEMLocalProblemSolver< FineDiscreteFunction, SubGridListType, MacroMicroGridSpecifierType,
                                   DiffusionModel > MsFEMLocalProblemSolverType;

protected:
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

public:
  typedef typename CoarseDiscreteFunctionSpace::GridPartType CoarseGridPart;
  typedef typename CoarseDiscreteFunctionSpace::GridType     CoarseGrid;

protected:
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

public:
  // ! ---------------- typedefs for the SubgridDiscreteFunctionSpace -----------------------
  // ( typedefs for the local grid and the corresponding local ('sub') )discrete space )

  // ! type of grid part
  typedef LeafGridPart< SubGridType > SubGridPart;

  // ! type of subgrid discrete function space
  typedef LagrangeDiscreteFunctionSpace< FunctionSpace, SubGridPart, 1 >  // 1=POLORDER
  LocalDiscreteFunctionSpace;

  // ! type of subgrid discrete function
  typedef AdaptiveDiscreteFunction< LocalDiscreteFunctionSpace > LocalDiscreteFunction;

  typedef typename LocalDiscreteFunctionSpace::IteratorType LocalGridIterator;

  typedef typename LocalGridIterator::Entity LocalGridEntity;

  typedef typename LocalGridEntity::EntityPointer LocalGridEntityPointer;

  typedef typename LocalDiscreteFunction::LocalFunctionType LocalGridLocalFunction;

  typedef typename LocalDiscreteFunctionSpace::LagrangePointSetType LGLagrangePointSet;

  typedef typename LocalDiscreteFunctionSpace::BaseFunctionSetType LocalGridBaseFunctionSet;

  typedef typename LocalGridEntity::Geometry LocalGridGeometry;

  typedef CachingQuadrature< SubGridPart, 0 > LocalGridQuadrature;

  // !-----------------------------------------------------------------------------------------

public:
  DiscreteEllipticMsFEMOperator(MacroMicroGridSpecifierType& specifier,
                                const CoarseDiscreteFunctionSpace& coarseDiscreteFunctionSpace,
                                // number of layers per coarse grid entity T:  U(T) is created by enrichting T with
                                // n(T)-layers:
                                SubGridListType& subgrid_list,
                                const DiffusionModel& diffusion_op,
                                std::ofstream& data_file,
                                std::string path = "")
    : specifier_(specifier)
      , coarseDiscreteFunctionSpace_(coarseDiscreteFunctionSpace)
      , subgrid_list_(subgrid_list)
      , diffusion_operator_(diffusion_op)
      , data_file_(&data_file)
      , path_(path) {
    bool silence = false;

    // coarseDiscreteFunctionSpace_ = specifier_.coarseSpace();
    // fineDiscreteFunctionSpace_ = specifier_.fineSpace();

    std::string local_path = path_ + "/local_problems/";

    MsFEMLocalProblemSolverType loc_prob_solver(
      specifier_.fineSpace(), specifier_, subgrid_list_, diffusion_operator_, data_file, local_path);

    loc_prob_solver.assemble_all(silence);
  }

private:
  DiscreteEllipticMsFEMOperator(const This&);

public:
  // dummy operator
  virtual void operator()(const CoarseDiscreteFunction& u, CoarseDiscreteFunction& w) const;

  template< class MatrixType >
  void assemble_matrix(MatrixType& global_matrix) const;

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

private:
  // create a hostgrid function from a subgridfunction
  // Note: the maximum gride levels for both underlying grids must be the same
  void subgrid_to_hostrid_function(const LocalDiscreteFunction& sub_func, FineDiscreteFunction& host_func);

private:
  MacroMicroGridSpecifierType& specifier_;

  const CoarseDiscreteFunctionSpace& coarseDiscreteFunctionSpace_;
  SubGridListType& subgrid_list_;
  const DiffusionModel& diffusion_operator_;

  // data file for saving information
  std::ofstream* data_file_;

  // path where to save the data output
  std::string path_;
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
    std::cout
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

    const unsigned int numBaseFunctions = sub_loc_value.baseFunctionSet().numBaseFunctions();
    for (unsigned int i = 0; i < numBaseFunctions; ++i)
    {
      host_loc_value[i] = sub_loc_value[i];
    }
  }
} // subgrid_to_hostrid_function

// dummy implementation of "operator()"
// 'w' = effect of the discrete operator on 'u'
template< class CoarseDiscreteFunctionImp, class MacroMicroGridSpecifierImp, class FineDiscreteFunctionImp,
          class DiffusionImp >
void DiscreteEllipticMsFEMOperator< CoarseDiscreteFunctionImp,
                                    MacroMicroGridSpecifierImp,
                                    FineDiscreteFunctionImp,
                                    DiffusionImp >::operator()(const CoarseDiscreteFunction& /*u*/,
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
  #ifdef PGF
  DSC_LOG_INFO << "Assembling Petrov-Galerkin-MsFEM Matrix." << std::endl;
  #else
  DSC_LOG_INFO << "Assembling MsFEM Matrix." << std::endl;
  #endif // ifdef PGF

  typedef typename MatrixType::LocalMatrixType LocalMatrix;

  global_matrix.reserve();
  global_matrix.clear();

  std::vector< typename CoarseBaseFunctionSet::JacobianRangeType > gradient_Phi(
    coarseDiscreteFunctionSpace_.mapper().maxNumDofs() );

  const CoarseGridLeafIndexSet& coarseGridLeafIndexSet = coarseDiscreteFunctionSpace_.gridPart().grid().leafIndexSet();

  // Coarse Entity Iterator
  const CoarseIterator coarse_grid_end = coarseDiscreteFunctionSpace_.end();
  for (CoarseIterator coarse_grid_it = coarseDiscreteFunctionSpace_.begin();
       coarse_grid_it != coarse_grid_end;
       ++coarse_grid_it)
  {
    // the coarse grid element T:
    const CoarseEntity& coarse_grid_entity = *coarse_grid_it;
    const CoarseGeometry& coarse_grid_geometry = coarse_grid_entity.geometry();
    assert(coarse_grid_entity.partitionType() == InteriorEntity);

    int global_index_entity = coarseGridLeafIndexSet.index(coarse_grid_entity);

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

    char location_lps[50];
    sprintf(location_lps, "/local_problems/_localProblemSolutions_%d", global_index_entity);
    std::string location_lps_s(location_lps);

    std::string local_solution_location;

    // the file/place, where we saved the solutions of the cell problems
    local_solution_location = path_ + location_lps_s;

    bool reader_is_open = false;
    // reader for the cell problem data file:
    DiscreteFunctionReader discrete_function_reader( (local_solution_location).c_str() );
    reader_is_open = discrete_function_reader.open();

    if (reader_is_open)
    {
      discrete_function_reader.read(0, local_problem_solution_e0);
    }

    if (reader_is_open)
    {
      discrete_function_reader.read(1, local_problem_solution_e1);
    }

    // 1 point quadrature!! We only need the gradient of the base function,
    // which is constant on the whole entity.
    CoarseQuadrature one_point_quadrature(coarse_grid_entity, 0);

    // the barycenter of the macro_grid_entity
    const typename CoarseQuadrature::CoordinateType& local_coarse_point
      = one_point_quadrature.point(0 /*=quadraturePoint*/);
    DomainType DUNE_UNUSED(coarse_entity_barycenter) = coarse_grid_geometry.global(local_coarse_point);

    // transposed of the the inverse jacobian
    const FieldMatrix< double, dimension, dimension >& inverse_jac
      = coarse_grid_geometry.jacobianInverseTransposed(local_coarse_point);

    for (unsigned int i = 0; i < numMacroBaseFunctions; ++i)
    {
      // jacobian of the base functions, with respect to the reference element
      typename CoarseBaseFunctionSet::JacobianRangeType gradient_Phi_ref_element;
      coarse_grid_baseSet.jacobian(i, one_point_quadrature[0], gradient_Phi_ref_element);

      // multiply it with transpose of jacobian inverse to obtain the jacobian with respect to the real entity
      inverse_jac.mv(gradient_Phi_ref_element[0], gradient_Phi[i][0]);
    }

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

          FineEntityPointer father_of_loc_grid_ent = localDiscreteFunctionSpace.grid().template getHostEntity< 0 >(
            local_grid_entity);

          for (int lev = 0; lev < specifier_.getLevelDifference(); ++lev)
            father_of_loc_grid_ent = father_of_loc_grid_ent->father();

          FineEntityPointer coarse_father_test = father_of_loc_grid_ent;

          bool father_found = false;
          while (father_found == false)
          {
            if (coarseGridLeafIndexSet.contains(*coarse_father_test) == true)
            {
              father_of_loc_grid_ent = coarse_father_test;
            }

            if (coarse_father_test->hasFather() == false)
            {
              father_found = true;
            } else {
              coarse_father_test = coarse_father_test->father();
            }
          }

          bool entities_identical = true;
          int number_of_nodes = (*coarse_grid_it).template count< 2 >();
          for (int k = 0; k < number_of_nodes; k += 1)
          {
            if ( !( coarse_grid_it->geometry().corner(k) == father_of_loc_grid_ent->geometry().corner(k) ) )
            {
              entities_identical = false;
            }
          }

          if (entities_identical == false)
          {
            // std :: cout << "coarse_grid_it->geometry().corner(0) = " << coarse_grid_it->geometry().corner(0) << std
            // :: endl;
            // std :: cout << "coarse_grid_it->geometry().corner(1) = " << coarse_grid_it->geometry().corner(1) << std
            // :: endl;
            // std :: cout << "coarse_grid_it->geometry().corner(2) = " << coarse_grid_it->geometry().corner(2) << std
            // :: endl;
            // std :: cout << "father_of_loc_grid_ent->geometry().corner(0) = " <<
            // father_of_loc_grid_ent->geometry().corner(0) << std :: endl;
            // std :: cout << "father_of_loc_grid_ent->geometry().corner(1) = " <<
            // father_of_loc_grid_ent->geometry().corner(1) << std :: endl;
            // std :: cout << "father_of_loc_grid_ent->geometry().corner(2) = " <<
            // father_of_loc_grid_ent->geometry().corner(2) << std :: endl << std :: endl;
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

            // if not Petrov-Galerkin:
            #ifndef PGF
            JacobianRangeType reconstruction_grad_phi_j(0.0);
            for (int k = 0; k < dimension; ++k)
            {
              reconstruction_grad_phi_j[0][k] += gradient_Phi[j][0][0] * grad_loc_sol_e0[0][k];
              reconstruction_grad_phi_j[0][k] += gradient_Phi[j][0][1] * grad_loc_sol_e1[0][k];
              reconstruction_grad_phi_j[0][k] += gradient_Phi[j][0][k];
            }

            local_integral += weight_local_quadrature * (diffusive_flux[0] * reconstruction_grad_phi_j[0]);
            #else // ifndef PGF
            local_integral += weight_local_quadrature * (diffusive_flux[0] * gradient_Phi[j][0]);
            #endif // ifndef PGF
          }
        }

        // add entries
        local_matrix.add(j, i, local_integral);
      }
    }
  }

  // discrete_function_reader.close();

  // boundary treatment
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

// ! ------------------------------------------------------------------------------------------------
// ! ------------------------------------------------------------------------------------------------
}

#endif // #ifndef DiscreteElliptic_HH
