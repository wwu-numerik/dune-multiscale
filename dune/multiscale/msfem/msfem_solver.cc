#ifdef HAVE_CMAKE_CONFIG
 #include "cmake_config.h"
#elif defined (HAVE_CONFIG_H)
 #include <config.h>
#endif // ifdef HAVE_CMAKE_CONFIG

#include <unordered_set>
#include "msfem_solver.hh"

#include <dune/multiscale/common/righthandside_assembler.hh>
#include <dune/multiscale/msfem/localproblems/subgrid-list.hh>
#include <dune/multiscale/msfem/elliptic_msfem_matrix_assembler.hh>
#include <dune/multiscale/tools/misc/linear-lagrange-interpolation.hh>
#include <dune/multiscale/msfem/msfem_grid_specifier.hh>
#include <dune/multiscale/msfem/elliptic_msfem_matrix_assembler.hh>
#include <dune/stuff/discretefunction/projection/heterogenous.hh>

namespace Dune {
namespace Multiscale {
namespace MsFEM {

Elliptic_MsFEM_Solver::Elliptic_MsFEM_Solver(const DiscreteFunctionSpace& discreteFunctionSpace)
  : discreteFunctionSpace_(discreteFunctionSpace)
{}

void Elliptic_MsFEM_Solver::subgrid_to_hostrid_projection(const SubgridDiscreteFunction& sub_func, DiscreteFunction& host_func) const {
  host_func.clear();

  const SubgridDiscreteFunctionSpace& subDiscreteFunctionSpace = sub_func.space();
  const SubGridType& subGrid = subDiscreteFunctionSpace.grid();

  typedef typename SubgridDiscreteFunctionSpace::IteratorType SubgridIterator;
  typedef typename SubgridIterator::Entity                    SubgridEntity;
  typedef typename SubgridDiscreteFunction::LocalFunctionType SubgridLocalFunction;

  const SubgridIterator sub_endit = subDiscreteFunctionSpace.end();
  for (SubgridIterator sub_it = subDiscreteFunctionSpace.begin(); sub_it != sub_endit; ++sub_it)
  {
    const SubgridEntity& sub_entity = *sub_it;

    const HostEntityPointer host_entity_pointer = subGrid.getHostEntity< 0 >(*sub_it);
    const HostEntity& host_entity = *host_entity_pointer;

    const SubgridLocalFunction sub_loc_value = sub_func.localFunction(sub_entity);
    LocalFunction host_loc_value = host_func.localFunction(host_entity);

    const auto numBaseFunctions = sub_loc_value.baseFunctionSet().size();
    for (unsigned int i = 0; i < numBaseFunctions; ++i)
    {
      host_loc_value[i] = sub_loc_value[i];
    }
  }
} // subgrid_to_hostrid_projection

void Elliptic_MsFEM_Solver::identify_coarse_scale_part( MacroMicroGridSpecifier& specifier,
                                 const DiscreteFunction& coarse_msfem_solution,
                                 DiscreteFunction& coarse_scale_part ) const
{

  DSC_LOG_INFO << "Indentifying coarse scale part of the MsFEM solution... ";

  coarse_scale_part.clear();
  Dune::Stuff::HeterogenousProjection<> projection;
  projection.project(coarse_msfem_solution, coarse_scale_part);
  DSC_LOG_INFO << " done." << std::endl;
}


void Elliptic_MsFEM_Solver::identify_fine_scale_part( MacroMicroGridSpecifier& specifier,
                                                        MsFEMTraits::SubGridListType& subgrid_list,
                                                        const DiscreteFunction& coarse_msfem_solution,
                                                        DiscreteFunction& fine_scale_part ) const
{

  fine_scale_part.clear();

  const HostGrid& grid = discreteFunctionSpace_.gridPart().grid();
  const GridPart& gridPart = discreteFunctionSpace_.gridPart();

  DiscreteFunctionSpace& coarse_space = specifier.coarseSpace();
  const HostGridLeafIndexSet& coarseGridLeafIndexSet = coarse_space.gridPart().grid().leafIndexSet();

  const int number_of_nodes = grid.size(2 /*codim*/);
  std::vector< std::vector< HostEntityPointer > > entities_sharing_same_node(number_of_nodes);

  // if the oversampling strategy is 1 or 2, we need identify the entities that share a certain node
  // (because a 'conforming projection' operator is required - for stragey 3, we directly get a conforming approximation)
  if ( ( specifier.getOversamplingStrategy() == 1 ) || ( specifier.getOversamplingStrategy() == 2 ) ) {
    // determine the entities that share a common global node with a given index
    // we need to iterate over the whole grid, not only from hostSpace_.begin() to
    // hostSpace_.end() for parallel runs!
    for (auto& hostEntity : DSC::viewRange(discreteFunctionSpace_.gridPart().grid().leafView())) {
      int number_of_nodes_in_entity = hostEntity.count< 2 >();
      for (int i = 0; i < number_of_nodes_in_entity; ++i) {
        const auto node              = hostEntity.subEntity< 2 >(i);
        const int             global_index_node = gridPart.indexSet().index(*node);

        entities_sharing_same_node[global_index_node].emplace_back(hostEntity);
      }
    }


  }

  DSC_LOG_INFO << "Indentifying fine scale part of the MsFEM solution... ";
  // traverse coarse space
  for (HostgridIterator coarse_it = coarse_space.begin(); coarse_it != coarse_space.end(); ++coarse_it) {
    // the coarse entity 'T': *coarse_it

    // only required for oversampling strategy 1 and 2, where we need to identify the correction for each
    DiscreteFunction correction_on_U_T("correction_on_U_T", discreteFunctionSpace_);
    correction_on_U_T.clear();

    const int index = coarseGridLeafIndexSet.index(*coarse_it);

    // the sub-grid U(T) that belongs to the coarse_grid_entity T
    SubGridType& sub_grid_U_T = subgrid_list.getSubGrid(index);
    SubGridPart subGridPart(sub_grid_U_T);

    const SubgridDiscreteFunctionSpace localDiscreteFunctionSpace(subGridPart);

    SubgridDiscreteFunction local_problem_solution_e0("Local problem Solution e_0", localDiscreteFunctionSpace);
    local_problem_solution_e0.clear();

    SubgridDiscreteFunction local_problem_solution_e1("Local problem Solution e_1", localDiscreteFunctionSpace);
    local_problem_solution_e1.clear();

    // --------- load local solutions -------
    // the file/place, where we saved the solutions of the cell problems
    const std::string local_solution_location = (boost::format("local_problems/_localProblemSolutions_%d_%d")
                                                % index % MPIManager::rank()).str();
    // reader for the cell problem data file:
    DiscreteFunctionReader discrete_function_reader(local_solution_location);
    discrete_function_reader.read(0, local_problem_solution_e0);
    discrete_function_reader.read(1, local_problem_solution_e1);

    LocalFunction local_coarse_part = coarse_msfem_solution.localFunction(*coarse_it);

    // 1 point quadrature!! We only need the gradient of the coarse scale part on the element, which is a constant.
    const CachingQuadrature< GridPart, 0 > one_point_quadrature(*coarse_it, 0);

    JacobianRangeType grad_coarse_msfem_on_entity;
    local_coarse_part.jacobian(one_point_quadrature[0], grad_coarse_msfem_on_entity);

    //!
    // std :: cout << "grad_coarse_msfem_on_entity[ 0 ][ 1 ] = " << grad_coarse_msfem_on_entity[ 0 ][ 1 ] << std ::
    // endl;
    // std :: cout << "grad_coarse_msfem_on_entity[ 0 ][ 0 ] = " << grad_coarse_msfem_on_entity[ 0 ][ 0 ] << std ::
    // endl;
    // get the coarse gradient on T, multiply it with the local correctors and sum it up.
    local_problem_solution_e0 *= grad_coarse_msfem_on_entity[0][0];
    local_problem_solution_e1 *= grad_coarse_msfem_on_entity[0][1];
    local_problem_solution_e0 += local_problem_solution_e1;

    // oneLinePrint( DSC_LOG_DEBUG, local_problem_solution_e0 );

    // oversampling strategy 3: just sum up the local correctors:
    if ( (specifier.getOversamplingStrategy() == 3) ) {
      subgrid_to_hostrid_projection(local_problem_solution_e0, correction_on_U_T);
    }

    // oversampling strategy 1 or 2: restrict the local correctors to the element T, sum them up and apply a conforming projection:
    if ( ( specifier.getOversamplingStrategy() == 1 ) || ( specifier.getOversamplingStrategy() == 2 ) ) {

      if ( sub_grid_U_T.maxLevel() != discreteFunctionSpace_.gridPart().grid().maxLevel() ) {
        DSC_LOG_ERROR << "Error: MaxLevel of SubGrid not identical to MaxLevel of FineGrid." << std::endl;
      }

      correction_on_U_T.clear();

      typedef typename SubgridDiscreteFunctionSpace::IteratorType SubgridIterator;
      typedef typename SubgridIterator::Entity                    SubgridEntity;
      typedef typename SubgridDiscreteFunction::LocalFunctionType SubgridLocalFunction;

      const SubgridIterator sub_endit = localDiscreteFunctionSpace.end();
      for (SubgridIterator sub_it = localDiscreteFunctionSpace.begin(); sub_it != sub_endit; ++sub_it) {
        const SubgridEntity& sub_entity = *sub_it;

        const HostEntityPointer fine_host_entity_pointer = sub_grid_U_T.getHostEntity< 0 >(*sub_it);
        const HostEntity& fine_host_entity = *fine_host_entity_pointer;

        const int hostFatherIndex = subgrid_list.getEnclosingMacroCellIndex(fine_host_entity_pointer);
        if (hostFatherIndex==index) {

          const SubgridLocalFunction sub_loc_value = local_problem_solution_e0.localFunction(sub_entity);
          LocalFunction host_loc_value = correction_on_U_T.localFunction(fine_host_entity);

          int number_of_nodes_entity = sub_it->count< 2 >();
          for (int i = 0; i < number_of_nodes_entity; ++i)
          {
            const typename HostEntity::Codim< 2 >::EntityPointer node = fine_host_entity.subEntity< 2 >(i);

            const int global_index_node = gridPart.indexSet().index(*node);

            // count the number of different coarse-grid-entities that share the above node
            std::unordered_set< int > coarse_entities;
            for (size_t j = 0; j < entities_sharing_same_node[global_index_node].size(); ++j) {
              const int innerIndex
                      = subgrid_list.getEnclosingMacroCellIndex(entities_sharing_same_node[global_index_node][j]);
              // the following will only add the entity index if it is not yet present
              coarse_entities.emplace(innerIndex);
            }

            host_loc_value[i] = ( sub_loc_value[i] / coarse_entities.size() );
          }
        }
      }
    }

    fine_scale_part += correction_on_U_T;
  }
  DSC_LOG_INFO << " done." << std::endl;
}


void Elliptic_MsFEM_Solver::solve_dirichlet_zero(const CommonTraits::DiffusionType& diffusion_op,
                          const CommonTraits::FirstSourceType& f,
                          // number of layers per coarse grid entity T:  U(T) is created by enrichting T with
                          // n(T)-layers.
                          MacroMicroGridSpecifier& specifier,
                          MsFEMTraits::SubGridListType& subgrid_list,
                          DiscreteFunction& coarse_scale_part,
                          DiscreteFunction& fine_scale_part,
                          DiscreteFunction& solution) const
{
  DSC::Profiler::ScopedTiming st("msfem.Elliptic_MsFEM_Solver.solve_dirichlet_zero");

  DiscreteFunctionSpace& coarse_space = specifier.coarseSpace();

  DiscreteFunction coarse_msfem_solution("Coarse Part MsFEM Solution", coarse_space);
  coarse_msfem_solution.clear();

  //! define the right hand side assembler tool
  // (for linear and non-linear elliptic and parabolic problems, for sources f and/or G )
  typedef RightHandSideAssembler< DiscreteFunction > RhsAssembler;

  //! define the discrete (elliptic) operator that describes our problem
  // discrete elliptic MsFEM operator (corresponds with MsFEM Matrix)
  // ( effect of the discretized differential operator on a certain discrete function )
  const DiscreteEllipticMsFEMOperator elliptic_msfem_op(specifier,
                                              coarse_space,
                                              subgrid_list,
                                              diffusion_op);
  // discrete elliptic operator (corresponds with FEM Matrix)

  //! (stiffness) matrix
  MsFEMMatrix msfem_matrix("MsFEM stiffness matrix", coarse_space, coarse_space);

  //! right hand side vector
  // right hand side for the finite element method:
  DiscreteFunction msfem_rhs("MsFEM right hand side", coarse_space);
  msfem_rhs.clear();

  DSC_LOG_INFO  << std::endl << "Solving MsFEM problem." << std::endl;
  DSC_LOG_INFO << "Solving linear problem with MsFEM and maximum coarse grid level "
              << coarse_space.gridPart().grid().maxLevel() << "." << std::endl;
  DSC_LOG_INFO << "------------------------------------------------------------------------------" << std::endl;

  // to assemble the computational time
  Dune::Timer assembleTimer;

  // assemble the MsFEM stiffness matrix
  elliptic_msfem_op.assemble_matrix(msfem_matrix);
  DSC_LOG_INFO << "Time to assemble MsFEM stiffness matrix: " << assembleTimer.elapsed() << "s" << std::endl;

  // assemble right hand side
  if ( DSC_CONFIG_GET("msfem.petrov_galerkin", 1 ) ) {
    RhsAssembler::assemble< 2* DiscreteFunctionSpace::polynomialOrder + 2 >(f, msfem_rhs);
  } else {
    RhsAssembler::assemble_for_MsFEM_symmetric< 2* DiscreteFunctionSpace::polynomialOrder + 2 >(f,
                                                                                                specifier,
                                                                                                subgrid_list,
                                                                                                msfem_rhs);
  }

  // oneLinePrint( DSC_LOG_DEBUG, fem_rhs );

  //! --- boundary treatment ---
  // set the dirichlet points to zero (in right hand side of the fem problem)
  const HostgridIterator endit = coarse_space.end();
  for (HostgridIterator it = coarse_space.begin(); it != endit; ++it) {
    IntersectionIterator iit = coarse_space.gridPart().ibegin(*it);
    const IntersectionIterator endiit = coarse_space.gridPart().iend(*it);
    for ( ; iit != endiit; ++iit) {
      if ( iit->boundary() ) {
        LocalFunction rhsLocal = msfem_rhs.localFunction(*it);

        const LagrangePointSet& lagrangePointSet
                = coarse_space.lagrangePointSet(*it);

        const int face = iit->indexInInside();

        auto faceIterator
                = lagrangePointSet.beginSubEntity< faceCodim >(face);
        const auto faceEndIterator
                = lagrangePointSet.endSubEntity< faceCodim >(face);
        for ( ; faceIterator != faceEndIterator; ++faceIterator)
          rhsLocal[*faceIterator] = 0;
      }
    }
  }
  msfem_rhs.communicate();
  //! --- end boundary treatment ---
  const InverseMsFEMMatrix msfem_biCGStab(msfem_matrix, 1e-8, 1e-8, 2000, true);
  msfem_biCGStab(msfem_rhs, coarse_msfem_solution);
  DSC_LOG_INFO << "---------------------------------------------------------------------------------" << std::endl;
  DSC_LOG_INFO << "MsFEM problem solved in " << assembleTimer.elapsed() << "s." << std::endl << std::endl
               << std::endl;

  // oneLinePrint( DSC_LOG_DEBUG, solution );
  // copy coarse grid function (defined on the subgrid) into a fine grid function
  solution.clear();

  //! copy coarse scale part of MsFEM solution into a function defined on the fine grid
  identify_coarse_scale_part( specifier, coarse_msfem_solution, coarse_scale_part );

  //! identify fine scale part of MsFEM solution (including the projection!)
  identify_fine_scale_part( specifier, subgrid_list, coarse_msfem_solution, fine_scale_part );

  // add coarse and fine scale part to solution
  solution += coarse_scale_part;
  solution += fine_scale_part;

  solution.communicate();
} // solve_dirichlet_zero

} //namespace MsFEM {
} //namespace Multiscale {
} //namespace Dune {
