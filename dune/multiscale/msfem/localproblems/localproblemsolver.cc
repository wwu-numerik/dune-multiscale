#ifdef HAVE_CMAKE_CONFIG
 #include "cmake_config.h"
#elif defined (HAVE_CONFIG_H)
 #include <config.h>
#endif // ifdef HAVE_CMAKE_CONFIG

#include "localproblemsolver.hh"

#include <dune/subgrid/subgrid.hh>
#include <dune/stuff/common/filesystem.hh>
#include <dune/stuff/fem/functions/checks.hh>
#include <dune/stuff/common/profiler.hh>
#include <dune/stuff/common/memory.hh>

#include <dune/multiscale/tools/subgrid_io.hh>
#include <dune/multiscale/hmm/cell_problem_numbering.hh>
#include <dune/multiscale/msfem/msfem_traits.hh>
#include <dune/multiscale/msfem/localproblems/localoperator.hh>

#include <dune/multiscale/tools/misc/uzawa.hh>
#include <dune/multiscale/tools/misc/weighted-clement-operator.hh>
#include <dune/multiscale/tools/discretefunctionwriter.hh>

#include <memory>
#include <vector>

namespace Dune {
namespace Multiscale {
namespace MsFEM {

LocalProblemDataOutputParameters::LocalProblemDataOutputParameters()
  : OutputParameters(DSC_CONFIG_GET("global.datadir", "data") + "/local_problems/")
{}

MsFEMLocalProblemSolver::MsFEMLocalProblemSolver(const HostDiscreteFunctionSpaceType& hostDiscreteFunctionSpace,
                        const MacroMicroGridSpecifierType& specifier,
                        SubGridList& subgrid_list,
                        const DiffusionOperatorType& diffusion_operator)
  : hostDiscreteFunctionSpace_(hostDiscreteFunctionSpace)
  , diffusion_(diffusion_operator)
  , specifier_(specifier)
  , subgrid_list_(subgrid_list)
  , ids_basis_functions_in_subgrid_(nullptr)
  , inverse_of_L1_norm_coarse_basis_funcs_(nullptr)
  , coarse_basis_(nullptr)
  , global_id_to_internal_id_(nullptr)
{}

MsFEMLocalProblemSolver::MsFEMLocalProblemSolver(const HostDiscreteFunctionSpaceType& hostDiscreteFunctionSpace,
                        const MacroMicroGridSpecifierType& specifier,
                        SubGridList& subgrid_list, std::vector< std::vector< int > >& ids_basis_functions_in_subgrid,
                        std::vector< double >& inverse_of_L1_norm_coarse_basis_funcs,
                        const DiffusionOperatorType& diffusion_operator,
                        const CoarseBasisFunctionListType& coarse_basis,
                        const std::map<int,int>& global_id_to_internal_id )
  : hostDiscreteFunctionSpace_(hostDiscreteFunctionSpace)
  , diffusion_(diffusion_operator)
  , specifier_(specifier)
  , subgrid_list_(subgrid_list)
  , ids_basis_functions_in_subgrid_( &ids_basis_functions_in_subgrid ) // ids of the coarse grid basis functions in the interior of the subgrid
  , inverse_of_L1_norm_coarse_basis_funcs_( &inverse_of_L1_norm_coarse_basis_funcs )
  , coarse_basis_( &coarse_basis )
  , global_id_to_internal_id_( &global_id_to_internal_id )
{}


//! ----------- method: solve the local MsFEM problem ------------------------------------------

void MsFEMLocalProblemSolver::solvelocalproblem(JacobianRangeType& e,
                       SubDiscreteFunctionType& local_problem_solution,
                       const int coarse_index /*= -1*/ ) const {

    typedef SparseRowMatrixTraits < SubDiscreteFunctionSpaceType, HostDiscreteFunctionSpaceType >
        WeightedClementMatrixObjectTraits;

    typedef WeightedClementOp< SubDiscreteFunctionType, HostDiscreteFunctionType, WeightedClementMatrixObjectTraits, CoarseBasisFunctionListType >
              WeightedClementOperatorType;

    // saddle point problem solver:
    typedef UzawaInverseOp< SubDiscreteFunctionType,
                            HostDiscreteFunctionType,
                            InverseLocProbFEMMatrix,
                            WeightedClementOperatorType >
       InverseUzawaOperatorType;

  // set solution equal to zero:
  local_problem_solution.clear();

  const SubDiscreteFunctionSpaceType& subDiscreteFunctionSpace = local_problem_solution.space();

  //! the matrix in our linear system of equations
  // in the non-linear case, it is the matrix for each iteration step
  LocProbFEMMatrix locprob_system_matrix("Local Problem System Matrix",
                                         subDiscreteFunctionSpace,
                                         subDiscreteFunctionSpace);

  //! define the discrete (elliptic) local MsFEM problem operator
  // ( effect of the discretized differential operator on a certain discrete function )
  LocalProblemOperator local_problem_op(subDiscreteFunctionSpace, diffusion_);

  const SubGridType& subGrid = subDiscreteFunctionSpace.grid();

  typedef typename SubDiscreteFunctionSpaceType::IteratorType SGIteratorType;
  typedef typename SubGridPartType::IntersectionIteratorType  SGIntersectionIteratorType;

  //! right hand side vector of the algebraic local MsFEM problem
  SubDiscreteFunctionType local_problem_rhs("rhs of local MsFEM problem", subDiscreteFunctionSpace);
  local_problem_rhs.clear();

  // NOTE:
  // is the right hand side of the local MsFEM problem equal to zero or almost identical to zero?
  // if yes, the solution of the local MsFEM problem is also identical to zero. The solver is getting a problem with
  // this situation, which is why we do not solve local msfem problems for zero-right-hand-side, since we already know
  // the result.

  switch ( specifier_.getOversamplingStrategy() )
  {
  case 1:
    local_problem_op.assemble_matrix(locprob_system_matrix);
    local_problem_op.assemble_local_RHS(e, local_problem_rhs);
    break;
  case 2:
    if ( coarse_index < 0 )
      DUNE_THROW(Dune::InvalidStateException, "Invalid coarse index: coarse_index < 0");
    local_problem_op.assemble_matrix(locprob_system_matrix, subgrid_list_.getCoarseNodeVector( coarse_index ) );
    local_problem_op.assemble_local_RHS(e,
            subgrid_list_.getCoarseNodeVector( coarse_index ),
            specifier_.getOversamplingStrategy(),
            local_problem_rhs );
    break;
  case 3:
    if ( coarse_index < 0 )
      DUNE_THROW(Dune::InvalidStateException, "Invalid coarse index: coarse_index < 0");
    if ( DSC_CONFIG_GET( "rigorous_msfem.oversampling_strategy", "Clement" ) == "Clement" ) {
      local_problem_op.assemble_matrix( locprob_system_matrix );
    } else {
      local_problem_op.assemble_matrix( locprob_system_matrix, subgrid_list_.getCoarseNodeVector( coarse_index ) );
    }
    local_problem_op.assemble_local_RHS(e,
            subgrid_list_.getCoarseNodeVector( coarse_index ),
            specifier_.getOversamplingStrategy(),
            local_problem_rhs );
    break;    
  default: DUNE_THROW(Dune::InvalidStateException, "Oversampling Strategy must be 1, 2 or 3.");
  }

  //! boundary treatment:
  typedef typename LocProbFEMMatrix::LocalMatrixType LocalMatrix;

  typedef typename SGLagrangePointSetType::Codim< faceCodim >::SubEntityIteratorType
      FaceDofIteratorType;

  const HostGridPartType& hostGridPart = hostDiscreteFunctionSpace_.gridPart();

  const SubgridIteratorType sg_end = subDiscreteFunctionSpace.end();
  for (SubgridIteratorType sg_it = subDiscreteFunctionSpace.begin(); sg_it != sg_end; ++sg_it)
  {
    const SubgridEntityType& subgrid_entity = *sg_it;

    HostEntityPointerType host_entity_pointer = subGrid.getHostEntity< 0 >(subgrid_entity);
    const HostEntityType& host_entity = *host_entity_pointer;

    LocalMatrix local_matrix = locprob_system_matrix.localMatrix(subgrid_entity, subgrid_entity);

    const SGLagrangePointSetType& lagrangePointSet = subDiscreteFunctionSpace.lagrangePointSet(subgrid_entity);

    const HostIntersectionIterator iend = hostGridPart.iend(host_entity);
    for (HostIntersectionIterator iit = hostGridPart.ibegin(host_entity); iit != iend; ++iit)
    {
      if ( iit->neighbor() ) // if there is a neighbor entity
      {
        // check if the neighbor entity is in the subgrid
        const HostEntityPointerType neighborHostEntityPointer = iit->outside();
        const HostEntityType& neighborHostEntity = *neighborHostEntityPointer;
        if ( subGrid.contains< 0 >(neighborHostEntity) )
        {
          continue;
        }
      }

      const int face = (*iit).indexInInside();
      const FaceDofIteratorType fdend = lagrangePointSet.endSubEntity< 1 >(face);
      for (FaceDofIteratorType fdit = lagrangePointSet.beginSubEntity< 1 >(face); fdit != fdend; ++fdit)
        local_matrix.unitRow(*fdit);
    }
  }


  // zero boundary condition for 'cell problems':
  // set Dirichlet Boundary to zero
  for (SubgridIteratorType sg_it = subDiscreteFunctionSpace.begin(); sg_it != sg_end; ++sg_it)
  {
    const SubgridEntityType& subgrid_entity = *sg_it;

    HostEntityPointerType host_entity_pointer = subGrid.getHostEntity< 0 >(subgrid_entity);
    const HostEntityType& host_entity = *host_entity_pointer;

    HostIntersectionIterator iit = hostGridPart.ibegin(host_entity);
    const HostIntersectionIterator endiit = hostGridPart.iend(host_entity);
    for ( ; iit != endiit; ++iit)
    {
      if ( iit->neighbor() ) // if there is a neighbor entity
      {
        // check if the neighbor entity is in the subgrid
        const HostEntityPointerType neighborHostEntityPointer = iit->outside();
        const HostEntityType& neighborHostEntity = *neighborHostEntityPointer;

        if ( subGrid.contains< 0 >(neighborHostEntity) )
        {
          continue;
        }
      }

      SubLocalFunctionType rhsLocal = local_problem_rhs.localFunction(subgrid_entity);
      const SGLagrangePointSetType& lagrangePointSet
          = subDiscreteFunctionSpace.lagrangePointSet(subgrid_entity);

      const int face = (*iit).indexInInside();

      FaceDofIteratorType faceIterator
          = lagrangePointSet.beginSubEntity< faceCodim >(face);
      const FaceDofIteratorType faceEndIterator
          = lagrangePointSet.endSubEntity< faceCodim >(face);
      for ( ; faceIterator != faceEndIterator; ++faceIterator)
        rhsLocal[*faceIterator] = 0;
    }
  }

  if ( !( local_problem_rhs.dofsValid() ) )
  {
    DUNE_THROW(Dune::InvalidStateException, "Local MsFEM Problem RHS invalid.");
  }

  if (local_problem_op.normRHS(local_problem_rhs) < /*1e-06*/ 1e-30)
  {
    local_problem_solution.clear();
    DSC_LOG_ERROR << "Local MsFEM problem with solution zero." << std::endl;
  }
  else
  {
    InverseLocProbFEMMatrix locprob_fem_biCGStab(locprob_system_matrix, 1e-8, 1e-8, 20000, DSC_CONFIG_GET("localproblemsolver_verbose", false));
    
    bool clement = false;
    if ( specifier_.getOversamplingStrategy() == 3 )
    { clement = (DSC_CONFIG_GET( "rigorous_msfem.oversampling_strategy", "Clement" ) == "Clement" ); }

    if ( clement )
    {
      HostDiscreteFunctionType zero("zero", specifier_.coarseSpace());
      zero.clear();
      const double dummy = 12345.67890;
      double solverEps = 1e-2;
      int maxIterations = 1000;

      // we want to solve the local problem with the constraint that the weighted Clement interpoltion
      // of the local problem solution is zero

      // implementation of a weighted Clement interpolation operator for our purpose:
      WeightedClementOperatorType clement_interpolation_op( subDiscreteFunctionSpace,
                                                            specifier_.coarseSpace(),
                                                            subgrid_list_.getCoarseNodeVector( coarse_index ),
                                                            *coarse_basis_, *global_id_to_internal_id_, specifier_ );
      //! NOTE TODO: implementation is not yet optimal, because the weighted Clement maps a function
      //! defined on the local subgrid to a function defined on the whole(!) coarse space.
      //! It would be better to implement a mapping to a localized coarse space, since
      //! the uzawa solver must treat ALL coarse grid nodes (expensive and worse convergence).

      //clement_interpolation_op.print();

      HostDiscreteFunctionType lagrange_multiplier("lagrange multiplier", specifier_.coarseSpace() );
      lagrange_multiplier.clear();

      // create inverse operator
      // saddle point problem solver with uzawa algorithm:
      {
        DSC::Profiler::ScopedTiming st("uzawa");
        InverseUzawaOperatorType uzawa( locprob_fem_biCGStab, clement_interpolation_op, dummy, solverEps, maxIterations, true);
        uzawa( local_problem_rhs, zero /*interpolation is zero*/, local_problem_solution, lagrange_multiplier );
      }
    }
    else {
      locprob_fem_biCGStab(local_problem_rhs, local_problem_solution);
    }
  }

  if ( !( local_problem_solution.dofsValid() ) ) {
    DUNE_THROW(Dune::InvalidStateException,"Current solution of the local msfem problem invalid!");
  }

  // oneLinePrint( DSC_LOG_DEBUG, local_problem_solution );
} // solvelocalproblem

// solve local problems for Local Orthogonal Decomposition Method (LOD) 
void MsFEMLocalProblemSolver::solvelocalproblems_lod(JacobianRangeType& e_0,
                       JacobianRangeType& e_1,
                       SubDiscreteFunctionType& local_problem_solution_0,
                       SubDiscreteFunctionType& local_problem_solution_1,
                       const int coarse_index /*= -1*/ ) const {

   //! if we do not sort out the coarse boundaries (on rigorous_msfem_solver.cc, line 175), the results get better

  // set solution equal to zero:
  local_problem_solution_0.clear();
  local_problem_solution_1.clear();

  if ( coarse_index < 0 )
      DUNE_THROW(Dune::InvalidStateException, "Invalid coarse index: coarse_index < 0");

  const SubDiscreteFunctionSpaceType& subDiscreteFunctionSpace = local_problem_solution_0.space();

  //! the matrix in our linear system of equations
  // in the non-linear case, it is the matrix for each iteration step
  LocProbFEMMatrix locprob_system_matrix("Local Problem System Matrix",
                                         subDiscreteFunctionSpace,
                                         subDiscreteFunctionSpace);

  //! define the discrete (elliptic) local MsFEM problem operator
  // ( effect of the discretized differential operator on a certain discrete function )
  LocalProblemOperator local_problem_op(subDiscreteFunctionSpace, diffusion_);

  const SubGridType& subGrid = subDiscreteFunctionSpace.grid();

  typedef typename SubDiscreteFunctionSpaceType::IteratorType SGIteratorType;
  typedef typename SubGridPartType::IntersectionIteratorType  SGIntersectionIteratorType;

  // NOTE:
  // is the right hand side of the local MsFEM problem equal to zero or almost identical to zero?
  // if yes, the solution of the local MsFEM problem is also identical to zero. The solver is getting a problem with
  // this situation, which is why we do not solve local msfem problems for zero-right-hand-side, since we already know
  // the result.
  
  switch ( specifier_.getOversamplingStrategy() )
  {
    case 3: break;
    default: DUNE_THROW(Dune::InvalidStateException, "method 'solvelocalproblems_lod' can be only used in combination with the LOD.");
  }

  bool clement = ( DSC_CONFIG_GET( "rigorous_msfem.oversampling_strategy", "Clement" ) == "Clement" );

  if ( !clement )
    DUNE_THROW(Dune::InvalidStateException, "method 'solvelocalproblems_lod' can be only used in combination with the LOD and Clement interpolation.");

Dune::Timer assembleTimer;  
  // assemble stiffness matrix
  local_problem_op.assemble_matrix( locprob_system_matrix );
DSC_LOG_INFO << "Time for assembling the locprob_system_matrix: " << assembleTimer.elapsed() << "s" << std::endl << std::endl; assembleTimer.reset();

  //! boundary treatment:
  // ----------------------------------------------------------------------------------------------------
  typedef typename LocProbFEMMatrix::LocalMatrixType LocalMatrix;

  typedef typename SGLagrangePointSetType::Codim< faceCodim >::SubEntityIteratorType
      FaceDofIteratorType;

  const HostGridPartType& hostGridPart = hostDiscreteFunctionSpace_.gridPart();

  const SubgridIteratorType sg_end = subDiscreteFunctionSpace.end();
  for (SubgridIteratorType sg_it = subDiscreteFunctionSpace.begin(); sg_it != sg_end; ++sg_it)
  {
    const SubgridEntityType& subgrid_entity = *sg_it;

    HostEntityPointerType host_entity_pointer = subGrid.getHostEntity< 0 >(subgrid_entity);
    const HostEntityType& host_entity = *host_entity_pointer;

    LocalMatrix local_matrix = locprob_system_matrix.localMatrix(subgrid_entity, subgrid_entity);

    const SGLagrangePointSetType& lagrangePointSet = subDiscreteFunctionSpace.lagrangePointSet(subgrid_entity);

    const HostIntersectionIterator iend = hostGridPart.iend(host_entity);
    for (HostIntersectionIterator iit = hostGridPart.ibegin(host_entity); iit != iend; ++iit)
    {
      if ( iit->neighbor() ) // if there is a neighbor entity
      {
        // check if the neighbor entity is in the subgrid
        const HostEntityPointerType neighborHostEntityPointer = iit->outside();
        const HostEntityType& neighborHostEntity = *neighborHostEntityPointer;
        if ( subGrid.contains< 0 >(neighborHostEntity) )
        {
          continue;
        }
      }

      const int face = (*iit).indexInInside();
      const FaceDofIteratorType fdend = lagrangePointSet.endSubEntity< 1 >(face);
      for (FaceDofIteratorType fdit = lagrangePointSet.beginSubEntity< 1 >(face); fdit != fdend; ++fdit)
        local_matrix.unitRow(*fdit);
    }
  }
  // ----------------------------------------------------------------------------------------------------
DSC_LOG_INFO << "Time for first boundary treatment in system matrix: " << assembleTimer.elapsed() << "s" << std::endl << std::endl; assembleTimer.reset();
  const InverseLocProbFEMMatrix locprob_inverse_system_matrix(locprob_system_matrix,
                                                        1e-8, 1e-8, 20000,
                                                        DSC_CONFIG_GET("localproblemsolver_verbose", false));

  //! Pre-processing step
  // For each coarse node j in the subgrid (local internal numbering, i.e. 0 <= j < M_subgrid), solve
  // for b_h_j with S_h b_h_j = C_h^T e_j, where C_h describes the algebraic version of the weighted
  // Clement interpolation operator
  // ----------------------------------------------------------------------------------------------------
  int number_of_interior_coarse_nodes_in_subgrid = (*ids_basis_functions_in_subgrid_)[ coarse_index ].size();
    
  std::vector<std::unique_ptr<SubDiscreteFunctionType>> b_h(number_of_interior_coarse_nodes_in_subgrid);
  std::vector<std::unique_ptr<SubDiscreteFunctionType>> rhs_Chj(number_of_interior_coarse_nodes_in_subgrid);
  for (int j = 0; j < number_of_interior_coarse_nodes_in_subgrid ; ++j)
  {
      b_h[j] = DSC::make_unique<SubDiscreteFunctionType>("q_h", subDiscreteFunctionSpace);
      rhs_Chj[j] = DSC::make_unique<SubDiscreteFunctionType>("rhs_Chj_h", subDiscreteFunctionSpace);
      b_h[j]->clear();
      rhs_Chj[j]->clear();
  }

  for (int j = 0; j < number_of_interior_coarse_nodes_in_subgrid ; ++j)
  {
      // clement_weight_j ( \psi_i, \Psi_j ), where
      // clement_weight_j = (*inverse_of_L1_norm_coarse_basis_funcs_)[interior_basis_func_id]
      
      // get the global id of all interior coarse basis functions (subgrid id, local id) -> (global interior id) 
      int interior_basis_func_id = (*ids_basis_functions_in_subgrid_)[coarse_index][j];
      local_problem_op.assemble_local_RHS_lg_problems( *((*coarse_basis_)[interior_basis_func_id]),
                                                       (*inverse_of_L1_norm_coarse_basis_funcs_)[interior_basis_func_id], *(rhs_Chj[j]) );
  }
DSC_LOG_INFO << "Time for assembling the the right hand sides rhs_Cj: " << assembleTimer.elapsed() << "s" << std::endl << std::endl; assembleTimer.reset();
  // zero boundary condition for 'rhs_Chj[j]':
  // set Dirichlet Boundary to zero
  for (SubgridIteratorType sg_it = subDiscreteFunctionSpace.begin(); sg_it != sg_end; ++sg_it)
  {
    const SubgridEntityType& subgrid_entity = *sg_it;

    HostEntityPointerType host_entity_pointer = subGrid.getHostEntity< 0 >(subgrid_entity);
    const HostEntityType& host_entity = *host_entity_pointer;

    HostIntersectionIterator iit = hostGridPart.ibegin(host_entity);
    const HostIntersectionIterator endiit = hostGridPart.iend(host_entity);
    for ( ; iit != endiit; ++iit)
    {
      if ( iit->neighbor() ) // if there is a neighbor entity
      {
        // check if the neighbor entity is in the subgrid
        const HostEntityPointerType neighborHostEntityPointer = iit->outside();
        const HostEntityType& neighborHostEntity = *neighborHostEntityPointer;

        if ( subGrid.contains< 0 >(neighborHostEntity) )
        {
          continue;
        }
      }

      //SubLocalFunctionType rhsLocal = rhs_Chj[j]->localFunction(subgrid_entity);
      const SGLagrangePointSetType& lagrangePointSet
          = subDiscreteFunctionSpace.lagrangePointSet(subgrid_entity);

      const int face = (*iit).indexInInside();

      FaceDofIteratorType faceIterator
          = lagrangePointSet.beginSubEntity< faceCodim >(face);
      const FaceDofIteratorType faceEndIterator
          = lagrangePointSet.endSubEntity< faceCodim >(face);

      for ( ; faceIterator != faceEndIterator; ++faceIterator)
        for (int j = 0; j < number_of_interior_coarse_nodes_in_subgrid ; ++j)
          ((rhs_Chj[j])->localFunction(subgrid_entity))[*faceIterator] = 0;
    }
  }
DSC_LOG_INFO << "Time for second boundary treatment C_h_j: " << assembleTimer.elapsed() << "s" << std::endl << std::endl; assembleTimer.reset();  
  // solve the pre-processing problems:
  for (int j = 0; j < number_of_interior_coarse_nodes_in_subgrid ; ++j)
     locprob_inverse_system_matrix( *(rhs_Chj[j]) , *(b_h[j]) );
DSC_LOG_INFO << "Time for solving all the local problems for b_h_j " << assembleTimer.elapsed() << "s" << std::endl << std::endl; assembleTimer.reset();
  // ----------------------------------------------------------------------------------------------------

  //! Solve the local problems without constraint
  // (without condition that the Lagrange interpolation must be zero)
  // ----------------------------------------------------------------------------------------------------

  // right hand side vectors of the algebraic local MsFEM problem
  SubDiscreteFunctionType local_problem_rhs_0("rhs of local MsFEM problem", subDiscreteFunctionSpace); // for e_0
  SubDiscreteFunctionType local_problem_rhs_1("rhs of local MsFEM problem", subDiscreteFunctionSpace); // for e_1
  
  local_problem_rhs_0.clear();
  local_problem_rhs_1.clear();

  // consider to make separate implementation of 'assemble_local_RHS' for the LOD method
  local_problem_op.assemble_local_RHS( e_0,
                                       subgrid_list_.getCoarseNodeVector( coarse_index ), /*coarse node vector is a dummy in this case*/
                                       specifier_.getOversamplingStrategy(), /*always '3' in this case */
                                       local_problem_rhs_0 );
  
  local_problem_op.assemble_local_RHS( e_1,
                                       subgrid_list_.getCoarseNodeVector( coarse_index ), /*coarse node vector is a dummy in this case*/
                                       specifier_.getOversamplingStrategy(), /*always three in this case*/
                                       local_problem_rhs_1 );
DSC_LOG_INFO << "Time for assembling the right hand sides: " << assembleTimer.elapsed() << "s" << std::endl << std::endl; assembleTimer.reset();
  local_problem_op.set_zero_boundary_condition_RHS( hostDiscreteFunctionSpace_ , local_problem_rhs_0 );
  local_problem_op.set_zero_boundary_condition_RHS( hostDiscreteFunctionSpace_ , local_problem_rhs_1 );
  
  //oneLinePrint( DSC_LOG_DEBUG, local_problem_rhs_0 );
  //oneLinePrint( DSC_LOG_DEBUG, local_problem_rhs_1 );
 
  locprob_inverse_system_matrix( local_problem_rhs_0 , local_problem_solution_0 );
  locprob_inverse_system_matrix( local_problem_rhs_1 , local_problem_solution_1 );
DSC_LOG_INFO << "Time for solving the problems for e_0 and e_1: " << assembleTimer.elapsed() << "s" << std::endl << std::endl; assembleTimer.reset();
  // ----------------------------------------------------------------------------------------------------
  
  
  //! Assemble and solve the problem for the lagrange multiplier v_h
  // (this is a low dimensional problem of size 'number_of_interior_coarse_nodes_in_subgrid')
  // ----------------------------------------------------------------------------------------------------
  
  // (stiffness) matrix for the lagrange multiplier (lm) problem  
  MatrixType lm_system_matrix( number_of_interior_coarse_nodes_in_subgrid, number_of_interior_coarse_nodes_in_subgrid );
  for (size_t i = 0; i != lm_system_matrix.N(); ++i) //rows
    for (size_t j = 0; j != lm_system_matrix.M(); ++j) //colums
      lm_system_matrix[i][j] = 0.0;

  // matrix with entries M[i][j] = weight_i ( b_h[j], coarse_basis_func[i] )_L2(\Omega)
  // 'i = row' and 'j = column'

  for (SubgridIteratorType sg_it = subDiscreteFunctionSpace.begin(); sg_it != sg_end; ++sg_it)
  {
    const SubgridEntityType& subgrid_entity = *sg_it;
    
    const SubGridEntityGeometry& sg_geometry = subgrid_entity.geometry();
    
    HostEntityPointerType host_entity_pointer = subGrid.getHostEntity< 0 >(subgrid_entity);
    const HostEntityType& host_entity = *host_entity_pointer;

    typedef CachingQuadrature< SubGridPartType, 0 > SubGridQuadrature;
    typedef CachingQuadrature< HostGridPartType, 0 > HostGridQuadrature;
    
    // exact for polynomials of degree 2:
    const SubGridQuadrature sg_quadrature( subgrid_entity, 2 * subDiscreteFunctionSpace.order() + 2);
    const HostGridQuadrature quadrature( host_entity, 2 * hostDiscreteFunctionSpace_.order() + 2);

    const HostGridEntityGeometry& geometry = host_entity.geometry();
    
    const size_t numQuadraturePoints = sg_quadrature.nop();
    for (size_t quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint)
    {
      const typename SubGridQuadrature::CoordinateType& local_point = sg_quadrature.point(quadraturePoint);

      const double quad_weight = sg_quadrature.weight(quadraturePoint)
              * sg_geometry.integrationElement(local_point);

      // global point in the subgrid
      const DomainType global_point = sg_geometry.global(local_point);
      
      // check for subgrid-hostgrid-compatibility
      assert( global_point == geometry.global( quadrature.point(quadraturePoint) ) );

      std::vector< RangeType > value_b, value_coarse_basis_func, clement_weight;
      value_b.resize(lm_system_matrix.M());
      value_coarse_basis_func.resize(lm_system_matrix.N());
      clement_weight.resize(lm_system_matrix.N());
      for (size_t j = 0; j != lm_system_matrix.M(); ++j) //rows
        ((b_h[j])->localFunction(subgrid_entity)).evaluate( sg_quadrature[quadraturePoint] , value_b[j]); 
      for (size_t i = 0; i != lm_system_matrix.N(); ++i) //rows
      {
        ((*coarse_basis_)[(*ids_basis_functions_in_subgrid_)[coarse_index][i]]
          ->localFunction(host_entity)).evaluate( quadrature[quadraturePoint] , value_coarse_basis_func[i]);
        clement_weight[i] = (*inverse_of_L1_norm_coarse_basis_funcs_)[(*ids_basis_functions_in_subgrid_)[coarse_index][i]];
      }
      
      for (size_t i = 0; i != lm_system_matrix.N(); ++i) //rows
        for (size_t j = 0; j != lm_system_matrix.M(); ++j) //colums
          lm_system_matrix[i][j] += quad_weight * clement_weight[i] * value_b[j] * value_coarse_basis_func[i];
    }
  } // lagrange multplier problem system matrix assembled
DSC_LOG_INFO << "Time for assembling the lm_system_matrix: " << assembleTimer.elapsed() << "s" << std::endl << std::endl; assembleTimer.reset();
  // print_matrix( lm_system_matrix );
  
  // right hand side vectors for the lagrange multiplier (lm) problems (for e_0 and e_1)
  // entries lm_rhs_0[i] = weight_i ( local_problem_solution_0, coarse_basis_func[i] )_L2(\Omega)
  VectorType lm_rhs_0( number_of_interior_coarse_nodes_in_subgrid );
  VectorType lm_rhs_1( number_of_interior_coarse_nodes_in_subgrid );
  for (size_t i = 0; i != number_of_interior_coarse_nodes_in_subgrid; ++i) //columns
  { lm_rhs_0[i] = 0.0; lm_rhs_1[i] = 0.0; }

  for (SubgridIteratorType sg_it = subDiscreteFunctionSpace.begin(); sg_it != sg_end; ++sg_it)
  {
    const SubgridEntityType& subgrid_entity = *sg_it;
    
    const SubGridEntityGeometry& sg_geometry = subgrid_entity.geometry();
    
    HostEntityPointerType host_entity_pointer = subGrid.getHostEntity< 0 >(subgrid_entity);
    const HostEntityType& host_entity = *host_entity_pointer;

    typedef CachingQuadrature< SubGridPartType, 0 > SubGridQuadrature;
    typedef CachingQuadrature< HostGridPartType, 0 > HostGridQuadrature;
    
    // exact for polynomials of degree 2:
    const SubGridQuadrature sg_quadrature( subgrid_entity, 2 * subDiscreteFunctionSpace.order() + 2);
    const HostGridQuadrature quadrature( host_entity, 2 * hostDiscreteFunctionSpace_.order() + 2);
    
    RangeType value_local_problem_solution_0, value_local_problem_solution_1, value_coarse_basis_func_i;

    const size_t numQuadraturePoints = sg_quadrature.nop();
    for (size_t quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint)
    {
      const typename SubGridQuadrature::CoordinateType& local_point = sg_quadrature.point(quadraturePoint);

      const double quad_weight = sg_quadrature.weight(quadraturePoint)
              * sg_geometry.integrationElement(local_point);
      
      // check for subgrid-hostgrid-compatibility
      assert( sg_geometry.global(local_point) == host_entity.geometry().global( quadrature.point(quadraturePoint) ) );

      const SubLocalFunctionType local_sol_0 = local_problem_solution_0.localFunction(subgrid_entity);
      const SubLocalFunctionType local_sol_1 = local_problem_solution_1.localFunction(subgrid_entity);

      local_sol_0.evaluate( sg_quadrature[quadraturePoint] , value_local_problem_solution_0);
      local_sol_1.evaluate( sg_quadrature[quadraturePoint] , value_local_problem_solution_1);

      for (size_t i = 0; i != number_of_interior_coarse_nodes_in_subgrid; ++i) //columns
      {

          int interior_coarse_basis_id_in_subgrid = (*ids_basis_functions_in_subgrid_)[coarse_index][i];
          HostLocalFunctionType local_coarse_basis_i
                 = (*coarse_basis_)[interior_coarse_basis_id_in_subgrid]->localFunction(host_entity);
          local_coarse_basis_i.evaluate( quadrature[quadraturePoint] , value_coarse_basis_func_i);

          double clement_weight = (*inverse_of_L1_norm_coarse_basis_funcs_)[interior_coarse_basis_id_in_subgrid];

          lm_rhs_0[i] += quad_weight * clement_weight * value_local_problem_solution_0 * value_coarse_basis_func_i;
          lm_rhs_1[i] += quad_weight * clement_weight * value_local_problem_solution_1 * value_coarse_basis_func_i;

      }
    }
  } // lagrange multplier problem system matrix assembled
DSC_LOG_INFO << "Time for assembling the lm_right_hand sides: " << assembleTimer.elapsed() << "s" << std::endl << std::endl; assembleTimer.reset();  
  //print_vector( lm_rhs_0 );
  //print_vector( lm_rhs_1 );
  
  MatrixOperatorType lm_matrix_op_0( lm_system_matrix );
  MatrixOperatorType lm_matrix_op_1( lm_system_matrix );
  
  PreconditionerType lm_preconditioner_0( lm_system_matrix, 100, 0.9 );
  PreconditionerType lm_preconditioner_1( lm_system_matrix, 100, 0.9 );
  
  Dune::InverseOperatorResult result_data_0, result_data_1;
  VectorType v_h_0( number_of_interior_coarse_nodes_in_subgrid );
  VectorType v_h_1( number_of_interior_coarse_nodes_in_subgrid );
  for (size_t col = 0; col != v_h_0.N(); ++col)
  {  v_h_0[col] = 0.0; v_h_1[col] = 0.0;  }
  
  typedef Dune::BiCGSTABSolver< VectorType > SolverType;

  double tol = DSC_CONFIG_GET("rigorous_msfem.local_micro_solver_tolerance", 1e-10 );
  int num_iterations = DSC_CONFIG_GET("rigorous_msfem.local_micro_solver_iterations", 10000 );
  
  SolverType lm_prob_solver_0( lm_matrix_op_0, lm_preconditioner_0, tol, num_iterations, false );
  SolverType lm_prob_solver_1( lm_matrix_op_1, lm_preconditioner_1, tol, num_iterations, false );
  
  lm_prob_solver_0.apply( v_h_0, lm_rhs_0, result_data_0);
  lm_prob_solver_1.apply( v_h_1, lm_rhs_1, result_data_1);
DSC_LOG_INFO << "Time for solving the lm problem: " << assembleTimer.elapsed() << "s" << std::endl << std::endl; assembleTimer.reset();    
  // set final_rhs_0 = \sum_j clement_weight_j v_h_0[j] coarse_basis_function[j]
  HostDiscreteFunctionType final_rhs_0("final rhs 0 local problem", hostDiscreteFunctionSpace_); // for e_0
  HostDiscreteFunctionType final_rhs_1("final rhs 1 local problem", hostDiscreteFunctionSpace_); // for e_1
  final_rhs_0.clear();
  final_rhs_1.clear();

  for (size_t i = 0; i != number_of_interior_coarse_nodes_in_subgrid; ++i) //columns
  {

     int interior_coarse_basis_id_in_subgrid = (*ids_basis_functions_in_subgrid_)[coarse_index][i];

     HostDiscreteFunctionType aux_func_0("auxilliary func 0", hostDiscreteFunctionSpace_);
     HostDiscreteFunctionType aux_func_1("auxilliary func 1", hostDiscreteFunctionSpace_);
     aux_func_0.clear();
     aux_func_1.clear();

     double coefficient_0 = (*inverse_of_L1_norm_coarse_basis_funcs_)[interior_coarse_basis_id_in_subgrid] * v_h_0[i];
     double coefficient_1 = (*inverse_of_L1_norm_coarse_basis_funcs_)[interior_coarse_basis_id_in_subgrid] * v_h_1[i];
     
     aux_func_0 += (*(*coarse_basis_)[interior_coarse_basis_id_in_subgrid]);
     aux_func_0 *= coefficient_0;
     final_rhs_0 += aux_func_0;

     aux_func_1 += (*(*coarse_basis_)[interior_coarse_basis_id_in_subgrid]);
     aux_func_1 *= coefficient_1;
     final_rhs_1 += aux_func_1;
  }
DSC_LOG_INFO << "Time for computing/assembling v_h : " << assembleTimer.elapsed() << "s" << std::endl << std::endl; assembleTimer.reset();  
  // right hand side vectors of the algebraic local MsFEM problem
  SubDiscreteFunctionType final_rhs_vector_0("final rhs of local MsFEM problem", subDiscreteFunctionSpace); // for e_0
  SubDiscreteFunctionType final_rhs_vector_1("final rhs of local MsFEM problem", subDiscreteFunctionSpace); // for e_1
  final_rhs_vector_0.clear();
  final_rhs_vector_1.clear();
  
  local_problem_op.assemble_local_RHS_lg_problems( final_rhs_0, 1.0, final_rhs_vector_0 );
  local_problem_op.assemble_local_RHS_lg_problems( final_rhs_1, 1.0, final_rhs_vector_1 );
  local_problem_op.set_zero_boundary_condition_RHS( hostDiscreteFunctionSpace_ , final_rhs_vector_0 );
  local_problem_op.set_zero_boundary_condition_RHS( hostDiscreteFunctionSpace_ , final_rhs_vector_1 );
  
  SubDiscreteFunctionType preliminary_solution_0("preliminary_solution_0", subDiscreteFunctionSpace); // for e_0
  SubDiscreteFunctionType preliminary_solution_1("preliminary_solution_1", subDiscreteFunctionSpace); // for e_1
  preliminary_solution_0.clear();
  preliminary_solution_1.clear();

  locprob_inverse_system_matrix( final_rhs_vector_0 , preliminary_solution_0 );
  locprob_inverse_system_matrix( final_rhs_vector_1 , preliminary_solution_1 );
DSC_LOG_INFO << "Time for solving the final problems: " << assembleTimer.elapsed() << "s" << std::endl << std::endl; assembleTimer.reset();    
  local_problem_solution_0 -= preliminary_solution_0;
  local_problem_solution_1 -= preliminary_solution_1;
DSC_LOG_INFO << "Time for finanlizing step: " << assembleTimer.elapsed() << "s" << std::endl << std::endl; assembleTimer.reset();  
#if 0


  // zero boundary condition for 'cell problems':
  // set Dirichlet Boundary to zero
  for (SubgridIteratorType sg_it = subDiscreteFunctionSpace.begin(); sg_it != sg_end; ++sg_it)
  {
    const SubgridEntityType& subgrid_entity = *sg_it;

    HostEntityPointerType host_entity_pointer = subGrid.getHostEntity< 0 >(subgrid_entity);
    const HostEntityType& host_entity = *host_entity_pointer;

    HostIntersectionIterator iit = hostGridPart.ibegin(host_entity);
    const HostIntersectionIterator endiit = hostGridPart.iend(host_entity);
    for ( ; iit != endiit; ++iit)
    {
      if ( iit->neighbor() ) // if there is a neighbor entity
      {
        // check if the neighbor entity is in the subgrid
        const HostEntityPointerType neighborHostEntityPointer = iit->outside();
        const HostEntityType& neighborHostEntity = *neighborHostEntityPointer;

        if ( subGrid.contains< 0 >(neighborHostEntity) )
        {
          continue;
        }
      }

      SubLocalFunctionType rhsLocal = local_problem_rhs.localFunction(subgrid_entity);
      const SGLagrangePointSetType& lagrangePointSet
          = subDiscreteFunctionSpace.lagrangePointSet(subgrid_entity);

      const int face = (*iit).indexInInside();

      FaceDofIteratorType faceIterator
          = lagrangePointSet.beginSubEntity< faceCodim >(face);
      const FaceDofIteratorType faceEndIterator
          = lagrangePointSet.endSubEntity< faceCodim >(face);
      for ( ; faceIterator != faceEndIterator; ++faceIterator)
        rhsLocal[*faceIterator] = 0;
    }
  }

  const double norm_rhs = local_problem_op.normRHS(local_problem_rhs);

  if ( !( local_problem_rhs.dofsValid() ) )
  {
    DUNE_THROW(Dune::InvalidStateException, "Local MsFEM Problem RHS invalid.");
  }

  if (norm_rhs < /*1e-06*/ 1e-30)
  {
    local_problem_solution.clear();
    DSC_LOG_ERROR << "Local MsFEM problem with solution zero." << std::endl;
  } else {

    InverseLocProbFEMMatrix locprob_fem_biCGStab(locprob_system_matrix, 1e-8, 1e-8, 20000, DSC_CONFIG_GET("localproblemsolver_verbose", false));
  
  

  
    bool clement = false;
    if ( specifier_.getOversamplingStrategy() == 3 )
    { clement = (DSC_CONFIG_GET( "rigorous_msfem.oversampling_strategy", "Clement" ) == "Clement" ); }

    if ( clement )
    {
           
//! UNDER CONSTRUCTION!
//! new implementation: direct inversion strategy.
#if 1
    // we solve the saddle point problem with a direct inversion strategy for the Schur complement
    
    // Let S_h(=locprob_system_matrix) denote the classical stiffness/system matrix corresponding to the fine basis functions in the subgrid 
    // r_h (=local_problem_rhs) is  the standard right hand side of the local problems ( - \int_T A e_k \nabla \phi_j )
    // First, we solve for q_h with S_h q_h = r_h
    SubDiscreteFunctionType q_h("q_h", subDiscreteFunctionSpace);
    q_h.clear();
    locprob_fem_biCGStab(local_problem_rhs, q_h );
    
    int number_of_interior_coarse_nodes_in_subgrid = ids_basis_functions_in_subgrid_[ coarse_index ].size();
    
    // For each coarse node j in the subgrid (local internal numbering, i.e. 0 <= j < M_subgrid), solve
    // for b_h_j with S_h b_h_j = C_h^T e_j, where C_h describes the algebraic version of the weighted
    // Clement interpolation operator
    SubDiscreteFunctionType** b_h = new SubDiscreteFunctionType* [number_of_interior_coarse_nodes_in_subgrid];
    SubDiscreteFunctionType** rhs_Chj = new SubDiscreteFunctionType* [number_of_interior_coarse_nodes_in_subgrid];
    for (int j = 0; j < number_of_interior_coarse_nodes_in_subgrid ; ++j)
    {
      b_h[j] = new SubDiscreteFunctionType("q_h", subDiscreteFunctionSpace);
      rhs_Chj[j] = new SubDiscreteFunctionType("q_h", subDiscreteFunctionSpace);
      b_h[j]->clear();
      rhs_Chj[j]->clear();
    }


#if 0
    std::vector< SubDiscreteFunctionType* > b_h;
    for (int j = 0; j < number_of_interior_coarse_nodes_in_subgrid ; ++j)
    { 
      SubDiscreteFunctionType b_h_j("q_h", subDiscreteFunctionSpace);
      b_h_j.clear();
      b_h.push_back( SubDiscreteFunctionType("q_h", subDiscreteFunctionSpace) );
    }
#endif
    
    // solve the local problems without using the dune-fem structure.
    // At the end: copy the solution vector into a dune-fem discrete function by using
     //const int global_dof_number = space.mapper().mapToGlobal(*it, loc_basis_number );

    for (int j = 0; j < number_of_interior_coarse_nodes_in_subgrid ; ++j)
    {  delete[] b_h[j]; delete[] rhs_Chj[j]; }
    delete[] b_h;
    delete[] rhs_Chj;
#endif

//! old implementation using uzawa solver:
#if 1
      HostDiscreteFunctionType zero("zero", specifier_.coarseSpace());
      zero.clear();
      const double dummy = 12345.67890;
      double solverEps = 1e-2;
      int maxIterations = 1000;

      // we want to solve the local problem with the constraint that the weighted Clement interpoltion
      // of the local problem solution is zero

      // implementation of a weighted Clement interpolation operator for our purpose:
      WeightedClementOperatorType clement_interpolation_op( subDiscreteFunctionSpace,
                                                            specifier_.coarseSpace(),
                                                            subgrid_list_.getCoarseNodeVector( coarse_index ),
                                                            *coarse_basis_, *global_id_to_internal_id_, specifier_ );
      //! NOTE TODO: implementation is not yet optimal, because the weighted Clement maps a function
      //! defined on the local subgrid to a function defined on the whole(!) coarse space.
      //! It would be better to implement a mapping to a localized coarse space, since
      //! the uzawa solver must treat ALL coarse grid nodes (expensive and worse convergence).

      //clement_interpolation_op.print();

      HostDiscreteFunctionType lagrange_multiplier("lagrange multiplier", specifier_.coarseSpace() );
      lagrange_multiplier.clear();

      // create inverse operator
      // saddle point problem solver with uzawa algorithm:
      {
        DSC::Profiler::ScopedTiming st("uzawa");
        InverseUzawaOperatorType uzawa( locprob_fem_biCGStab, clement_interpolation_op, dummy, solverEps, maxIterations, true);
        uzawa( local_problem_rhs, zero /*interpolation is zero*/, local_problem_solution, lagrange_multiplier );
      }
#endif
    }
    else {
      locprob_fem_biCGStab(local_problem_rhs, local_problem_solution);
    }
  }

  if ( !( local_problem_solution.dofsValid() ) ) {
    DUNE_THROW(Dune::InvalidStateException,"Current solution of the local msfem problem invalid!");
  }

  // oneLinePrint( DSC_LOG_DEBUG, local_problem_solution );
#endif

} // solvelocalproblem

void MsFEMLocalProblemSolver::subgrid_to_hostrid_function(const SubDiscreteFunctionType& sub_func,
                                 HostDiscreteFunctionType& host_func) {
  host_func.clear();

  const SubDiscreteFunctionSpaceType& subDiscreteFunctionSpace = sub_func.space();
  const SubGridType& subGrid = subDiscreteFunctionSpace.grid();

  SubgridIteratorType sub_endit = subDiscreteFunctionSpace.end();
  for (SubgridIteratorType sub_it = subDiscreteFunctionSpace.begin(); sub_it != sub_endit; ++sub_it)
  {
    const SubgridEntityType& sub_entity = *sub_it;

    HostEntityPointerType host_entity_pointer = subGrid.getHostEntity< 0 >(*sub_it);
    const HostEntityType& host_entity = *host_entity_pointer;

    SubLocalFunctionType sub_loc_value = sub_func.localFunction(sub_entity);
    HostLocalFunctionType host_loc_value = host_func.localFunction(host_entity);

    const auto numBaseFunctions = sub_loc_value.baseFunctionSet().size();
    for (unsigned int i = 0; i < numBaseFunctions; ++i)
    {
      host_loc_value[i] = sub_loc_value[i];
    }
  }
} // subgrid_to_hostrid_function

void MsFEMLocalProblemSolver::output_local_solution(const int coarse_index, const int which,
                           const HostDiscreteFunctionType& host_local_solution) const
{
  if (!DSC_CONFIG_GET("global.local_solution_vtk_output", false))
    return;
  typedef tuple< const HostDiscreteFunctionType* >      IOTupleType;
  typedef DataOutput< HostGridType, IOTupleType > DataOutputType;

  // general output parameters
  LocalProblemDataOutputParameters outputparam;
  // --------- data output local solution --------------

  // create and initialize output class
  IOTupleType local_solution_series(&host_local_solution);

  const std::string ls_name_s = (boost::format("/local_problem_solution_e%d_%d") % which % coarse_index).str();

  outputparam.set_prefix(ls_name_s);
  DataOutputType localsol_dataoutput(
    hostDiscreteFunctionSpace_.gridPart().grid(), local_solution_series, outputparam);
  localsol_dataoutput.writeData( 1.0 /*dummy*/, (boost::format("local-problem-solution-%d") % which).str() );
}

void MsFEMLocalProblemSolver::assemble_all(bool /*silent*/) {
  std::string local_path = DSC_CONFIG_GET("global.datadir", "data") + "/local_problems/";
  Dune::Stuff::Common::testCreateDirectory(local_path);

  enum { dimension = CommonTraits::GridType::dimension };
  enum { maxnumOfBaseFct = 100 };

  JacobianRangeType e[dimension];
  for (int i = 0; i < dimension; ++i)
    for (int j = 0; j < dimension; ++j) {
      if (i == j) {
        e[i][0][j] = 1.0;
      } else {
        e[i][0][j] = 0.0;
      }
    }

  // number of coarse grid entities (of codim 0).
  int number_of_coarse_grid_entities = specifier_.getNumOfCoarseEntities();

  DSC_LOG_INFO << "in method 'assemble_all': number_of_coarse_grid_entities = " << number_of_coarse_grid_entities
            << std::endl;
  DSC_PROFILER.startTiming("msfem.localproblemsolver.assemble_all");

  // we want to determine minimum, average and maxiumum time for solving a local msfem problem in the current method
  Dune::Stuff::Common::MinMaxAvg<double> cell_time;

  std::vector<int> coarse_indices;
  const HostDiscreteFunctionSpaceType& coarseSpace = specifier_.coarseSpace();
  const HostGridLeafIndexSet& coarseGridLeafIndexSet = coarseSpace.gridPart().grid().leafIndexSet();
  for (HostGridEntityIteratorType coarse_it = coarseSpace.begin(); coarse_it != coarseSpace.end(); ++coarse_it)
  {
    coarse_indices.push_back(coarseGridLeafIndexSet.index(*coarse_it));
  }


//  const auto& comm = Dune::MPIHelper::getCollectiveCommunication();
//  int slice = coarse_indices.size() / comm.size();
//  for(int gc = comm.rank() * slice; gc < std::min(long(comm.rank() +1)* slice, long(coarse_indices.size())); ++gc)
  for (int gc=0; gc<coarse_indices.size(); ++gc)
  {
    const int coarse_index = coarse_indices[gc];

    DSC_LOG_INFO << "-------------------------" << std::endl
                 << "Coarse index " << coarse_index << std::endl;

    const std::string locprob_solution_location =
        (boost::format("local_problems/_localProblemSolutions_%d_%d") % coarse_index % MPIManager::rank()).str();

    DiscreteFunctionWriter dfw(locprob_solution_location);

    SubGridType& subGrid = subgrid_list_.getSubGrid(coarse_index);

    SubGridPartType subGridPart(subGrid);

    const SubDiscreteFunctionSpaceType subDiscreteFunctionSpace(subGridPart);

    const std::string name_local_solution = (boost::format("Local Problem Solution %d") % coarse_index).str();

    //! only for dimension 2!
    SubDiscreteFunctionType local_problem_solution_0(name_local_solution, subDiscreteFunctionSpace);
    local_problem_solution_0.clear();

    SubDiscreteFunctionType local_problem_solution_1(name_local_solution, subDiscreteFunctionSpace);
    local_problem_solution_1.clear();

    Dune::Timer assembleTimer;
    bool uzawa = DSC_CONFIG_GET( "rigorous_msfem.uzawa_solver", false );
    bool clement = ( DSC_CONFIG_GET( "rigorous_msfem.oversampling_strategy", "Clement" ) == "Clement" );
    if ( (!uzawa) && (specifier_.getOversamplingStrategy() == 3) && clement ) {

      // requires a pre-processing step (that is the same for both directions e_0 and e_1)
      // one method for both solutions to half the computational complexity 
      solvelocalproblems_lod( e[0], e[1], local_problem_solution_0,
                              local_problem_solution_1, coarse_index ); 

    }
    else {

      DSC_LOG_INFO  << std::endl
                    << "Number of the local problem: " << dimension * coarse_index << " (of "
                    << (dimension * number_of_coarse_grid_entities) - 1 << " problems in total)" << std::endl
                    << "   Subgrid " << coarse_index << " contains " << subGrid.size(0) << " elements and "
                    << subGrid.size(2) << " nodes." << std::endl;

      // take time
      DSC_PROFILER.startTiming("none.local_problem_solution");

      // solve the problems
      solvelocalproblem(e[0], local_problem_solution_0, coarse_index);

      cell_time(DSC_PROFILER.stopTiming("none.local_problem_solution") / 1000.f);
      DSC_PROFILER.resetTiming("none.local_problem_solution");

      DSC_LOG_INFO  << std::endl
                    << "Number of the local problem: "
                    << (dimension * coarse_index) + 1 << " (of "
                    << (dimension * number_of_coarse_grid_entities) - 1 << " problems in total)" << std::endl
                    << "   Subgrid " << coarse_index << " contains " << subGrid.size(0) << " elements and "
                    << subGrid.size(2) << " nodes." << std::endl;

      // take time
      DSC_PROFILER.startTiming("none.local_problem_solution");

      // solve the problems
      solvelocalproblem(e[1], local_problem_solution_1, coarse_index);

      // min/max time
      cell_time(DSC_PROFILER.stopTiming("none.local_problem_solution") / 1000.f);
      DSC_PROFILER.resetTiming("none.local_problem_solution");
    }

    DSC_LOG_INFO << "Total time for solving all local problems for the current subgrid: " << assembleTimer.elapsed() << "s" << std::endl << std::endl;
//abort();
    dfw.append(local_problem_solution_0);
    dfw.append(local_problem_solution_1);

    HostDiscreteFunctionType host_local_solution(name_local_solution, hostDiscreteFunctionSpace_);
    subgrid_to_hostrid_function(local_problem_solution_0, host_local_solution);
    output_local_solution(coarse_index, 0, host_local_solution);
      
    subgrid_to_hostrid_function(local_problem_solution_1, host_local_solution);
    output_local_solution(coarse_index, 1, host_local_solution);
  } //for

  const auto total_time = DSC_PROFILER.stopTiming("msfem.localproblemsolver.assemble_all")/1000.f;
  DSC_LOG_INFO << std::endl;
  DSC_LOG_INFO << "In method: assemble_all." << std::endl << std::endl;
  DSC_LOG_INFO << "MsFEM problems solved for " << number_of_coarse_grid_entities << " coarse grid entities."
                << std::endl;
  DSC_LOG_INFO << dimension * number_of_coarse_grid_entities << " local MsFEM problems solved in total."
                << std::endl;
  DSC_LOG_INFO << "Minimum time for solving a local problem = " << cell_time.min() << "s." << std::endl;
  DSC_LOG_INFO << "Maximum time for solving a localproblem = " << cell_time.max() << "s." << std::endl;
  DSC_LOG_INFO << "Average time for solving a localproblem = "
                << cell_time.average()<< "s." << std::endl;
  DSC_LOG_INFO << "Total time for computing and saving the localproblems = "
                  << total_time << "s," << std::endl << std::endl;
} // assemble_all

} //namespace MsFEM {
} //namespace Multiscale {
} //namespace Dune {
