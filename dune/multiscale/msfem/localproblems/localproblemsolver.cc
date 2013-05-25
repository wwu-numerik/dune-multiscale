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

#include <dune/multiscale/tools/subgrid_io.hh>
#include <dune/multiscale/hmm/cell_problem_numbering.hh>
#include <dune/multiscale/msfem/msfem_traits.hh>
#include <dune/multiscale/msfem/localproblems/localoperator.hh>

#include <dune/multiscale/tools/misc/uzawa.hh>
#include <dune/multiscale/tools/misc/weighted-clement-operator.hh>
#include <dune/multiscale/tools/discretefunctionwriter.hh>

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
  case 1: break;
  case 2: break;
  case 3: break;
  default: DUNE_THROW(Dune::InvalidStateException, "Oversampling Strategy must be 1 or 2.");
  }

  // assemble the stiffness matrix
  if ( specifier_.getOversamplingStrategy() == 1 ) {
    local_problem_op.assemble_matrix(locprob_system_matrix);
  }

  if ( specifier_.getOversamplingStrategy() == 2 ) {
    if ( coarse_index < 0 )
      DUNE_THROW(Dune::InvalidStateException, "Invalid coarse index: coarse_index < 0");
    local_problem_op.assemble_matrix(locprob_system_matrix, subgrid_list_.getCoarseNodeVector( coarse_index ) );
  }

  if ( specifier_.getOversamplingStrategy() == 3 ) {
    if ( coarse_index < 0 )
      DUNE_THROW(Dune::InvalidStateException, "Invalid coarse index: coarse_index < 0");
    bool clement = ( DSC_CONFIG_GET( "rigorous_msfem.oversampling_strategy", "Clement" ) == "Clement" );

    if ( clement ) {
      local_problem_op.assemble_matrix( locprob_system_matrix );
    } else {
      local_problem_op.assemble_matrix( locprob_system_matrix, subgrid_list_.getCoarseNodeVector( coarse_index ) );
    }
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


  // assemble right hand side of algebraic local msfem problem
  if ( specifier_.getOversamplingStrategy() == 1 ) {
    local_problem_op.assemble_local_RHS(e, local_problem_rhs);
  } else
  if ( ( specifier_.getOversamplingStrategy() == 2 ) || ( specifier_.getOversamplingStrategy() == 3 ) ) {
    if ( coarse_index < 0 )
      DUNE_THROW(Dune::InvalidStateException, "Invalid coarse index: coarse_index < 0");
    local_problem_op.assemble_local_RHS(e,
            subgrid_list_.getCoarseNodeVector( coarse_index ),
            specifier_.getOversamplingStrategy(),
            local_problem_rhs );
  } else
    DUNE_THROW(Dune::InvalidStateException, "Oversampling Strategy must be 1, 2 or 3!");
  //oneLinePrint( DSC_LOG_DEBUG, local_problem_rhs );

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

      //! old implementation LODM with Clement using uzawa solver:
      // ------------------------------------------------------------
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
      // ------------------------------------------------------------
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

  // assemble stiffness matrix
  local_problem_op.assemble_matrix( locprob_system_matrix );


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

  InverseLocProbFEMMatrix locprob_inverse_system_matrix(locprob_system_matrix,
                                                        1e-8, 1e-8, 20000,
                                                        DSC_CONFIG_GET("localproblemsolver_verbose", false));

  //! Pre-processing step
  // For each coarse node j in the subgrid (local internal numbering, i.e. 0 <= j < M_subgrid), solve
  // for b_h_j with S_h b_h_j = C_h^T e_j, where C_h describes the algebraic version of the weighted
  // Clement interpolation operator
  // ----------------------------------------------------------------------------------------------------
  int number_of_interior_coarse_nodes_in_subgrid = (*ids_basis_functions_in_subgrid_)[ coarse_index ].size();
    
  SubDiscreteFunctionType** b_h = new SubDiscreteFunctionType* [number_of_interior_coarse_nodes_in_subgrid];
  SubDiscreteFunctionType** rhs_Chj = new SubDiscreteFunctionType* [number_of_interior_coarse_nodes_in_subgrid];
  for (int j = 0; j < number_of_interior_coarse_nodes_in_subgrid ; ++j)
  {
      b_h[j] = new SubDiscreteFunctionType("q_h", subDiscreteFunctionSpace);
      rhs_Chj[j] = new SubDiscreteFunctionType("rhs_Chj_h", subDiscreteFunctionSpace);
      b_h[j]->clear();
      rhs_Chj[j]->clear();
  }

  for (int j = 0; j < number_of_interior_coarse_nodes_in_subgrid ; ++j)
  {
      int interior_basis_func_id = (*ids_basis_functions_in_subgrid_)[coarse_index][j];
      local_problem_op.assemble_local_RHS_pre_processing( *((*coarse_basis_)[interior_basis_func_id]),
                                                          (*inverse_of_L1_norm_coarse_basis_funcs_)[interior_basis_func_id], *(rhs_Chj[j]) );
  }

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
  
  // solve the pre-processing problems:
  for (int j = 0; j < number_of_interior_coarse_nodes_in_subgrid ; ++j)
     locprob_inverse_system_matrix( *(rhs_Chj[j]) , *(b_h[j]) );

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

  local_problem_op.set_zero_boundary_condition_RHS( hostDiscreteFunctionSpace_ , local_problem_rhs_0 );
  local_problem_op.set_zero_boundary_condition_RHS( hostDiscreteFunctionSpace_ , local_problem_rhs_1 );
  
  //oneLinePrint( DSC_LOG_DEBUG, local_problem_rhs_0 );
  //oneLinePrint( DSC_LOG_DEBUG, local_problem_rhs_1 );
 
  locprob_inverse_system_matrix( local_problem_rhs_0 , local_problem_solution_0 );
  locprob_inverse_system_matrix( local_problem_rhs_1 , local_problem_solution_1 );
  
  // ----------------------------------------------------------------------------------------------------
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
  
  
// assemble_local_RHS_pre_processing(
  // std::vector< std::vector< int > > ids_basis_functions_in_subgrid,
  // std::vector< double >& inverse_of_L1_norm_coarse_basis_funcs,
  // CoarseBasisFunctionListType& coarse_basis
    
  
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
    for (int j = 0; j < number_of_interior_coarse_nodes_in_subgrid ; ++j)
    {  delete b_h[j]; delete rhs_Chj[j]; }
    delete[] b_h;
    delete[] rhs_Chj;

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

    bool clement = ( DSC_CONFIG_GET( "rigorous_msfem.oversampling_strategy", "Clement" ) == "Clement" );
    if ( (specifier_.getOversamplingStrategy() == 3) && clement ) {

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
