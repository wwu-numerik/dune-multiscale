#ifndef DiscreteEllipticMsFEMLocalProblem_HH
#define DiscreteEllipticMsFEMLocalProblem_HH

#include <vector>

#include <dune/common/fmatrix.hh>
#include <dune/fem/operator/matrix/spmatrix.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/gridpart/common/gridpart.hh>
#include <dune/fem/operator/2order/lagrangematrixsetup.hh>
#include <dune/subgrid/subgrid.hh>
#include <dune/stuff/common/filesystem.hh>
#include <dune/stuff/fem/functions/checks.hh>
#include <dune/multiscale/tools/subgrid_io.hh>
#include <dune/multiscale/tools/disc_func_writer/discretefunctionwriter.hh>
#include <dune/multiscale/tools/misc/outputparameter.hh>
#include <dune/multiscale/msfem/msfem_traits.hh>
#include <dune/multiscale/tools/solver/MsFEM/msfem_localproblems/localoperator.hh>
#include <dune/multiscale/tools/misc/uzawa.hh>
#include <dune/multiscale/tools/misc/weighted-clement-operator.hh>

namespace Dune {
namespace Multiscale {
namespace MsFEM {

/** \brief define output parameters for local problems
 *  appends "local_problems" for path
 **/
struct LocalProblemDataOutputParameters
  : public myDataOutputParameters
{
public:
  explicit LocalProblemDataOutputParameters()
    : myDataOutputParameters(DSC_CONFIG_GET("global.datadir", "data") + "/local_problems/")
  {}
};

//! --------------------- the essential local msfem problem solver class ---------------------------
template< class HostDiscreteFunctionType,
          class SubGridListType,
          class MacroMicroGridSpecifierType,
          class DiffusionOperatorType,
          class CoarseBasisFunctionListType >
class MsFEMLocalProblemSolver
{
private:
  //! type of discrete function space
  typedef typename HostDiscreteFunctionType::DiscreteFunctionSpaceType HostDiscreteFunctionSpaceType;
  //! type of (non-discrete )function space
  typedef typename HostDiscreteFunctionSpaceType::FunctionSpaceType FunctionSpaceType;
  typedef typename HostDiscreteFunctionSpaceType::GridPartType HostGridPartType;
  typedef typename HostDiscreteFunctionSpaceType::GridType HostGridType;
  typedef typename HostDiscreteFunctionSpaceType::RangeType RangeType;
  //! type of value of a gradient of a function
  typedef typename HostDiscreteFunctionSpaceType::JacobianRangeType JacobianRangeType;
  typedef typename HostDiscreteFunctionSpaceType::DomainType DomainType;

  typedef typename HostGridType::Traits::LeafIndexSet HostGridLeafIndexSet;
  typedef typename HostDiscreteFunctionSpaceType::IteratorType HostGridEntityIteratorType;
  typedef typename HostGridEntityIteratorType::Entity HostEntityType;
  typedef typename HostEntityType::EntityPointer HostEntityPointerType;
  typedef typename HostGridType::template Codim< 0 >::Geometry HostGridEntityGeometry;
  typedef typename HostDiscreteFunctionType::LocalFunctionType HostLocalFunctionType;
  typedef typename HostGridPartType::IntersectionIteratorType HostIntersectionIterator;

  //! ---------------- typedefs for the SubgridDiscreteFunctionSpace -----------------------
  // ( typedefs for the local grid and the corresponding local ('sub') )discrete space )

  //! type of grid
  typedef typename SubGridListType::SubGridType SubGridType;
  //! type of grid part
  typedef LeafGridPart< SubGridType > SubGridPartType;
  //! type of subgrid discrete function space
  typedef LagrangeDiscreteFunctionSpace< FunctionSpaceType, SubGridPartType, 1 >  // 1=POLORDER
    SubDiscreteFunctionSpaceType;

  //! type of subgrid discrete function
  typedef AdaptiveDiscreteFunction< SubDiscreteFunctionSpaceType > SubDiscreteFunctionType;
  typedef typename SubDiscreteFunctionSpaceType::IteratorType SubgridIteratorType;
  typedef typename SubgridIteratorType::Entity SubgridEntityType;
  typedef typename SubgridEntityType::EntityPointer SubgridEntityPointerType;
  typedef typename SubDiscreteFunctionType::LocalFunctionType SubLocalFunctionType;
  typedef typename SubDiscreteFunctionSpaceType::LagrangePointSetType SGLagrangePointSetType;
  typedef typename SubDiscreteFunctionSpaceType::LagrangePointSetType SubgridLagrangePointSetType;

  enum { faceCodim = 1 };
  typedef typename SubgridLagrangePointSetType::template Codim< faceCodim >::SubEntityIteratorType
    SubgridFaceDofIteratorType;

  //! polynomial order of base functions
  enum { polynomialOrder = SubDiscreteFunctionSpaceType::polynomialOrder };

  struct LocProbMatrixTraits
  {
    typedef SubDiscreteFunctionSpaceType                          RowSpaceType;
    typedef SubDiscreteFunctionSpaceType                          ColumnSpaceType;
    typedef LagrangeMatrixSetup< false >                          StencilType;
    typedef ParallelScalarProduct< SubDiscreteFunctionSpaceType > ParallelScalarProductType;

    template< class M >
    struct Adapter
    {
      typedef LagrangeParallelMatrixAdapter< M > MatrixAdapterType;
    };
  };

  typedef SparseRowMatrixOperator< SubDiscreteFunctionType, SubDiscreteFunctionType,
                                   LocProbMatrixTraits > LocProbFEMMatrix;

  #ifdef SYMMETRIC_DIFFUSION_MATRIX
  typedef CGInverseOperator< SubDiscreteFunctionType, LocProbFEMMatrix > InverseLocProbFEMMatrix;
  #else
  // OEMGMRESOp //OEMBICGSQOp // OEMBICGSTABOp
  typedef OEMBICGSTABOp< SubDiscreteFunctionType, LocProbFEMMatrix > InverseLocProbFEMMatrix;
  #endif // ifdef SYMMETRIC_DIFFUSION_MATRIX

  // discrete elliptic operator describing the elliptic local msfem problems
  typedef LocalProblemOperator< SubDiscreteFunctionType, DiffusionOperatorType > LocalProblemOperatorType;

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
     
  const HostDiscreteFunctionSpaceType& hostDiscreteFunctionSpace_;
  const DiffusionOperatorType& diffusion_;
  const MacroMicroGridSpecifierType& specifier_;
  SubGridListType& subgrid_list_;

  const CoarseBasisFunctionListType* coarse_basis_;
  const std::map<int,int>* global_id_to_internal_id_;
  
public:
  /** \brief constructor - with diffusion operator A^{\epsilon}(x)
   * \param subgrid_list cannot be const because Dune::Fem does not provide Gridparts that can be build on a const grid
   * \param DSC_LOG_INFO does not take ownership
   **/
  MsFEMLocalProblemSolver(const HostDiscreteFunctionSpaceType& hostDiscreteFunctionSpace,
                          const MacroMicroGridSpecifierType& specifier,
                          SubGridListType& subgrid_list,
                          const DiffusionOperatorType& diffusion_operator)
    : hostDiscreteFunctionSpace_(hostDiscreteFunctionSpace)
      , diffusion_(diffusion_operator)
      , specifier_(specifier)
      , subgrid_list_(subgrid_list)
  {}

  MsFEMLocalProblemSolver(const HostDiscreteFunctionSpaceType& hostDiscreteFunctionSpace,
                          const MacroMicroGridSpecifierType& specifier,
                          SubGridListType& subgrid_list,
                          const DiffusionOperatorType& diffusion_operator,
                          const CoarseBasisFunctionListType& coarse_basis,
                          const std::map<int,int>& global_id_to_internal_id )
    : hostDiscreteFunctionSpace_(hostDiscreteFunctionSpace)
      , diffusion_(diffusion_operator)
      , specifier_(specifier)
      , subgrid_list_(subgrid_list)
      , coarse_basis_( &coarse_basis )
      , global_id_to_internal_id_( &global_id_to_internal_id )
  {}
  
  template< class Stream >
  void oneLinePrint(Stream& stream, const SubDiscreteFunctionType& func) const {
    typedef typename SubDiscreteFunctionType::ConstDofIteratorType
    DofIteratorType;
    DofIteratorType it = func.dbegin();
    stream << "\n" << func.name() << ": [ ";
    for ( ; it != func.dend(); ++it)
      stream << std::setw(5) << *it << std::endl;// "  ";

    stream << " ] " << std::endl;
  } // oneLinePrint

  //! ----------- method: solve the local MsFEM problem ------------------------------------------

  void solvelocalproblem(JacobianRangeType& e,
                         SubDiscreteFunctionType& local_problem_solution,
                         const int coarse_index = -1 ) const {
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
    LocalProblemOperatorType local_problem_op(subDiscreteFunctionSpace, diffusion_);

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
    if ( specifier_.getOversamplingStrategy() == 1 )
    { local_problem_op.assemble_matrix(locprob_system_matrix); }

    if ( specifier_.getOversamplingStrategy() == 2 )
    { if ( coarse_index < 0 )
        DUNE_THROW(Dune::InvalidStateException, "Invalid coarse index: coarse_index < 0");
      local_problem_op.assemble_matrix(locprob_system_matrix, subgrid_list_.getCoarseNodeVector( coarse_index ) );
    }

    if ( specifier_.getOversamplingStrategy() == 3 )
    {
      if ( coarse_index < 0 )
        DUNE_THROW(Dune::InvalidStateException, "Invalid coarse index: coarse_index < 0");
      bool clement = ( DSC_CONFIG_GET( "rigorous_msfem.oversampling_strategy", "Clement" ) == "Clement" );

      if ( clement )
      { local_problem_op.assemble_matrix( locprob_system_matrix ); }
      else
      { local_problem_op.assemble_matrix( locprob_system_matrix, subgrid_list_.getCoarseNodeVector( coarse_index ) ); }
    }
    // can be deleted (just to check the coarse node vector)
    /*
    for ( int coarse_node_local_id = 0; coarse_node_local_id < subgrid_list_.getCoarseNodeVector( coarse_index ).size(); ++coarse_node_local_id )
       {
         std::cout << coarse_node_local_id+1 << " : " << subgrid_list_.getCoarseNodeVector( coarse_index )[coarse_node_local_id] << std :: endl << std :: endl << std :: endl;
       }
*/
    
    //! boundary treatment:
    typedef typename LocProbFEMMatrix::LocalMatrixType LocalMatrix;

    typedef typename SGLagrangePointSetType::template Codim< faceCodim >::SubEntityIteratorType
        FaceDofIteratorType;

    const HostGridPartType& hostGridPart = hostDiscreteFunctionSpace_.gridPart();

    const SubgridIteratorType sg_end = subDiscreteFunctionSpace.end();
    for (SubgridIteratorType sg_it = subDiscreteFunctionSpace.begin(); sg_it != sg_end; ++sg_it)
    {
      const SubgridEntityType& subgrid_entity = *sg_it;

      HostEntityPointerType host_entity_pointer = subGrid.template getHostEntity< 0 >(subgrid_entity);
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
          if ( subGrid.template contains< 0 >(neighborHostEntity) )
          {
            continue;
          }
        }

        const int face = (*iit).indexInInside();
        const FaceDofIteratorType fdend = lagrangePointSet.template endSubEntity< 1 >(face);
        for (FaceDofIteratorType fdit = lagrangePointSet.template beginSubEntity< 1 >(face); fdit != fdend; ++fdit)
          local_matrix.unitRow(*fdit);
      }
    }


    // assemble right hand side of algebraic local msfem problem
    if ( specifier_.getOversamplingStrategy() == 1 )
    { local_problem_op.assemble_local_RHS(e, local_problem_rhs); }
    else if ( ( specifier_.getOversamplingStrategy() == 2 ) ||
              ( specifier_.getOversamplingStrategy() == 3 ) )
    { if ( coarse_index < 0 )
        DUNE_THROW(Dune::InvalidStateException, "Invalid coarse index: coarse_index < 0");
      local_problem_op.assemble_local_RHS(e, subgrid_list_.getCoarseNodeVector( coarse_index ), specifier_.getOversamplingStrategy(), local_problem_rhs ); }
    else
      DUNE_THROW(Dune::InvalidStateException, "Oversampling Strategy must be 1, 2 or 3!");
    //oneLinePrint( DSC_LOG_DEBUG, local_problem_rhs );
    
    // zero boundary condition for 'cell problems':
    // set Dirichlet Boundary to zero
    for (SubgridIteratorType sg_it = subDiscreteFunctionSpace.begin(); sg_it != sg_end; ++sg_it)
    {
      const SubgridEntityType& subgrid_entity = *sg_it;

      HostEntityPointerType host_entity_pointer = subGrid.template getHostEntity< 0 >(subgrid_entity);
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

          if ( subGrid.template contains< 0 >(neighborHostEntity) )
          {
            continue;
          }
        }

        SubLocalFunctionType rhsLocal = local_problem_rhs.localFunction(subgrid_entity);
        const SGLagrangePointSetType& lagrangePointSet
            = subDiscreteFunctionSpace.lagrangePointSet(subgrid_entity);

        const int face = (*iit).indexInInside();

        FaceDofIteratorType faceIterator
            = lagrangePointSet.template beginSubEntity< faceCodim >(face);
        const FaceDofIteratorType faceEndIterator
            = lagrangePointSet.template endSubEntity< faceCodim >(face);
        for ( ; faceIterator != faceEndIterator; ++faceIterator)
          rhsLocal[*faceIterator] = 0;
      }
    }

    // After boundary treatment:
    // oneLinePrint( DSC_LOG_DEBUG, local_problem_rhs );

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
      InverseLocProbFEMMatrix locprob_fem_biCGStab(locprob_system_matrix, 1e-8, 1e-8, 20000, DSC_CONFIG_GET("LOCPROBLEMSOLVER_VERBOSE", false));
      
      bool clement = false;
      if ( specifier_.getOversamplingStrategy() == 3 )
      { clement = (DSC_CONFIG_GET( "rigorous_msfem.oversampling_strategy", "Clement" ) == "Clement" ); }

      if ( clement )
      {

        HostDiscreteFunctionType zero("zero", specifier_.coarseSpace());
        zero.clear();
        const double dummy = 12345.67890;
        double solverEps = 1e-8 ;
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
      else
      { locprob_fem_biCGStab(local_problem_rhs, local_problem_solution); }
    }

    if ( !( local_problem_solution.dofsValid() ) )
    {
      DUNE_THROW(Dune::InvalidStateException,"Current solution of the local msfem problem invalid!");
    }

    // oneLinePrint( DSC_LOG_DEBUG, local_problem_solution );
  } // solvelocalproblem

  //! ----------- end method: solve local MsFEM problem ------------------------------------------

  // create a hostgrid function from a subgridfunction
  void subgrid_to_hostrid_function(const SubDiscreteFunctionType& sub_func,
                                   HostDiscreteFunctionType& host_func) {
    host_func.clear();

    const SubDiscreteFunctionSpaceType& subDiscreteFunctionSpace = sub_func.space();
    const SubGridType& subGrid = subDiscreteFunctionSpace.grid();

    SubgridIteratorType sub_endit = subDiscreteFunctionSpace.end();
    for (SubgridIteratorType sub_it = subDiscreteFunctionSpace.begin(); sub_it != sub_endit; ++sub_it)
    {
      const SubgridEntityType& sub_entity = *sub_it;

      HostEntityPointerType host_entity_pointer = subGrid.template getHostEntity< 0 >(*sub_it);
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

  void output_local_solution(const int coarse_index, const int which,
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

  // method for solving and saving the solutions of the local msfem problems
  // for the whole set of macro-entities and for every unit vector e_i

  //! ---- method: solve and save the whole set of local msfem problems -----

  // Use the host-grid entities of Level 'computational_level' as computational domains for the subgrid computations
  void assemble_all(bool /*silent*/ = true /* state information on subgrids */) {
    std::string local_path = DSC_CONFIG_GET("global.datadir", "data") + "/local_problems/";
    Dune::Stuff::Common::testCreateDirectory(local_path);

    enum { dimension = MsFEMTraits::GridType::dimension };
    enum { maxnumOfBaseFct = 100 };

    JacobianRangeType e[dimension];
    for (int i = 0; i < dimension; ++i)
      for (int j = 0; j < dimension; ++j)
      {
        if (i == j)
        { e[i][0][j] = 1.0; } else
        { e[i][0][j] = 0.0; }
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


    const auto& comm = Dune::MPIHelper::getCollectiveCommunication();
    int slice = coarse_indices.size() / comm.size();
    for(int gc = comm.rank() * slice; gc < std::min(long(comm.rank() +1)* slice, long(coarse_indices.size())); ++gc)
    {
      const int coarse_index = coarse_indices[gc];

      DSC_LOG_INFO << "-------------------------" << std::endl
                   << "Coarse index " << coarse_index << std::endl;

      const std::string locprob_solution_location =
          (boost::format("local_problems/_localProblemSolutions_%d") % coarse_index).str();

      DiscreteFunctionWriter dfw(locprob_solution_location);

      SubGridType& subGrid = subgrid_list_.getSubGrid(coarse_index);

      SubGridPartType subGridPart(subGrid);

      DSC_LOG_INFO  << std::endl
                    << "Number of the local problem: " << dimension * coarse_index << " (of "
                    << (dimension * number_of_coarse_grid_entities) - 1 << " problems in total)" << std::endl
                    << "   Subgrid " << coarse_index << " contains " << subGrid.size(0) << " elements and "
                    << subGrid.size(2) << " nodes." << std::endl;

      const SubDiscreteFunctionSpaceType subDiscreteFunctionSpace(subGridPart);

      const std::string name_local_solution = (boost::format("Local Problem Solution %d") % coarse_index).str();

      //! only for dimension 2!
      SubDiscreteFunctionType local_problem_solution_0(name_local_solution, subDiscreteFunctionSpace);
      local_problem_solution_0.clear();

      SubDiscreteFunctionType local_problem_solution_1(name_local_solution, subDiscreteFunctionSpace);
      local_problem_solution_1.clear();

      // take time
      DSC_PROFILER.startTiming("none.local_problem_solution");

      // solve the problems
      solvelocalproblem(e[0], local_problem_solution_0, coarse_index);

      cell_time(DSC_PROFILER.stopTiming("none.local_problem_solution") / 1000.f);
      DSC_PROFILER.resetTiming("none.local_problem_solution");

      dfw.append(local_problem_solution_0);

      HostDiscreteFunctionType host_local_solution(name_local_solution, hostDiscreteFunctionSpace_);
      subgrid_to_hostrid_function(local_problem_solution_0, host_local_solution);
      output_local_solution(coarse_index, 0, host_local_solution);

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

      dfw.append(local_problem_solution_1);

      subgrid_to_hostrid_function(local_problem_solution_1, host_local_solution);
      output_local_solution(coarse_index, 1, host_local_solution);
    } //for

    const auto total_time = DSC_PROFILER.stopTiming("msfem.localproblemsolver.assemble_all");
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
}; // end class

} //namespace MsFEM {
} //namespace Multiscale {
} //namespace Dune {

#endif // #ifndef DiscreteEllipticMsFEMLocalProblem_HH
