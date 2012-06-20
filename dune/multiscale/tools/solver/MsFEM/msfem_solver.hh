#ifndef Elliptic_MSEM_Solver_HH
#define Elliptic_MSEM_Solver_HH

#include <dune/common/fmatrix.hh>

#include <dune/subgrid/subgrid.hh>

#include <dune/fem/solver/oemsolver/oemsolver.hh>
#include <dune/fem/solver/inverseoperators.hh>
#include <dune/fem/operator/discreteoperatorimp.hh>

#include <dune/fem/operator/2order/lagrangematrixsetup.hh>
#include <dune/fem/operator/matrix/spmatrix.hh>

#include <dune/fem/space/common/adaptmanager.hh>

#include <dune/multiscale/tools/assembler/righthandside_assembler.hh>
#include <dune/multiscale/tools/solver/MsFEM/msfem_localproblems/subgrid-list.hh>

#include <dune/multiscale/tools/assembler/matrix_assembler/elliptic_msfem_matrix_assembler.hh>

#include <dune/multiscale/tools/misc/linear-lagrange-interpolation.hh>

// / done

namespace Dune {
template< class DiscreteFunctionSpaceType >
class MacroMicroGridSpecifier
{
  typedef typename DiscreteFunctionSpaceType::RangeType RangeType;

public:
  MacroMicroGridSpecifier(DiscreteFunctionSpaceType& coarse_scale_space,
                          DiscreteFunctionSpaceType& fine_scale_space)
    : coarse_scale_space_(coarse_scale_space)
      , fine_scale_space_(fine_scale_space)
      , coarse_level_fine_level_difference_( fine_scale_space.gridPart().grid().maxLevel()
                                             - coarse_scale_space.gridPart().grid().maxLevel() )
      , // number of coarse space entities:
      number_of_level_host_entities_( coarse_scale_space.gridPart().grid().size(0 /*codim*/) ) {
    int number_of_coarse_entities = coarse_scale_space.gridPart().grid().size(0 /*codim*/);

    for (int i = 0; i < number_of_coarse_entities; i += 1)
    {
      // initialize with 0 layers:
      number_of_layers.push_back(0);
    }
  }

  // get number of coarse grid entities
  int getNumOfCoarseEntities() {
    return number_of_level_host_entities_;
  }

  void setLayer(int i, int number_of_layers_for_entity) {
    if (i < number_of_level_host_entities_)
    { number_of_layers[i] = number_of_layers_for_entity; } else {
      std::cout << "Error. Assertion (i < number_of_level_host_entities_) not filfilled." << std::endl;
      abort();
    }
  } // setLayer

  int getLayer(int i) {
    if (i < number_of_level_host_entities_)
    { return number_of_layers[i]; } else {
      std::cout << "Error. Assertion (i < number_of_level_host_entities_) not filfilled." << std::endl;
      abort();
    }
    return 0;
  } // getLayer

  // difference between coarse and fine level
  int getLevelDifference() {
    return coarse_level_fine_level_difference_;
  }

  // the coarse space
  DiscreteFunctionSpaceType& coarseSpace() {
    return coarse_scale_space_;
  }

  // the coarse space
  DiscreteFunctionSpaceType& fineSpace() {
    return fine_scale_space_;
  }

  void initialize_local_error_manager() {
    for (int i = 0; i < number_of_level_host_entities_; ++i)
    {
      loc_coarse_residual_.push_back(0.0);
      loc_projection_error_.push_back(0.0);
      loc_coarse_grid_jumps_.push_back(0.0);
      loc_conservative_flux_jumps_.push_back(0.0);
      loc_approximation_error_.push_back(0.0);
      loc_fine_grid_jumps_.push_back(0.0);
    }
  } // initialize_local_error_manager

  void set_loc_coarse_residual(int& index, RangeType& loc_coarse_residual) {
    loc_coarse_residual_[index] = loc_coarse_residual;
  }

  void set_loc_coarse_grid_jumps(int& index, RangeType& loc_coarse_grid_jumps) {
    loc_coarse_grid_jumps_[index] = loc_coarse_grid_jumps;
  }

  void set_loc_projection_error(int& index, RangeType& loc_projection_error) {
    loc_projection_error_[index] = loc_projection_error;
  }

  void set_loc_conservative_flux_jumps(int& index, RangeType& loc_conservative_flux_jumps) {
    loc_conservative_flux_jumps_[index] = loc_conservative_flux_jumps;
  }

  void set_loc_approximation_error(int& index, RangeType& loc_approximation_error) {
    loc_approximation_error_[index] = loc_approximation_error;
  }

  void set_loc_fine_grid_jumps(int& index, RangeType& loc_fine_grid_jumps) {
    loc_fine_grid_jumps_[index] = loc_fine_grid_jumps;
  }

  RangeType get_loc_coarse_residual(int& index) {
    if (loc_coarse_residual_.size() == 0)
    {
      std::cout << "Error! Use: initialize_local_error_manager()-method for the grid specifier first!" << std::endl;
    } else
    { return loc_coarse_residual_[index]; }
  } // get_loc_coarse_residual

  RangeType get_loc_coarse_grid_jumps(int& index) {
    if (loc_coarse_grid_jumps_.size() == 0)
    {
      std::cout << "Error! Use: initialize_local_error_manager()-method for the grid specifier first!" << std::endl;
    } else
    { return loc_coarse_grid_jumps_[index]; }
  } // get_loc_coarse_grid_jumps

  RangeType get_loc_projection_error(int& index) {
    if (loc_projection_error_.size() == 0)
    {
      std::cout << "Error! Use: initialize_local_error_manager()-method for the grid specifier first!" << std::endl;
    } else
    { return loc_projection_error_[index]; }
  } // get_loc_projection_error

  RangeType get_loc_conservative_flux_jumps(int& index) {
    if (loc_conservative_flux_jumps_.size() == 0)
    {
      std::cout << "Error! Use: initialize_local_error_manager()-method for the grid specifier first!" << std::endl;
    } else
    { return loc_conservative_flux_jumps_[index]; }
  } // get_loc_conservative_flux_jumps

  RangeType get_loc_approximation_error(int& index) {
    if (loc_approximation_error_.size() == 0)
    {
      std::cout << "Error! Use: initialize_local_error_manager()-method for the grid specifier first!" << std::endl;
    } else
    { return loc_approximation_error_[index]; }
  } // get_loc_approximation_error

  RangeType get_loc_fine_grid_jumps(int& index) {
    if (loc_fine_grid_jumps_.size() == 0)
    {
      std::cout << "Error! Use: initialize_local_error_manager()-method for the grid specifier first!" << std::endl;
    } else
    { return loc_fine_grid_jumps_[index]; }
  } // get_loc_fine_grid_jumps

private:
  DiscreteFunctionSpaceType& coarse_scale_space_;
  DiscreteFunctionSpaceType& fine_scale_space_;

  // level difference bettween coarse grid level and fine grid level
  int coarse_level_fine_level_difference_;

  // number of coarse grid entities
  int number_of_level_host_entities_;

  // layers for each coarse grid entity
  std::vector< int > number_of_layers;

  // ----- local error indicators (for each coarse grid element T) -------------

  // local coarse residual, i.e. H ||f||_{L^2(T)}
  std::vector< RangeType > loc_coarse_residual_;

  // local coarse grid jumps (contribute to the total coarse residual)
  std::vector< RangeType > loc_coarse_grid_jumps_;

  // local projection error (we project to get a globaly continous approximation)
  std::vector< RangeType > loc_projection_error_;

  // local jump in the conservative flux
  std::vector< RangeType > loc_conservative_flux_jumps_;

  // local approximation error
  std::vector< RangeType > loc_approximation_error_;

  // local sum over the fine grid jumps (for a fixed subgrid that cooresponds with a coarse entity T)
  std::vector< RangeType > loc_fine_grid_jumps_;
};

template< class DiscreteFunctionType >
class Elliptic_MsFEM_Solver
{
public:
  typedef DiscreteFunctionType DiscreteFunction;

  typedef typename DiscreteFunction::FunctionSpaceType FunctionSpace;

  typedef typename DiscreteFunction::DiscreteFunctionSpaceType DiscreteFunctionSpace;

  typedef typename DiscreteFunction::LocalFunctionType LocalFunction;

  typedef typename DiscreteFunctionSpace::LagrangePointSetType LagrangePointSet;

  typedef typename DiscreteFunctionSpace::GridPartType GridPart;

  typedef typename DiscreteFunctionSpace::GridType HostGrid;

  typedef typename HostGrid::Traits::LeafIndexSet HostGridLeafIndexSet;

  typedef typename HostGrid::Traits::LeafIndexSet CoarseGridLeafIndexSet;

  typedef typename DiscreteFunctionSpace::DomainType        DomainType;
  typedef typename DiscreteFunctionSpace::RangeType         RangeType;
  typedef typename DiscreteFunctionSpace::JacobianRangeType JacobianRangeType;

  // typedef typename HostGrid ::template Codim< 0 > :: template Partition< All_Partition > :: LevelIterator
  // LevelEntityIteratorType;

  typedef typename DiscreteFunctionSpace::IteratorType HostgridIterator;

  typedef typename HostgridIterator::Entity HostEntity;

  typedef typename HostEntity::EntityPointer HostEntityPointer;

  // typedef typename HostGrid :: template Codim< 0 > :: template Partition< All_Partition > :: LevelIterator
  // HostGridLevelEntityIterator;

  enum { faceCodim = 1 };

  typedef typename GridPart::IntersectionIteratorType IntersectionIterator;

  typedef typename LagrangePointSet::template Codim< faceCodim >
    ::SubEntityIteratorType
  FaceDofIterator;

  // --------------------------- subgrid typedefs ------------------------------------

  typedef SubGrid< HostGrid::dimension, HostGrid > SubGridType;

  typedef LeafGridPart< SubGridType > SubGridPart;

  // typedef typename SubGridType ::template Codim< 0 > :: template Partition< All_Partition > :: LevelIterator
  // SubGridLevelEntityIteratorType;

  typedef LagrangeDiscreteFunctionSpace< FunctionSpace, SubGridPart, 1 >  // 1=POLORDER
  SubgridDiscreteFunctionSpace;

  typedef AdaptiveDiscreteFunction< SubgridDiscreteFunctionSpace > SubgridDiscreteFunction;

  typedef typename SubgridDiscreteFunctionSpace::IteratorType CoarseGridIterator;

  typedef typename CoarseGridIterator::Entity CoarseGridEntity;

  typedef typename CoarseGridEntity::EntityPointer CoarseGridEntityPointer;

  typedef typename SubgridDiscreteFunction::LocalFunctionType CoarseGridLocalFunction;

  typedef typename SubgridDiscreteFunctionSpace::LagrangePointSetType
  CoarseGridLagrangePointSet;

  typedef typename CoarseGridLagrangePointSet::template Codim< faceCodim >
    ::SubEntityIteratorType
  CoarseGridFaceDofIterator;

  // !-----------------------------------------------------------------------------------------

  // ! --------------------- the standard matrix traits -------------------------------------

  struct MatrixTraits
  {
    typedef SubgridDiscreteFunctionSpace                          RowSpaceType;
    typedef SubgridDiscreteFunctionSpace                          ColumnSpaceType;
    typedef LagrangeMatrixSetup< false >                          StencilType;
    typedef ParallelScalarProduct< SubgridDiscreteFunctionSpace > ParallelScalarProductType;

    template< class M >
    struct Adapter
    {
      typedef LagrangeParallelMatrixAdapter< M > MatrixAdapterType;
    };
  };

  // ! --------------------------------------------------------------------------------------

  // ! --------------------- type of fem stiffness matrix -----------------------------------

  typedef SparseRowMatrixOperator< DiscreteFunction, DiscreteFunction, MatrixTraits > MsFEMMatrix;

  // ! --------------------------------------------------------------------------------------

  // ! --------------- solver for the linear system of equations ----------------------------

  // use Bi CG Stab [OEMBICGSTABOp] or GMRES [OEMGMRESOp] for non-symmetric matrices and CG [CGInverseOp] for symmetric
  // ones. GMRES seems to be more stable, but is extremely slow!
  typedef /*OEMBICGSQOp*//*CGInverseOp*/ OEMBICGSTABOp< DiscreteFunction, MsFEMMatrix > InverseMsFEMMatrix;

  // ! --------------------------------------------------------------------------------------

private:
  const DiscreteFunctionSpace& discreteFunctionSpace_;

  std::ofstream* data_file_;

  // path where to save the data output
  std::string path_;

public:
  Elliptic_MsFEM_Solver(const DiscreteFunctionSpace& discreteFunctionSpace, std::string path = "")
    : discreteFunctionSpace_(discreteFunctionSpace)
      , data_file_(NULL)
  { path_ = path; }

  Elliptic_MsFEM_Solver(const DiscreteFunctionSpace& discreteFunctionSpace,
                        std::ofstream& data_file,
                        std::string path = "")
    : discreteFunctionSpace_(discreteFunctionSpace)
      , data_file_(&data_file)
  { path_ = path; }

  template< class Stream >
  void oneLinePrint(Stream& stream, const DiscreteFunction& func) {
    typedef typename DiscreteFunction::ConstDofIteratorType
    DofIteratorType;
    DofIteratorType it = func.dbegin();
    stream << "\n" << func.name() << ": [ ";
    for ( ; it != func.dend(); ++it)
      stream << std::setw(5) << *it << "  ";

    stream << " ] " << std::endl;
  } // oneLinePrint

  template< class Stream >
  void oneLinePrint(Stream& stream, const SubgridDiscreteFunction& func) {
    typedef typename SubgridDiscreteFunction::ConstDofIteratorType
    DofIteratorType;
    DofIteratorType it = func.dbegin();
    stream << "\n" << func.name() << ": [ ";
    for ( ; it != func.dend(); ++it)
      stream << std::setw(5) << *it << "  ";

    stream << " ] " << std::endl;
  } // oneLinePrint

  // create a hostgrid function from a subgridfunction (projection for global continuity)
  // Note: the maximum gride levels for both underlying grids must be the same
  void subgrid_to_hostrid_projection(const SubgridDiscreteFunction& sub_func, DiscreteFunction& host_func) {
    #if 0
    if ( sub_func.space().gridPart().grid().maxLevel() != host_func.space().gridPart().grid().maxLevel() )
    {
      std::cout << "Error in method 'subgrid_to_hostrid_function': MaxLevel of SubGrid = "
                << sub_func.space().gridPart().grid().maxLevel() << " not identical to MaxLevel of FineGrid = "
                << host_func.space().gridPart().grid().maxLevel() << "." << std::endl;
    }
    #endif // if 0

    host_func.clear();

    const SubgridDiscreteFunctionSpace& subDiscreteFunctionSpace = sub_func.space();
    const SubGridType& subGrid = subDiscreteFunctionSpace.grid();

    typedef typename SubgridDiscreteFunctionSpace::IteratorType SubgridIterator;
    typedef typename SubgridIterator::Entity                    SubgridEntity;
    typedef typename SubgridDiscreteFunction::LocalFunctionType SubgridLocalFunction;

    SubgridIterator sub_endit = subDiscreteFunctionSpace.end();
    for (SubgridIterator sub_it = subDiscreteFunctionSpace.begin(); sub_it != sub_endit; ++sub_it)
    {
      const SubgridEntity& sub_entity = *sub_it;

      HostEntityPointer host_entity_pointer = subGrid.template getHostEntity< 0 >(*sub_it);
      const HostEntity& host_entity = *host_entity_pointer;

      SubgridLocalFunction sub_loc_value = sub_func.localFunction(sub_entity);
      LocalFunction host_loc_value = host_func.localFunction(host_entity);

      const unsigned int numBaseFunctions = sub_loc_value.baseFunctionSet().size();
      for (unsigned int i = 0; i < numBaseFunctions; ++i)
      {
        host_loc_value[i] = sub_loc_value[i];
      }
    }
  } // subgrid_to_hostrid_projection

  // - ∇ (A(x,∇u)) + b ∇u + c u = f - divG
  // then:
  // A --> diffusion operator ('DiffusionOperatorType')
  // b --> advective part ('AdvectionTermType')
  // c --> reaction part ('ReactionTermType')
  // f --> 'first' source term, scalar ('SourceTermType')
  // G --> 'second' source term, vector valued ('SecondSourceTermType')

  // homogenous Dirchilet boundary condition!:
  template< class DiffusionOperator, class SourceTerm, class SubGridListType >
  void solve_dirichlet_zero(const DiffusionOperator& diffusion_op,
                            const SourceTerm& f,
                            // number of layers per coarse grid entity T:  U(T) is created by enrichting T with
                            // n(T)-layers.
                            MacroMicroGridSpecifier< DiscreteFunctionSpace >& specifier,
                            SubGridListType& subgrid_list,
                            DiscreteFunction& coarse_scale_part,
                            DiscreteFunction& fine_scale_part,
                            DiscreteFunction& solution) {
    // discrete elliptic MsFEM operator (corresponds with MsFEM Matrix)
    typedef DiscreteEllipticMsFEMOperator< DiscreteFunction /*type of coarse space*/,
                                           MacroMicroGridSpecifier< DiscreteFunctionSpace >,
                                           DiscreteFunction /*type of fine space*/,
                                           DiffusionOperator > EllipticMsFEMOperatorType;

    DiscreteFunctionSpace& coarse_space = specifier.coarseSpace();

    HostgridIterator coarse_iterator_end = coarse_space.end();
    HostgridIterator coarse_iterator_begin = coarse_space.begin();

    #if 0
    int i = 0;
    for (HostgridIterator coarse_it = coarse_iterator_begin; coarse_it != coarse_iterator_end; ++coarse_it)
    {
      std::cout << "i = " << i << std::endl;
      std::cout << "coarse element corner(0) = " << coarse_it->geometry().corner(0) << std::endl;
      std::cout << "coarse element corner(1) = " << coarse_it->geometry().corner(1) << std::endl;
      std::cout << "coarse element corner(2) = " << coarse_it->geometry().corner(2) << std::endl << std::endl;
      i += 1;
    }
    #endif // if 0

    HostGrid& grid = discreteFunctionSpace_.gridPart().grid();
    const GridPart& gridPart = discreteFunctionSpace_.gridPart();

    // ------------------------------------------------------------

    DiscreteFunction coarse_msfem_solution("Coarse Part MsFEM Solution", coarse_space);
    coarse_msfem_solution.clear();

    // ! create subgrids:
    bool DUNE_UNUSED(silence) = false;

    // ! define the right hand side assembler tool
    // (for linear and non-linear elliptic and parabolic problems, for sources f and/or G )
    RightHandSideAssembler< DiscreteFunction > rhsassembler;

    // ! define the discrete (elliptic) operator that describes our problem
    // ( effect of the discretized differential operator on a certain discrete function )
    EllipticMsFEMOperatorType elliptic_msfem_op(specifier,
                                                coarse_space,
                                                subgrid_list,
                                                diffusion_op, *data_file_, path_);
    // discrete elliptic operator (corresponds with FEM Matrix)

    // ! (stiffness) matrix
    MsFEMMatrix msfem_matrix("MsFEM stiffness matrix", coarse_space, coarse_space);

    // ! right hand side vector
    // right hand side for the finite element method:
    DiscreteFunction msfem_rhs("MsFEM right hand side", coarse_space);
    msfem_rhs.clear();

    std::cout << std::endl;
    std::cout << "Solving MsFEM problem." << std::endl;

    if (data_file_)
    {
      if ( data_file_->is_open() )
      {
        *data_file_ << "Solving linear problem with MsFEM and maximum coarse grid level "
                    << coarse_space.gridPart().grid().maxLevel() << "." << std::endl;
        *data_file_ << "------------------------------------------------------------------------------" << std::endl;
      }
    }

    // to assemble the computational time
    Dune::Timer assembleTimer;

    // assemble the MsFEM stiffness matrix
    elliptic_msfem_op.assemble_matrix(msfem_matrix);   // einbinden!

    std::cout << "Time to assemble MsFEM stiffness matrix: " << assembleTimer.elapsed() << "s" << std::endl;

    if (data_file_)
    {
      if ( data_file_->is_open() )
      {
        *data_file_ << "Time to assemble MsFEM stiffness matrix: " << assembleTimer.elapsed() << "s" << std::endl;
      }
    }

    // assemble right hand side
    rhsassembler.template assemble< 2* DiscreteFunctionSpace::polynomialOrder + 2 >(f, msfem_rhs);

    // oneLinePrint( std::cout , fem_rhs );

    #if 1
    // --- boundary treatment ---
    // set the dirichlet points to zero (in righ hand side of the fem problem)
    HostgridIterator endit = coarse_space.end();
    for (HostgridIterator it = coarse_space.begin(); it != endit; ++it)
    {
      IntersectionIterator iit = coarse_space.gridPart().ibegin(*it);
      const IntersectionIterator endiit = coarse_space.gridPart().iend(*it);
      for ( ; iit != endiit; ++iit)
      {
        if ( !(*iit).boundary() )
          continue;

        LocalFunction rhsLocal = msfem_rhs.localFunction(*it);

        const LagrangePointSet& lagrangePointSet
          = coarse_space.lagrangePointSet(*it);

        const int face = (*iit).indexInInside();

        FaceDofIterator faceIterator
          = lagrangePointSet.template beginSubEntity< faceCodim >(face);
        const FaceDofIterator faceEndIterator
          = lagrangePointSet.template endSubEntity< faceCodim >(face);
        for ( ; faceIterator != faceEndIterator; ++faceIterator)
          rhsLocal[*faceIterator] = 0;
      }
    }
    // --- end boundary treatment ---
    #endif // if 1

    InverseMsFEMMatrix msfem_biCGStab(msfem_matrix, 1e-8, 1e-8, 20000, true /*VERBOSE*/);
    msfem_biCGStab(msfem_rhs, coarse_msfem_solution);

    if (data_file_)
    {
      if ( data_file_->is_open() )
      {
        *data_file_ << "---------------------------------------------------------------------------------" << std::endl;
        *data_file_ << "MsFEM problem solved in " << assembleTimer.elapsed() << "s." << std::endl << std::endl
                    << std::endl;
      }
    }

    // oneLinePrint( std::cout , solution );

    // copy coarse grid function (defined on the subgrid) into a fine grid function
    solution.clear();

    const HostGridLeafIndexSet& coarseGridLeafIndexSet = coarse_space.gridPart().grid().leafIndexSet();

    #if 1
    std::cout << "Indentifying coarse scale part of the MsFEM solution... ";

    // ! copy coarse scale part of MsFEM solution into a function defined on the fine grid
    // ------------------------------------------------------------------------------------
    typedef typename HostEntity::template Codim< 2 >::EntityPointer HostNodePointer;

    typedef typename GridPart::IntersectionIteratorType HostIntersectionIterator;

    for (HostgridIterator it = discreteFunctionSpace_.begin(); it != discreteFunctionSpace_.end(); ++it)
    {
      typedef typename HostEntity::template Codim< 0 >::EntityPointer
      HostEntityPointer;
      HostEntityPointer coarse_father(*it);
      for (int lev = 0; lev < specifier.getLevelDifference(); ++lev)
        coarse_father = coarse_father->father();

      // ! new version:
      #if 1
      typename HostEntity::template Codim< 0 >::EntityPointer coarse_father_test = coarse_father;

      bool father_found = false;
      while (father_found == false)
      {
        if (coarseGridLeafIndexSet.contains(*coarse_father_test) == true)
        {
          coarse_father = coarse_father_test;
        }

        if (coarse_father_test->hasFather() == false)
        { father_found = true; } else
        { coarse_father_test = coarse_father_test->father(); }
      }
      #endif // if 1

      // ! old version
      #if 0
      bool father_found = false;
      while (father_found == false)
      {
        if (coarseGridLeafIndexSet.contains(*coarse_father) == true)
        {
          if (coarse_father->hasFather() == false)
          { father_found = true; } else {
            if (coarseGridLeafIndexSet.contains( *( coarse_father->father() ) ) == false)
            { father_found = true; }
          }
        }

        if (father_found == false)
        { coarse_father = coarse_father->father(); }
      }
      #endif // if 0

      #if 0
      int index_ce = coarseGridLeafIndexSet.index(*coarse_father);
      std::cout << "index_ce = " << index_ce << std::endl;
      std::cout << "fine element corner(0) = " << it->geometry().corner(0) << std::endl;
      std::cout << "fine element corner(1) = " << it->geometry().corner(1) << std::endl;
      std::cout << "fine element corner(2) = " << it->geometry().corner(2) << std::endl << std::endl;
      std::cout << "coarse element corner(0) = " << coarse_father->geometry().corner(0) << std::endl;
      std::cout << "coarse element corner(1) = " << coarse_father->geometry().corner(1) << std::endl;
      std::cout << "coarse element corner(2) = " << coarse_father->geometry().corner(2) << std::endl << std::endl;
      #endif // if 0

      LinearLagrangeFunction2D< DiscreteFunctionSpace > interpolation_coarse(coarse_father);

      interpolation_coarse.set_corners(coarse_msfem_solution);

      LocalFunction host_loc_value = coarse_scale_part.localFunction(*it);

      const int number_of_nodes = (*it).template count< 2 >();
      if ( !( number_of_nodes == int( host_loc_value.baseFunctionSet().size() ) ) )
      { std::cout << "Error! Inconsistency in 'msfem_solver.hh'." << std::endl; }

      for (int i = 0; i < number_of_nodes; i += 1)
      {
        const HostNodePointer node = (*it).template subEntity< 2 >(i);

        DomainType coordinates_of_node = node->geometry().corner(0);
        if ( !( coordinates_of_node == it->geometry().corner(i) ) )
        { std::cout << "Error! Inconsistency in 'msfem_solver.hh'." << std::endl; }

        RangeType coarse_value(0.0);
        interpolation_coarse.evaluate(coordinates_of_node, coarse_value);

        // int global_index_node = gridPart.indexSet().index( *node );
        host_loc_value[i] = coarse_value;
      }
    }
    std::cout << " done." << std::endl;
    #endif // if 1
    // ------------------------------------------------------------------------------------

    fine_scale_part.clear();

    int number_of_nodes = grid.size(2 /*codim*/);
    std::vector< std::vector< HostEntityPointer > > entities_sharing_same_node(number_of_nodes);

    for (HostgridIterator it = discreteFunctionSpace_.begin(); it != discreteFunctionSpace_.end(); ++it)
    {
      int number_of_nodes_in_entity = (*it).template count< 2 >();
      for (int i = 0; i < number_of_nodes_in_entity; i += 1)
      {
        const typename HostEntity::template Codim< 2 >::EntityPointer node = (*it).template subEntity< 2 >(i);
        int global_index_node = gridPart.indexSet().index(*node);

        entities_sharing_same_node[global_index_node].push_back( HostEntityPointer(*it) );
      }
    }

    // ! indentify fine scale part of MsFEM solution (including the projection!)
    // ------------------------------------------------------------------------------------

    std::cout << "Indentifying fine scale part of the MsFEM solution... ";
    // iterator ueber coarse space
    for (HostgridIterator coarse_it = coarse_space.begin(); coarse_it != coarse_space.end(); ++coarse_it)
    {
      // the coarse entity 'T': *coarse_it

      DiscreteFunction correction_on_U_T("correction_on_U_T", discreteFunctionSpace_);
      correction_on_U_T.clear();

      int index = coarseGridLeafIndexSet.index(*coarse_it);

      // the sub grid U(T) that belongs to the coarse_grid_entity T
      SubGridType& sub_grid_U_T = subgrid_list.getSubGrid(index);
      SubGridPart subGridPart(sub_grid_U_T);

      SubgridDiscreteFunctionSpace localDiscreteFunctionSpace(subGridPart);

      SubgridDiscreteFunction local_problem_solution_e0("Local problem Solution e_0", localDiscreteFunctionSpace);
      local_problem_solution_e0.clear();

      SubgridDiscreteFunction local_problem_solution_e1("Local problem Solution e_1", localDiscreteFunctionSpace);
      local_problem_solution_e1.clear();

      // --------- load local solutions -------

      char location_lps[50];
      sprintf(location_lps, "/local_problems/_localProblemSolutions_%d", index);
      std::string location_lps_s(location_lps);

      std::string local_solution_location;

      // the file/place, where we saved the solutions of the cell problems
      local_solution_location = path_ + location_lps_s;

      bool reader_is_open = false;
      // reader for the cell problem data file:
      DiscreteFunctionReader discrete_function_reader( (local_solution_location).c_str() );
      reader_is_open = discrete_function_reader.open();

      if (reader_is_open)
      { discrete_function_reader.read(0, local_problem_solution_e0); }

      if (reader_is_open)
      { discrete_function_reader.read(1, local_problem_solution_e1); }

      // oneLinePrint( std::cout , local_problem_solution_e0 );
      // oneLinePrint( std::cout , local_problem_solution_e1 );

      LocalFunction local_coarse_part = coarse_msfem_solution.localFunction(*coarse_it);

      // 1 point quadrature!! We only need the gradient of the coarse scale part on the element, which is a constant.
      CachingQuadrature< GridPart, 0 > one_point_quadrature(*coarse_it, 0);

      JacobianRangeType grad_coarse_msfem_on_entity;
      local_coarse_part.jacobian(one_point_quadrature[0], grad_coarse_msfem_on_entity);

      // !
      // std :: cout << "grad_coarse_msfem_on_entity[ 0 ][ 1 ] = " << grad_coarse_msfem_on_entity[ 0 ][ 1 ] << std ::
      // endl;
      // std :: cout << "grad_coarse_msfem_on_entity[ 0 ][ 0 ] = " << grad_coarse_msfem_on_entity[ 0 ][ 0 ] << std ::
      // endl;
      local_problem_solution_e0 *= grad_coarse_msfem_on_entity[0][0];
      local_problem_solution_e1 *= grad_coarse_msfem_on_entity[0][1];
      local_problem_solution_e0 += local_problem_solution_e1;

      // oneLinePrint( std::cout , local_problem_solution_e0 );

      subgrid_to_hostrid_projection(local_problem_solution_e0, correction_on_U_T);
      // hol die den Gradient und addiere.
      #if 1
      if ( sub_grid_U_T.maxLevel() != discreteFunctionSpace_.gridPart().grid().maxLevel() )
      { std::cout << "Error: MaxLevel of SubGrid not identical to MaxLevel of FineGrid." << std::endl; }

      correction_on_U_T.clear();

      typedef typename SubgridDiscreteFunctionSpace::IteratorType SubgridIterator;
      typedef typename SubgridIterator::Entity                    SubgridEntity;
      typedef typename SubgridDiscreteFunction::LocalFunctionType SubgridLocalFunction;

      SubgridIterator sub_endit = localDiscreteFunctionSpace.end();
      for (SubgridIterator sub_it = localDiscreteFunctionSpace.begin(); sub_it != sub_endit; ++sub_it)
      {
        const SubgridEntity& sub_entity = *sub_it;

        HostEntityPointer fine_host_entity_pointer = sub_grid_U_T.template getHostEntity< 0 >(*sub_it);
        const HostEntity& fine_host_entity = *fine_host_entity_pointer;

        HostEntityPointer father = fine_host_entity_pointer;
        for (int lev = 0; lev < specifier.getLevelDifference(); ++lev)
          father = father->father();

        // ! new version:
        #if 1
        typename HostEntity::template Codim< 0 >::EntityPointer coarse_father_test = father;

        bool father_found = false;
        while (father_found == false)
        {
          if (coarseGridLeafIndexSet.contains(*coarse_father_test) == true)
          {
            father = coarse_father_test;
          }

          if (coarse_father_test->hasFather() == false)
          { father_found = true; } else
          { coarse_father_test = coarse_father_test->father(); }
        }
        #endif // if 1

        // ! old version:
        #if 0
        bool father_found = false;
        while (father_found == false)
        {
          if (coarseGridLeafIndexSet.contains(*father) == true)
          {
            if (father->hasFather() == false)
            { father_found = true; } else {
              if (coarseGridLeafIndexSet.contains( *( father->father() ) ) == false)
              { father_found = true; }
            }
          }

          if (father_found == false)
          { father = father->father(); }
        }
        #endif // if 0

        bool entities_identical = true;
        int number_of_nodes = (*father).template count< 2 >();
        for (int k = 0; k < number_of_nodes; k += 1)
        {
          if ( !( father->geometry().corner(k) == (*coarse_it).geometry().corner(k) ) )
          { entities_identical = false; }
        }

        if (entities_identical == false)
        {
          // std :: cout << "father->geometry().corner(0) = " << father->geometry().corner(0) << std :: endl;
          // std :: cout << "father->geometry().corner(1) = " << father->geometry().corner(1) << std :: endl;
          // std :: cout << "father->geometry().corner(2) = " << father->geometry().corner(2) << std :: endl;
          // std :: cout << "coarse_host_entity.geometry().corner(0) = " << coarse_host_entity.geometry().corner(0) <<
          // std :: endl;
          // std :: cout << "coarse_host_entity.geometry().corner(1) = " << coarse_host_entity.geometry().corner(1) <<
          // std :: endl;
          // std :: cout << "coarse_host_entity.geometry().corner(2) = " << coarse_host_entity.geometry().corner(2) <<
          // std :: endl << std :: endl;
          continue;
        }

        SubgridLocalFunction sub_loc_value = local_problem_solution_e0.localFunction(sub_entity);
        LocalFunction host_loc_value = correction_on_U_T.localFunction(fine_host_entity);

        int number_of_nodes_entity = (*sub_it).template count< 2 >();
        for (int i = 0; i < number_of_nodes_entity; i += 1)
        {
          const typename HostEntity::template Codim< 2 >::EntityPointer node = fine_host_entity.template subEntity< 2 >(
            i);

          int global_index_node = gridPart.indexSet().index(*node);

          // vector of coarse entities that share the above node
          std::vector< HostEntityPointer > coarse_entities;

          // count the number of different coarse-grid-entities that share the above node
          for (size_t j = 0; j < entities_sharing_same_node[global_index_node].size(); j += 1)
          {
            HostEntityPointer inner_it = entities_sharing_same_node[global_index_node][j];
            for (int lev = 0; lev < specifier.getLevelDifference(); ++lev)
              inner_it = inner_it->father();

            // ! new version:
            #if 1
            typename HostEntity::template Codim< 0 >::EntityPointer coarse_father_test = inner_it;

            bool father_found = false;
            while (father_found == false)
            {
              if (coarseGridLeafIndexSet.contains(*coarse_father_test) == true)
              { inner_it = coarse_father_test; }

              if (coarse_father_test->hasFather() == false)
              { father_found = true; } else
              { coarse_father_test = coarse_father_test->father(); }
            }
            #endif // if 1

            // ! old version:
            #if 0
            bool found = coarseGridLeafIndexSet.contains(*inner_it);
            while (found == false)
            {
              inner_it = inner_it->father();
              found = coarseGridLeafIndexSet.contains(*inner_it);
            }
            #endif // if 0

            bool new_entity_found = true;
            for (size_t k = 0; k < coarse_entities.size(); k += 1)
            {
              if (coarse_entities[k] == inner_it)
              { new_entity_found = false; }
            }
            if (new_entity_found == true)
            { coarse_entities.push_back(inner_it); }
          }

          host_loc_value[i] = ( sub_loc_value[i] / coarse_entities.size() );
        }
      }

      #endif // if 1

      fine_scale_part += correction_on_U_T;
    }
    std::cout << " done." << std::endl;

    // ------------------------------------------------------------------------------------

    // Auf Grobskalen MsFEM Anteil noch Feinksalen MsFEM Anteil aufaddieren.
    solution += coarse_scale_part;
    solution += fine_scale_part;
  } // solve_dirichlet_zero

  // ! the following methods are not yet implemented, however note that the required tools are
  // ! already available via 'righthandside_assembler.hh' and 'elliptic_fem_matrix_assembler.hh'!

  template< class DiffusionOperatorType, class ReactionTermType, class SourceTermType >
  void solve() {
    std::cout << "No implemented!" << std::endl;
  }

  template< class DiffusionOperatorType, class ReactionTermType, class SourceTermType, class SecondSourceTermType >
  void solve() {
    std::cout << "No implemented!" << std::endl;
  }
};
}

#endif // #ifndef Elliptic_MSEM_Solver_HH
