#include <config.h>

#include "lod_solver.hh"

#include <dune/common/fmatrix.hh>

#include <dune/fem/solver/oemsolver/oemsolver.hh>
#include <dune/fem/operator/discreteoperatorimp.hh>
#include <dune/fem/operator/2order/lagrangematrixsetup.hh>
#include <dune/fem/operator/matrix/spmatrix.hh>
#include <dune/fem/space/common/adaptmanager.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/function/common/function.hh>

#include <dune/multiscale/common/righthandside_assembler.hh>
#include <dune/multiscale/msfem/localproblems/subgrid-list.hh>
#include <dune/multiscale/tools/misc/linear-lagrange-interpolation.hh>
#include <dune/multiscale/msfem/localproblems/localproblemsolver.hh>
#include <dune/multiscale/msfem/msfem_grid_specifier.hh>
#include <dune/multiscale/common/output_traits.hh>
#include <dune/multiscale/problems/base.hh>
#include <dune/multiscale/problems/selector.hh>

#include <dune/istl/matrix.hh>
#include <dune/stuff/fem/functions/checks.hh>

namespace Dune {
namespace Multiscale {
namespace MsFEM {

Elliptic_Rigorous_MsFEM_Solver::Elliptic_Rigorous_MsFEM_Solver(const DiscreteFunctionSpace& discreteFunctionSpace)
  : discreteFunctionSpace_(discreteFunctionSpace) {}

// create a hostgrid function from a subgridfunction (projection for global continuity)
// Note: the maximum gride levels for both underlying grids must be the same
void Elliptic_Rigorous_MsFEM_Solver::subgrid_to_hostrid_projection(const SubGridDiscreteFunctionType& sub_func,
                                                                   DiscreteFunction& host_func) const {
  host_func.clear();

  const SubGridDiscreteFunctionSpaceType& subDiscreteFunctionSpace = sub_func.space();
  const SubGridType& subGrid = subDiscreteFunctionSpace.grid();

  const SubGridIteratorType sub_endit = subDiscreteFunctionSpace.end();
  for (SubGridIteratorType sub_it = subDiscreteFunctionSpace.begin(); sub_it != sub_endit; ++sub_it) {
    const SubGridEntityType& sub_entity = *sub_it;

    const HostEntityPointer host_entity_pointer = subGrid.getHostEntity<0>(*sub_it);
    const HostEntity& host_entity = *host_entity_pointer;

    const SubGridLocalFunctionType sub_loc_value = sub_func.localFunction(sub_entity);
    auto host_loc_value = host_func.localFunction(host_entity);

    const auto numBaseFunctions = sub_loc_value.basisFunctionSet().size();
    for (unsigned int i = 0; i < numBaseFunctions; ++i) {
      host_loc_value[i] = sub_loc_value[i];
    }
  }
} // subgrid_to_hostrid_projection

// vtk visualization of msfem basis functions
void Elliptic_Rigorous_MsFEM_Solver::vtk_output(MsFEMBasisFunctionType& msfem_basis_function_list,
                                                std::string basis_name) const {

  // general output parameters
  Dune::Multiscale::OutputParameters outputparam(DSC_CONFIG_GET("global.datadir", "data") + "/msfem_basis");

  typedef typename OutputTraits::IOTupleType IOTType;
  const auto& gridPart = msfem_basis_function_list[0]->space().gridPart();

  for (size_t i = 0; i < msfem_basis_function_list.size(); i += 1) {

    IOTType msfem_basis_series(&(*msfem_basis_function_list[i]));

    const std::string ls_name_s = "/" + basis_name + (boost::format("_%d") % i).str();
    outputparam.set_prefix(ls_name_s);

    std::string outstring = basis_name;

    OutputTraits::DataOutputType msfem_basis_dataoutput(gridPart.grid(), msfem_basis_series, outputparam);
    msfem_basis_dataoutput.writeData(1.0 /*dummy*/, outstring);
  }

  std::cout << "VTK Output for MsFEM basis functions successful." << std::endl << std::endl;
}

// ------------------------------------------------------------------------------------------------------
//! for each subgrid, store the vector of basis functions ids that correspond to interior coarse grid nodes in the
// subgrid
// information stored in 'std::vector< std::vector< int > >'
void Elliptic_Rigorous_MsFEM_Solver::assemble_interior_basis_ids(
    MacroMicroGridSpecifier& specifier, SubGridListType& subgrid_list, std::map<int, int>& global_id_to_internal_id,
    std::map<OrderedDomainType, int>& coordinates_to_global_coarse_node_id,
    // all coarse nodes, where the corresponding coarse (nodal) basis function has a support that intersects with the
    // subgrid
    std::vector<std::vector<int>>& ids_basis_function_in_extended_subgrid,
    // all coarse nodes that are in the subgrid including its boundary but excluding the Dirichlet boundary of Omega
    std::vector<std::vector<int>>& ids_basis_function_in_subgrid,
    // all coarse nodes that are in the _interior_ of the subgrid
    std::vector<std::vector<int>>& ids_basis_function_in_interior_subgrid) const {
  // First: determine the index set (internal_coarse_nodes numbering) for the coarse nodes in the interior of U(T)
  // when assembling the local clement operator for a given subgrid U(T), we need to know the standard coarse
  // basis functions that belong to the interior(!) coarse nodes in U(T). Coarse nodes on the boundary of
  // are not relevant.

  int number_of_subgrids = subgrid_list.getNumberOfSubGrids();
  DiscreteFunctionSpace& fine_space = specifier.fineSpace();

  ids_basis_function_in_subgrid.resize(number_of_subgrids);
  ids_basis_function_in_extended_subgrid.resize(number_of_subgrids);

  typedef std::vector<DomainType> CoarseNodeVectorType;
  typedef std::vector<CoarseNodeVectorType> CoarseGridNodeStorageType;
  std::vector<std::vector<int>> coarse_node_ids_in_subgrid;
  std::vector<std::vector<int>> coarse_boundary_node_ids_in_subgrid; // including the boundary
  std::vector<std::vector<int>> coarse_interior_node_ids_in_subgrid;
  coarse_node_ids_in_subgrid.resize(number_of_subgrids);
  coarse_boundary_node_ids_in_subgrid.resize(number_of_subgrids);
  coarse_interior_node_ids_in_subgrid.resize(number_of_subgrids);

  for (unsigned int sg_id = 0; sg_id < number_of_subgrids; sg_id += 1) {
    CoarseNodeVectorType coarse_nodes_in_subgrid = subgrid_list.getCoarseNodeVector(sg_id);
    CoarseNodeVectorType coarse_nodes_in_extended_subgrid = subgrid_list.getExtendedCoarseNodeVector(sg_id);

    for (unsigned int cn = 0; cn < coarse_nodes_in_subgrid.size(); ++cn) {

      int global_coarse_node_id = coordinates_to_global_coarse_node_id[coarse_nodes_in_subgrid[cn]];
      coarse_node_ids_in_subgrid[sg_id].push_back(global_coarse_node_id);

      // sort out the Dirichlet boundary nodes on Omega
      if (specifier.is_coarse_dirichlet_node(global_coarse_node_id) == false) {
        ids_basis_function_in_subgrid[sg_id].push_back(global_id_to_internal_id[global_coarse_node_id]);
      }
    }

    for (unsigned int cn = 0; cn < coarse_nodes_in_extended_subgrid.size(); ++cn) {

      int global_coarse_node_id = coordinates_to_global_coarse_node_id[coarse_nodes_in_extended_subgrid[cn]];
      // sort out the Dirichlet boundary nodes on Omega
      if (specifier.is_coarse_dirichlet_node(global_coarse_node_id) == false) {
        ids_basis_function_in_extended_subgrid[sg_id].push_back(global_id_to_internal_id[global_coarse_node_id]);
      }
    }
  }

  for (unsigned int sg_id = 0; sg_id < number_of_subgrids; sg_id += 1) {
    auto subGridPart = subgrid_list.gridPart(sg_id);
    const SubGridDiscreteFunctionSpaceType subDiscreteFunctionSpace(subGridPart);

    CoarseNodeVectorType coarse_nodes_in_subgrid = subgrid_list.getCoarseNodeVector(sg_id);

    const SubGridIteratorType sg_end = subDiscreteFunctionSpace.end();
    for (SubGridIteratorType sg_it = subDiscreteFunctionSpace.begin(); sg_it != sg_end; ++sg_it) {
      //! MARK actual subgrid usage
      const HostEntityPointer host_entity_pointer = subGridPart.grid().getHostEntity<0>(*sg_it);
      const HostEntity& host_entity = *host_entity_pointer;

      const auto iend = fine_space.gridPart().iend(host_entity);
      for (auto iit = fine_space.gridPart().ibegin(host_entity); iit != iend; ++iit) {

        // if it is not a boundary element of the subgrid: continue;
        bool is_subgrid_boundary_face = true;
        if (iit->neighbor()) // if there is a neighbor entity
        {
          // check if the neighbor entity is in the subgrid
          const HostEntityPointer neighborHostEntityPointer = iit->outside();
          const HostEntity& neighborHostEntity = *neighborHostEntityPointer;
          //! MARK actual subgrid usage
          if (subGridPart.grid().contains<0>(neighborHostEntity))
            is_subgrid_boundary_face = false;
        }

        if (is_subgrid_boundary_face == false)
          continue;

        for (int c = 0; c < iit->geometry().corners(); ++c) {
          DomainType fine_corner = iit->geometry().corner(c);
          for (unsigned int cn = 0; cn < coarse_nodes_in_subgrid.size(); ++cn) {
            if (coarse_nodes_in_subgrid[cn] == fine_corner) {
              OrderedDomainType coarse_node_coordinates = fine_corner;
              int global_coarse_node_id = coordinates_to_global_coarse_node_id[coarse_node_coordinates];
              if (std::find(coarse_boundary_node_ids_in_subgrid[sg_id].begin(),
                            coarse_boundary_node_ids_in_subgrid[sg_id].end(),
                            global_coarse_node_id) == coarse_boundary_node_ids_in_subgrid[sg_id].end()) {
                coarse_boundary_node_ids_in_subgrid[sg_id].push_back(global_coarse_node_id);
              }
            }
          }
        }
      }
    }
  }

  for (unsigned int sg_id = 0; sg_id < number_of_subgrids; sg_id += 1) {
    for (unsigned int i_all = 0; i_all < coarse_node_ids_in_subgrid[sg_id].size(); ++i_all) {

      if (std::find(coarse_boundary_node_ids_in_subgrid[sg_id].begin(),
                    coarse_boundary_node_ids_in_subgrid[sg_id].end(),
                    coarse_node_ids_in_subgrid[sg_id][i_all]) == coarse_boundary_node_ids_in_subgrid[sg_id].end()) {
        coarse_interior_node_ids_in_subgrid[sg_id].push_back(coarse_node_ids_in_subgrid[sg_id][i_all]);
      }
    }
  }

  // Finalize: for each subgrid, store the vector of basis functions ids that correspond to interior coarse grid nodes
  // in the subgrid
  ids_basis_function_in_interior_subgrid.resize(number_of_subgrids);
  for (unsigned int sg_id = 0; sg_id < number_of_subgrids; sg_id += 1)
    for (unsigned int i = 0; i < coarse_interior_node_ids_in_subgrid[sg_id].size(); ++i)
      ids_basis_function_in_interior_subgrid[sg_id]
          .push_back(global_id_to_internal_id[coarse_interior_node_ids_in_subgrid[sg_id][i]]);
}
// ------------------------------------------------------------------------------------------------------

//! create standard coarse grid basis functions as discrete functions defined on the fine grid
// ------------------------------------------------------------------------------------
void
Elliptic_Rigorous_MsFEM_Solver::add_coarse_basis_contribution(MacroMicroGridSpecifier& specifier,
                                                              std::map<int, int>& global_id_to_internal_id,
                                                              MsFEMBasisFunctionType& msfem_basis_function_list) const {

  DSC_LOG_INFO << "Create standard coarse grid basis functions as discrete functions on the fine grid... ";

  typedef typename HostEntity::Codim<HostGrid::dimension>::EntityPointer HostNodePointer;
  typedef typename GridPart::IntersectionIteratorType HostIntersectionIterator;
  typedef typename DiscreteFunctionSpace::BasisFunctionSetType CoarseBaseFunctionSet;

  const HostGridLeafIndexSet& coarseGridLeafIndexSet = specifier.coarseSpace().gridPart().grid().leafIndexSet();

  for (HostgridIterator it = discreteFunctionSpace_.begin(); it != discreteFunctionSpace_.end(); ++it) {
    typedef typename HostEntity::Codim<0>::EntityPointer HostEntityPointer;
    HostEntityPointer coarse_father =
        Stuff::Grid::make_father(coarseGridLeafIndexSet, HostEntityPointer(*it), specifier.getLevelDifference());

    const CoarseBaseFunctionSet coarseBaseSet = specifier.coarseSpace().basisFunctionSet(*coarse_father);
    const auto numBaseFunctions = coarseBaseSet.size();

    const auto& lagrangepoint_set = specifier.coarseSpace().lagrangePointSet(*coarse_father);
    const auto& coarse_geometry = (*coarse_father).geometry();

    const auto number_of_points = lagrangepoint_set.nop();

    std::vector<RangeType> phi(numBaseFunctions);
    std::vector<std::size_t> indices;
    specifier.coarseSpace().mapper().map(*coarse_father, indices);

    //! TODO switch loops for more efficiency
    for (size_t loc_basis_number = 0; loc_basis_number < numBaseFunctions; ++loc_basis_number) {
      const int global_dof_number = indices[loc_basis_number];
      if (specifier.is_coarse_dirichlet_node(global_dof_number)) {
        continue;
      }

      const int global_interior_dof_number = global_id_to_internal_id[global_dof_number];

      // only implemented for 3 Lagrange Points, i.e. piecewise linear functions
      assert(number_of_points == 3);
      std::vector<RangeType> phi_i(number_of_points);
      std::vector<DomainType> corners(number_of_points);

      for (size_t loc_point = 0; loc_point < number_of_points; ++loc_point) {

        coarseBaseSet.evaluateAll(lagrangepoint_set.point(loc_point), phi);
        phi_i[loc_point] = phi[loc_basis_number];
        corners[loc_point] = coarse_geometry.global(lagrangepoint_set.point(loc_point));
      }

      // LinearLagrangeInterpolation2D should be eventually replaced by
      // LinearLagrangeFunction2D< DiscreteFunctionSpace > coarse_basis_interpolation
      LinearLagrangeInterpolation2D<DiscreteFunctionSpace> coarse_basis_interpolation(corners[0], phi_i[0], corners[1],
                                                                                      phi_i[1], corners[2], phi_i[2]);

      LocalFunction loc_coarse_basis_function =
          (msfem_basis_function_list[global_interior_dof_number])->localFunction(*it);

      const int number_of_nodes_in_fine_entity = it->count<HostGrid::dimension>();
      if (!(number_of_nodes_in_fine_entity == int(loc_coarse_basis_function.basisFunctionSet().size()))) {
        DSC_LOG_ERROR << "Error! Inconsistency in 'rigorous_msfem_solver.hh'." << std::endl;
      }

      for (int i = 0; i < number_of_nodes_in_fine_entity; i += 1) {
        const HostNodePointer node = it->subEntity<HostGrid::dimension>(i);

        const DomainType coordinates_of_node = node->geometry().corner(0);
        if (!(coordinates_of_node == it->geometry().corner(i))) {
          DSC_LOG_ERROR << "Error! Inconsistency in 'rigorous_msfem_solver.hh'." << std::endl;
        }

        RangeType coarse_value(0.0);
        coarse_basis_interpolation.evaluate(coordinates_of_node, coarse_value);

        // int global_index_node = gridPart.indexSet().index( *node );
        loc_coarse_basis_function[i] = coarse_value;
      }
    }
  }
  DSC_LOG_INFO << " done." << std::endl;
}

//! add corrector part to MsFEM basis functions
void Elliptic_Rigorous_MsFEM_Solver::add_corrector_contribution(
    MacroMicroGridSpecifier& specifier, std::map<int, int>& global_id_to_internal_id, SubGridListType& subgrid_list,
    MsFEMBasisFunctionType& msfem_basis_function_list) const {

  DSC_LOG_INFO << "Add global corrector to create MsFEM basis functions from standard FEM basis functions... ";

  const HostGridLeafIndexSet& coarseGridLeafIndexSet = specifier.coarseSpace().gridPart().grid().leafIndexSet();

  typedef typename DiscreteFunctionSpace::IteratorType CoarseIterator;
  typedef typename CoarseIterator::Entity CoarseEntity;
  typedef typename CoarseEntity::Geometry CoarseGeometry;

  typedef typename DiscreteFunctionSpace::BasisFunctionSetType CoarseBaseFunctionSet;

  std::vector<JacobianRangeType> gradient_Phi(specifier.coarseSpace().mapper().maxNumDofs());

  for (const CoarseEntity& coarse_grid_entity : specifier.coarseSpace()) {
    assert(coarse_grid_entity.partitionType() == InteriorEntity);

    const int global_index_entity = coarseGridLeafIndexSet.index(coarse_grid_entity);

    const CoarseBaseFunctionSet coarseBaseSet = specifier.coarseSpace().basisFunctionSet(coarse_grid_entity);
    const auto numBaseFunctions = coarseBaseSet.size();

    // the sub grid U(T) that belongs to the coarse_grid_entity T
    auto subGridPart = subgrid_list.gridPart(global_index_entity);

    const SubGridDiscreteFunctionSpaceType localDiscreteFunctionSpace(subGridPart);

    SubGridDiscreteFunctionType local_problem_solution_e0("Local problem Solution e_0", localDiscreteFunctionSpace);
    local_problem_solution_e0.clear();

    SubGridDiscreteFunctionType local_problem_solution_e1("Local problem Solution e_1", localDiscreteFunctionSpace);
    local_problem_solution_e1.clear();

    // --------- load local solutions -------
    // the file/place, where we saved the solutions of the cell problems
    const std::string local_solution_location =
        (boost::format("local_problems/_localProblemSolutions_%d") %
         specifier.coarseSpace().gridPart().grid().globalIdSet().id(coarse_grid_entity)).str();
    // reader for the cell problem data file:
    DiscreteFunctionReader discrete_function_reader(local_solution_location);
    // std::cout<< "... reading local problem solution " << global_index_entity << "/" << 0 << std::endl;
    discrete_function_reader.read(0, local_problem_solution_e0);
    // std::cout<< "... reading local problem solution " << global_index_entity << "/" << 1 << std::endl;
    discrete_function_reader.read(1, local_problem_solution_e1);

    // 1 point quadrature!! We only need the gradient of the base function,
    // which is constant on the whole entity.
    const CoarseQuadrature one_point_quadrature(coarse_grid_entity, 0);
    coarseBaseSet.jacobianAll(one_point_quadrature[0], gradient_Phi);

    std::vector<std::size_t> indices;
    specifier.coarseSpace().mapper().map(coarse_grid_entity, indices);
    for (unsigned int i = 0; i < numBaseFunctions; ++i) {
      int global_dof_number = indices[i];
      if (specifier.is_coarse_dirichlet_node(global_dof_number)) {
        continue;
      }

      int global_interior_dof_number = global_id_to_internal_id[global_dof_number];

      DiscreteFunction correction_on_U_T("correction_on_U_T", discreteFunctionSpace_);

      correction_on_U_T.clear();
      subgrid_to_hostrid_projection(local_problem_solution_e0, correction_on_U_T);
      correction_on_U_T *= gradient_Phi[i][0][0];
      (*(msfem_basis_function_list[global_interior_dof_number])) += correction_on_U_T;

      correction_on_U_T.clear();
      subgrid_to_hostrid_projection(local_problem_solution_e1, correction_on_U_T);
      correction_on_U_T *= gradient_Phi[i][0][1];
      (*(msfem_basis_function_list[global_interior_dof_number])) += correction_on_U_T;
    }
  }

  DSC_LOG_INFO << " done." << std::endl;
}

//! assemble global dirichlet corrector
void Elliptic_Rigorous_MsFEM_Solver::assemble_global_dirichlet_corrector(
    MacroMicroGridSpecifier& specifier, MsFEMTraits::SubGridListType& subgrid_list,
    Elliptic_Rigorous_MsFEM_Solver::DiscreteFunction& global_dirichlet_corrector) const {

  DSC_LOG_INFO << "Sum up local Dirichlet boundary correctors to create the global Dirichlet Corrector... ";

  global_dirichlet_corrector.clear();
  const HostGridLeafIndexSet& coarseGridLeafIndexSet = specifier.coarseSpace().gridPart().grid().leafIndexSet();

  typedef typename DiscreteFunctionSpace::IteratorType CoarseIterator;
  typedef typename CoarseIterator::Entity CoarseEntity;
  typedef typename CoarseEntity::Geometry CoarseGeometry;

  for (const CoarseEntity& coarse_grid_entity : specifier.coarseSpace()) {

    assert(coarse_grid_entity.partitionType() == InteriorEntity);

    const int global_index_entity = coarseGridLeafIndexSet.index(coarse_grid_entity);

    std::vector<std::size_t> indices;
    specifier.coarseSpace().mapper().map(coarse_grid_entity, indices);

    // is a Dirichlet boundary corrector available for this coarse entity (i.e. for the corresponding subgrid)
    bool boundary_corrector_available = false;
    for (const auto& intersection :
         Dune::Stuff::Common::intersectionRange(specifier.coarseSpace().gridPart(), coarse_grid_entity)) {

      // boundaryId 1 = Dirichlet face; boundaryId 2 = Neumann face;
      if (intersection.boundary() && (intersection.boundaryId() == 1)) {
        boundary_corrector_available = true;
      } else {
        const auto& lagrangePointSet = specifier.coarseSpace().lagrangePointSet(coarse_grid_entity);

        const int face = intersection.indexInInside();
        auto faceIterator = lagrangePointSet.beginSubEntity<faceCodim>(face);
        const auto faceEndIterator = lagrangePointSet.endSubEntity<faceCodim>(face);

        for (; faceIterator != faceEndIterator; ++faceIterator) {
          if (specifier.is_coarse_dirichlet_node(indices[*faceIterator])) {
            boundary_corrector_available = true;
            continue;
          }
        }
        continue;
      }
    }

    if (boundary_corrector_available) {
      // the sub grid U(T) that belongs to the coarse_grid_entity T
      auto subGridPart = subgrid_list.gridPart(global_index_entity);
      const SubGridDiscreteFunctionSpaceType localDiscreteFunctionSpace(subGridPart);

      SubGridDiscreteFunctionType local_dirichlet_corrector("Local Dirichlet corrector", localDiscreteFunctionSpace);
      local_dirichlet_corrector.clear();

      // --------- load local Dirichlet corrector -------
      // the file/place, where we saved the corrector
      const std::string local_solution_location =
          (boost::format("local_problems/_dirichletBoundaryCorrector_%d") %
           specifier.coarseSpace().gridPart().grid().globalIdSet().id(coarse_grid_entity)).str();

      // reader for the cell problem data file:
      DiscreteFunctionReader discrete_function_reader(local_solution_location);
      discrete_function_reader.read(0, local_dirichlet_corrector);

      DiscreteFunction local_dirichlet_corrector_aux("local_dirichlet_corrector_aux", discreteFunctionSpace_);
      local_dirichlet_corrector_aux.clear();
      subgrid_to_hostrid_projection(local_dirichlet_corrector, local_dirichlet_corrector_aux);
      global_dirichlet_corrector += local_dirichlet_corrector_aux;
    }
  }

  DSC_LOG_INFO << " done." << std::endl;
}

//! assemble global neumann corrector
void Elliptic_Rigorous_MsFEM_Solver::assemble_global_neumann_corrector(
    MacroMicroGridSpecifier& specifier, MsFEMTraits::SubGridListType& subgrid_list,
    Elliptic_Rigorous_MsFEM_Solver::DiscreteFunction& global_neumann_corrector) const {

  DSC_LOG_INFO << "Sum up local Neumann boundary correctors to create the global Neumann Corrector... ";

  global_neumann_corrector.clear();
  const HostGridLeafIndexSet& coarseGridLeafIndexSet = specifier.coarseSpace().gridPart().grid().leafIndexSet();

  typedef typename DiscreteFunctionSpace::IteratorType CoarseIterator;
  typedef typename CoarseIterator::Entity CoarseEntity;
  typedef typename CoarseEntity::Geometry CoarseGeometry;

  for (const CoarseEntity& coarse_grid_entity : specifier.coarseSpace()) {

    assert(coarse_grid_entity.partitionType() == InteriorEntity);

    const int global_index_entity = coarseGridLeafIndexSet.index(coarse_grid_entity);

    // is a neumann boundary corrector available for this coarse entity (i.e. for the corresponding subgrid)
    bool boundary_corrector_available = false;
    for (const auto& intersection :
         Dune::Stuff::Common::intersectionRange(specifier.coarseSpace().gridPart(), coarse_grid_entity)) {

      // boundaryId 1 = Dirichlet face; boundaryId 2 = Neumann face;
      if (intersection.boundary() && (intersection.boundaryId() == 2))
        boundary_corrector_available = true;
    }

    if (boundary_corrector_available) {
      // the sub grid U(T) that belongs to the coarse_grid_entity T
      auto subGridPart = subgrid_list.gridPart(global_index_entity);
      const SubGridDiscreteFunctionSpaceType localDiscreteFunctionSpace(subGridPart);

      SubGridDiscreteFunctionType local_neumann_corrector("Local Neumann corrector", localDiscreteFunctionSpace);
      local_neumann_corrector.clear();

      // --------- load local neumann corrector -------
      // the file/place, where we saved the corrector
      const std::string local_solution_location =
          (boost::format("local_problems/_neumannBoundaryCorrector_%d") %
           specifier.coarseSpace().gridPart().grid().globalIdSet().id(coarse_grid_entity)).str();

      // reader for the cell problem data file:
      DiscreteFunctionReader discrete_function_reader(local_solution_location);
      discrete_function_reader.read(0, local_neumann_corrector);

      DiscreteFunction local_neumann_corrector_aux("local_neumann_corrector_aux", discreteFunctionSpace_);
      local_neumann_corrector_aux.clear();
      subgrid_to_hostrid_projection(local_neumann_corrector, local_neumann_corrector_aux);
      global_neumann_corrector += local_neumann_corrector_aux;
    }
  }

  DSC_LOG_INFO << " done." << std::endl;
}

void Elliptic_Rigorous_MsFEM_Solver::solve(
    const CommonTraits::DiffusionType& diffusion_op, const CommonTraits::FirstSourceType& f,
    const CommonTraits::DiscreteFunctionType& dirichlet_extension, const CommonTraits::NeumannBCType& neumann_bc,
    // number of layers per coarse grid entity T:  U(T) is created by enrichting T with
    // n(T)-layers.
    MacroMicroGridSpecifier& specifier, SubGridListType& subgrid_list, DiscreteFunction& coarse_scale_part,
    DiscreteFunction& fine_scale_part, DiscreteFunction& solution) const {
  DSC::Profiler::ScopedTiming st("msfem.Elliptic_Rigorous_MsFEM_Solver.solve_dirichlet_zero");

  specifier.setOversamplingStrategy(3); // for rigorous MsFEM!

  DiscreteFunctionSpace& coarse_space = specifier.coarseSpace();
  DiscreteFunctionSpace& fine_space = specifier.fineSpace();

  const int number_of_relevant_coarse_nodes = coarse_space.size() - specifier.get_number_of_coarse_dirichlet_nodes();

  // mapper: global_id_of_node -> new_id_of_node
  // ('new' means that we only count the internal nodes and the non-Dirichlet boundary nodes,
  //  Dirichlet boundary nodes do not receive an id)
  std::map<int, int> global_id_to_internal_id;
  for (int internal_id = 0, global_id = 0; global_id < coarse_space.size(); ++global_id) {
    if (!specifier.is_coarse_dirichlet_node(global_id)) {
      global_id_to_internal_id[global_id] = internal_id;
      ++internal_id;
    }
  }

  MsFEMBasisFunctionType msfem_basis_function;
  for (int internal_id = 0; internal_id < number_of_relevant_coarse_nodes; ++internal_id) {
    msfem_basis_function.emplace_back(new DiscreteFunction("MsFEM basis function", fine_space));
    msfem_basis_function[internal_id]->clear();
  }

  // determine (\int_\Omega coarse standard basis function )^{-1}, which will be a weight in the weighted
  // clement interpolation
  // ------------------------------------------------------------------------------------------------------
  // coefficients in the matrix that describes the weighted Clement interpolation,
  // i.e. coff[c] = (\int_{\Omega} \Phi_j)^{-1}
  std::vector<double> coff(number_of_relevant_coarse_nodes, 0.0);

  for (CoarsegridIterator it = coarse_space.begin(); it != coarse_space.end(); ++it) {

    HostEntity& entity = *it;

    assert(entity.partitionType() == InteriorEntity);

    std::vector<RangeType> phi(coarse_space.mapper().maxNumDofs());

    // get base function set
    const auto& coarse_baseSet = coarse_space.basisFunctionSet(entity);
    const auto numBaseFunctions = coarse_baseSet.size();

    // create quadrature of appropriate order
    CoarseQuadrature quadrature(entity, 2 * DiscreteFunctionSpace::polynomialOrder + 2);

    // loop over all quadrature points
    const auto numQuadraturePoints = quadrature.nop();
    for (size_t quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint) {
      const typename CoarseQuadrature::CoordinateType& local_point = quadrature.point(quadraturePoint);

      const double weight = quadrature.weight(quadraturePoint) * entity.geometry().integrationElement(local_point);

      coarse_baseSet.evaluateAll(quadrature[quadraturePoint], phi);

      std::vector<std::size_t> indices;
      coarse_space.mapper().map(entity, indices);
      for (unsigned int i = 0; i < numBaseFunctions; ++i) {
        const int global_dof_number = indices[i];
        if (!specifier.is_coarse_dirichlet_node(global_dof_number))
          coff[global_id_to_internal_id[global_dof_number]] += weight * phi[i];
      }
    }
  }

  for (size_t c = 0; c < coff.size(); ++c) {
    if (coff[c] != 0.0)
      coff[c] = 1.0 / coff[c];
  }
  // ------------------------------------------------------------------------------------------------------

  //! Determine the support of each ms basis function and the intersection domain of two ms basis functions
  // save the corresponding entity seeds in support_of_ms_basis_func_intersection
  // ------------------------------------------------------------------------------------------------------

  // the support of the interaction of two ms basis functions
  std::vector<std::vector<std::vector<FineGridEntitySeed>>> support_of_ms_basis_func_intersection;

  // support of the global Dirichlet and the global Neumann boundary corrector
  std::vector<FineGridEntitySeed> support_global_dirichlet_corrector;
  std::vector<FineGridEntitySeed> support_global_neumann_corrector;

  // relevant constellations of two ms basis functions (i.e. when they have a non-empty intersection of their supports)
  std::vector<std::tuple<unsigned int, unsigned int>> relevant_constellations;
  // we only store the tuples relevant_constellations[i][j] for 'j<=i' the rest is obtained by symmetry
  // the reason for storing only the values for 'j<=i' is that we can use (if given) the symmetry of the diffusion
  // matrix,
  // in the sense that we only compute the entries of the stiffness matrix for 'j<=i' and then symmetrize the matrix

  // map the coordinates of a coarse node to its global (coarse node) index
  std::map<OrderedDomainType, int> coordinates_to_global_coarse_node_id;

  support_of_ms_basis_func_intersection.resize(number_of_relevant_coarse_nodes);
  for (int k = 0; k < number_of_relevant_coarse_nodes; ++k)
    support_of_ms_basis_func_intersection[k].resize(number_of_relevant_coarse_nodes);

  // for each subgrid, determine all ms basis functions that were constructed using the corresponding subgrid corrector
  // for each subgrid id, we save the vector of the (internal) id's of the corresponding ms basis functions
  std::vector<std::vector<int>> subgrid_id_to_ms_basis_func_ids;
  subgrid_id_to_ms_basis_func_ids.resize(coarse_space.gridPart().grid().size(0));
  for (CoarsegridIterator it = coarse_space.begin(); it != coarse_space.end(); ++it) {
    const CoarseGridLeafIndexSet& coarseGridLeafIndexSet = coarse_space.gridPart().grid().leafIndexSet();
    int subgrid_id = coarseGridLeafIndexSet.index(*it);

    auto intersection_it = coarse_space.gridPart().ibegin(*it);
    const auto endiit = coarse_space.gridPart().iend(*it);
    for (; intersection_it != endiit; ++intersection_it) {

      const auto& lagrangePointSet = coarse_space.lagrangePointSet(*it);

      const int face = (*intersection_it).indexInInside();

      std::vector<std::size_t> indices;
      coarse_space.mapper().map(*it, indices);

      auto faceIterator = lagrangePointSet.beginSubEntity<faceCodim>(face);
      const auto faceEndIterator = lagrangePointSet.endSubEntity<faceCodim>(face);
      for (; faceIterator != faceEndIterator; ++faceIterator) {
        const int global_id_node = indices[*faceIterator];

        // create map entry: ' global coord of node <--> global_id_node '
        OrderedDomainType coord = (it->geometry()).global(lagrangePointSet.point(*faceIterator));
        coordinates_to_global_coarse_node_id[coord] = global_id_node;

        if (specifier.is_coarse_dirichlet_node(global_id_node))
          continue;

        const int internal_id_node = global_id_to_internal_id[global_id_node];
        subgrid_id_to_ms_basis_func_ids[subgrid_id].push_back(internal_id_node);
      }
    }
  }

  const HostGridLeafIndexSet& hostGridLeafIndexSet = fine_space.gridPart().grid().leafIndexSet();
  for (HostgridIterator it = fine_space.begin(); it != fine_space.end(); ++it) {

    int fine_entity_id = hostGridLeafIndexSet.index(*it);

    bool add_to_support_dirichlet_corrector = false;
    bool add_to_support_neumann_corrector = false;

    // IDs of the ms basis functions that containt the fine grid entity 'it'
    std::vector<int> ms_basis_funcs_that_contain_entity;

    // first iterate of the subgrids that contain the fine grid entity
    for (unsigned int m = 0; m < subgrid_list.getSubgridIDs_that_contain_entity(fine_entity_id).size(); ++m) {
      int subgrid_id = subgrid_list.getSubgridIDs_that_contain_entity(fine_entity_id)[m];

      HostEntityPointer coarse_it = coarse_space.grid().entityPointer(subgrid_list.get_coarse_entity_seed(subgrid_id));
      HostEntity& coarse_entity = *coarse_it;
      for (const auto& coarse_intersection :
           Dune::Stuff::Common::intersectionRange(coarse_space.gridPart(), coarse_entity)) {

        // boundaryId 1 = Dirichlet face; boundaryId 2 = Neumann face;
        if (coarse_intersection.boundary() && (coarse_intersection.boundaryId() == 1))
          add_to_support_dirichlet_corrector = true;
        if (coarse_intersection.boundary() && (coarse_intersection.boundaryId() == 2))
          add_to_support_neumann_corrector = true;
      }

      // now iterate over the ms basis functions that were assembled using a corrector that belongs to the current
      // subgrid
      for (unsigned int l = 0; l < subgrid_id_to_ms_basis_func_ids[subgrid_id].size(); ++l) {
        int ms_basis_func_id = subgrid_id_to_ms_basis_func_ids[subgrid_id][l];
        if (std::find(ms_basis_funcs_that_contain_entity.begin(), ms_basis_funcs_that_contain_entity.end(),
                      ms_basis_func_id) == ms_basis_funcs_that_contain_entity.end()) {
          ms_basis_funcs_that_contain_entity.push_back(ms_basis_func_id);
        }
      }
    }

    for (unsigned int mid1 = 0; mid1 < ms_basis_funcs_that_contain_entity.size(); ++mid1) {
      int ms_basis_id_1 = ms_basis_funcs_that_contain_entity[mid1];
      for (unsigned int mid2 = 0; mid2 < ms_basis_funcs_that_contain_entity.size(); ++mid2) {
        int ms_basis_id_2 = ms_basis_funcs_that_contain_entity[mid2];
        // check if the tuple was already added to the relevant_constellations-vector
        if ((support_of_ms_basis_func_intersection[ms_basis_id_1][ms_basis_id_2].size() == 0) &&
            (ms_basis_id_2 <= ms_basis_id_1)) {
          std::tuple<unsigned int, unsigned int> tup = make_tuple(ms_basis_id_1, ms_basis_id_2);
          relevant_constellations.push_back(tup);
        }
        // add the seed to the support of the insection between ms basis ms_basis_id_1 and ms_basis_id_2
        support_of_ms_basis_func_intersection[ms_basis_id_1][ms_basis_id_2].push_back((*it).seed());
      }
    }

    if (add_to_support_dirichlet_corrector == true)
      support_global_dirichlet_corrector.push_back((*it).seed());

    if (add_to_support_neumann_corrector == true)
      support_global_neumann_corrector.push_back((*it).seed());
  }
  // ------------------------------------------------------------------------------------------------------

  MsFEMBasisFunctionType standard_basis_function;
  for (int internal_id = 0; internal_id < number_of_relevant_coarse_nodes; internal_id += 1) {
    standard_basis_function.emplace_back(new DiscreteFunction("Standard basis function", fine_space));
    standard_basis_function[internal_id]->clear();
  }

  add_coarse_basis_contribution(specifier, global_id_to_internal_id, msfem_basis_function);
  add_coarse_basis_contribution(specifier, global_id_to_internal_id, standard_basis_function);

  // for each subgrid, store the vector of basis functions ids that correspond to coarse grid nodes
  // in the subgrid WITHOUT Dirichlet boundary nodes of Omega
  std::vector<std::vector<int>> ids_basis_functions_in_subgrid;
  // for each subgrid, store the vector of basis functions ids that correspond to interior coarse grid nodes in the
  // subgrid
  std::vector<std::vector<int>> ids_basis_functions_in_interior_subgrid;
  // for each subgrid, store the vector of basis functions ids where the support of the bassis function intersects with
  // the subgrid
  std::vector<std::vector<int>> ids_basis_functions_in_extended_subgrid;
  assemble_interior_basis_ids(specifier, subgrid_list, global_id_to_internal_id, coordinates_to_global_coarse_node_id,
                              ids_basis_functions_in_extended_subgrid, ids_basis_functions_in_subgrid,
                              ids_basis_functions_in_interior_subgrid);

  //! assemble all local problems (within constructor!)
  MsFEMLocalProblemSolver loc_prob_solver(specifier.fineSpace(), specifier, subgrid_list,
                                          // ids_basis_functions_in_interior_subgrid,
                                          // ids_basis_functions_in_subgrid,
                                          ids_basis_functions_in_extended_subgrid, coff, diffusion_op,
                                          standard_basis_function, global_id_to_internal_id, neumann_bc,
                                          dirichlet_extension);
  loc_prob_solver.assemble_all(/*silence=*/false);

  // define the discrete (elliptic) operator that describes our problem
  add_corrector_contribution(specifier, global_id_to_internal_id, subgrid_list, msfem_basis_function);

  DiscreteFunction global_dirichlet_corrector("Global Dirichlet Corrector", fine_space);
  DiscreteFunction global_neumann_corrector("Global Neumann Corrector", fine_space);
  // assemble the global boundary correctors
  assemble_global_dirichlet_corrector(specifier, subgrid_list, global_dirichlet_corrector);
  assemble_global_neumann_corrector(specifier, subgrid_list, global_neumann_corrector);

  // just for VTK output for the basis function correctors
  /*
  MsFEMBasisFunctionType corrector_basis_function;
  for (int internal_id = 0; internal_id < number_of_relevant_coarse_nodes; internal_id += 1 )
   {
    corrector_basis_function.emplace_back(new DiscreteFunction("Corrector basis function", fine_space));
    corrector_basis_function[internal_id]->clear();
   }

  add_corrector_contribution( specifier, global_id_to_internal_id, subgrid_list, corrector_basis_function );
  vtk_output( corrector_basis_function, "corrector_basis_function" );
  */

  if (DSC_CONFIG_GET("rigorous_msfem.msfem_basis_vtk_output", 0)) {
    vtk_output(msfem_basis_function);
    vtk_output(standard_basis_function, "standard_basis_function");
  }

  VectorType solution_vector(number_of_relevant_coarse_nodes);
  for (size_t col = 0; col != solution_vector.N(); ++col)
    solution_vector[col] = 0.0;

  double tol = DSC_CONFIG_GET("rigorous_msfem.macro_solver_tolerance", 1e-10);
  int num_iterations = DSC_CONFIG_GET("rigorous_msfem.macro_solver_iterations", 10000);

  // linear version
  if (DSC_CONFIG_GET("problem.linear", true)) {

    DSC_LOG_INFO << "Start assembling the stiffness matrix of the global problems.." << std::endl;
    Dune::Timer assembleTimer;

    //! (stiffness) matrix
    MatrixType system_matrix(number_of_relevant_coarse_nodes, number_of_relevant_coarse_nodes);

    if (DSC_CONFIG_GET("lod.petrov_galerkin", true)) {
      assemble_matrix(diffusion_op, msfem_basis_function, standard_basis_function,
                      support_of_ms_basis_func_intersection, relevant_constellations, system_matrix);
    } else {
      assemble_matrix(diffusion_op, msfem_basis_function, msfem_basis_function, support_of_ms_basis_func_intersection,
                      relevant_constellations, system_matrix);
    }
    // NOTE: in the case that we use the Petrov Galerkin version of the method 'support_of_ms_basis_func_intersection'
    // is not yet optimally assembled (it is a little larger as required, since we still determine the intersection of
    // two ms basis functions,
    // whereas the support of the classical basis function is typically smaller). It is correct, but not optimal!

    DSC_LOG_INFO << ".. assembling of the stiffness matrix done." << std::endl;
    DSC_LOG_INFO << "Time to assemble Rigorous MsFEM stiffness matrix: " << assembleTimer.elapsed() << "s" << std::endl;

    VectorType rhs(number_of_relevant_coarse_nodes);
    if (DSC_CONFIG_GET("lod.petrov_galerkin", true)) {
      assemble_rhs(f, diffusion_op, dirichlet_extension, neumann_bc, global_dirichlet_corrector,
                   global_neumann_corrector,
                   //                     support_global_dirichlet_corrector, support_global_neumann_corrector,
                   standard_basis_function, support_of_ms_basis_func_intersection, rhs);
    } else {
      assemble_rhs(f, diffusion_op, dirichlet_extension, neumann_bc, global_dirichlet_corrector,
                   global_neumann_corrector,
                   //                     support_global_dirichlet_corrector, support_global_neumann_corrector,
                   msfem_basis_function, support_of_ms_basis_func_intersection, rhs);
    }

    MatrixOperatorType matrix_op(system_matrix);
    PreconditionerType preconditioner(system_matrix, 100, 0.9);

    Dune::InverseOperatorResult result_data;

    if (DSC_CONFIG_GET("lod.algebraic_solver", "bcgs") == "bcgs") {
      typedef Dune::BiCGSTABSolver<VectorType> SolverType;

      SolverType solver(matrix_op, preconditioner, tol, num_iterations, true);
      solver.apply(solution_vector, rhs, result_data);
    } else if (DSC_CONFIG_GET("lod.algebraic_solver", "bcgs") == "cg") {
      if (DSC_CONFIG_GET("lod.petrov_galerkin", true)) {
        std::cout << "Warning! Key 'lod.petrov_galerkin' is set true (leading to an unsymetric system matrix)."
                  << "This is incompatible with the key "
                  << "'lod.algebraic_solver' which is set to 'cg' (only for symetric systems). "
                  << "Automatically switched the key to 'bcgs'." << std::endl;
        typedef Dune::BiCGSTABSolver<VectorType> SolverType;
        SolverType solver(matrix_op, preconditioner, tol, num_iterations, true);
        solver.apply(solution_vector, rhs, result_data);
      } else {
        typedef Dune::CGSolver<VectorType> SolverType;
        SolverType solver(matrix_op, preconditioner, tol, num_iterations, true);
        solver.apply(solution_vector, rhs, result_data);
      }
    } else {
      std::cout << "Unknown key for 'lod.algebraic_solver'." << std::endl;
      abort();
    }

    coarse_scale_part.clear();
    for (int internal_id = 0; internal_id < number_of_relevant_coarse_nodes; internal_id += 1) {
      DiscreteFunction aux("auxilliary function", fine_space);
      aux.clear();
      aux += *(standard_basis_function[internal_id]);
      aux *= solution_vector[internal_id];
      coarse_scale_part += aux;
    }

    solution.clear();
    for (int internal_id = 0; internal_id < number_of_relevant_coarse_nodes; internal_id += 1) {
      DiscreteFunction aux("auxilliary function", fine_space);
      aux.clear();
      aux += *(msfem_basis_function[internal_id]);
      aux *= solution_vector[internal_id];
      solution += aux;
    }
    solution += global_dirichlet_corrector;
    solution -= global_neumann_corrector;

    fine_scale_part.assign(solution);
    fine_scale_part -= coarse_scale_part;
  } else // if nonlinear
  {

    //! nonlinear version:
    // montone nonlinear problem (nonlinearity only in the lower order terms)
    // solve: - div( A \grad u ) + F( x, u , \grad u ) = f
    // implementation not yet efficient, since the re-assemblation of the system matrices is quite expensive
    // (always need to got down to the micro-scale for evaluations)
    // nonlinearity = F

    DSC_LOG_INFO << std::endl << "Starting Newton iterations." << std::endl;

    // Dune::Multiscale::Problem::Eleven::Nonlinearity nonlinear_term;
    auto nonlinear_term_ptr = Dune::Multiscale::Problem::getLowerOrderTerm();
    const auto& nonlinear_term = *nonlinear_term_ptr;
    solution.clear();

    VectorType newton_solution_vector(number_of_relevant_coarse_nodes);
    VectorType newton_step_rhs(number_of_relevant_coarse_nodes);
    for (size_t col = 0; col != newton_solution_vector.N(); ++col) {
      newton_solution_vector[col] = 0.0;
      newton_step_rhs[col] = 0.0;
    }

    double previous_newton_error = 10000.0;

    bool first_cycle = true;

    int iteration = 0;
    bool stop_newton_cycle = false;
    while (!stop_newton_cycle) {

      DSC_LOG_INFO << std::endl << "Newton iteration " << iteration << std::endl;
      double newton_tolerance = 1e-10;
      double damping_parameter = 1.0;

      //! Newton (stiffness) matrix
      MatrixType newton_system_matrix(number_of_relevant_coarse_nodes, number_of_relevant_coarse_nodes);

      for (size_t row = 0; row != newton_system_matrix.N(); ++row)
        for (size_t col = 0; col != newton_system_matrix.M(); ++col)
          newton_system_matrix[row][col] = 0.0;

      VectorType newton_step_solution_vector(number_of_relevant_coarse_nodes);
      for (size_t col = 0; col != solution_vector.N(); ++col)
        newton_step_solution_vector[col] = 0.0;

      if (DSC_CONFIG_GET("lod.petrov_galerkin", true)) {
        std::cout << "lod.petrov_galerkin not implemented for nonlinear problems!" << std::endl;
        abort();
      }

      for (unsigned int t = 0; t < relevant_constellations.size(); ++t) {
        unsigned int row = get<0>(relevant_constellations[t]);
        unsigned int col = get<1>(relevant_constellations[t]);

        int switch_row_col = -1;
        if (row == col) {
          switch_row_col = 0;
        }

        while (switch_row_col < 1) {

          int polOrder = 2 * DiscreteFunctionSpace::polynomialOrder + 2;
          for (int it_id = 0; it_id < support_of_ms_basis_func_intersection[row][col].size(); ++it_id) {

            HostEntityPointer it =
                discreteFunctionSpace_.grid().entityPointer(support_of_ms_basis_func_intersection[row][col][it_id]);

            LocalFunction loc_func_1 = msfem_basis_function[row]->localFunction(*it);
            LocalFunction loc_func_2 = msfem_basis_function[col]->localFunction(*it);
            LocalFunction local_dirichlet_extension = dirichlet_extension.localFunction(*it);
            LocalFunction glob_dirichlet_corrector_localized = global_dirichlet_corrector.localFunction(*it);
            LocalFunction glob_neumann_corrector_localized = global_neumann_corrector.localFunction(*it);

            const auto& geometry = (*it).geometry();

            if ((row == col) && (first_cycle == true)) {
              for (const auto& intersection :
                   Dune::Stuff::Common::intersectionRange(discreteFunctionSpace_.gridPart(), (*it))) {
                if (!intersection.boundary())
                  continue;
                // boundaryId 1 = Dirichlet face; boundaryId 2 = Neumann face;
                if (intersection.boundary() && (intersection.boundaryId() != 2))
                  continue;

                const HostFaceQuadrature faceQuadrature(discreteFunctionSpace_.gridPart(), intersection, polOrder,
                                                        HostFaceQuadrature::INSIDE);
                const int numFaceQuadraturePoints = faceQuadrature.nop();

                static const int faceCodim = 1;
                for (int faceQuadraturePoint = 0; faceQuadraturePoint < numFaceQuadraturePoints;
                     ++faceQuadraturePoint) {
                  RangeType func_in_x;
                  loc_func_1.evaluate(faceQuadrature[faceQuadraturePoint], func_in_x);

                  const auto local_point_entity = faceQuadrature.point(faceQuadraturePoint);
                  const auto global_point = geometry.global(local_point_entity);
                  const auto local_point_face = intersection.geometry().local(global_point);

                  RangeType neumann_value(0.0);
                  neumann_bc.evaluate(global_point, neumann_value);

                  const double face_weight = intersection.geometry().integrationElement(local_point_face) *
                                             faceQuadrature.weight(faceQuadraturePoint);

                  newton_step_rhs[col] += face_weight * (func_in_x * neumann_value);
                }
              }
            }

            const Fem::CachingQuadrature<GridPart, 0> quadrature(*it, polOrder);
            const int numQuadraturePoints = quadrature.nop();
            for (int quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint) {

              DomainType global_point = geometry.global(quadrature.point(quadraturePoint));

              // weight
              double weight = geometry.integrationElement(quadrature.point(quadraturePoint));
              weight *= quadrature.weight(quadraturePoint);

              // gradients of func1 and func2
              JacobianRangeType grad_func_1, grad_func_2, grad_previous_solution;
              loc_func_1.jacobian(quadrature[quadraturePoint], grad_func_1);
              loc_func_2.jacobian(quadrature[quadraturePoint], grad_func_2);

              RangeType value_func_1, value_func_2, value_previous_solution;
              loc_func_1.evaluate(quadrature[quadraturePoint], value_func_1);
              loc_func_2.evaluate(quadrature[quadraturePoint], value_func_2);

              solution.localFunction(*it).jacobian(quadrature[quadraturePoint], grad_previous_solution);
              solution.localFunction(*it).evaluate(quadrature[quadraturePoint], value_previous_solution);

              // A \nabla func1
              JacobianRangeType diffusive_flux(0.0);
              diffusion_op.diffusiveFlux(global_point, grad_func_1, diffusive_flux);

              RangeType value_local_dirichlet_extension;
              local_dirichlet_extension.evaluate(quadrature[quadraturePoint], value_local_dirichlet_extension);

              JacobianRangeType grad_local_dirichlet_extension;
              local_dirichlet_extension.jacobian(quadrature[quadraturePoint], grad_local_dirichlet_extension);

              JacobianRangeType grad_global_dirichlet_corrector;
              glob_dirichlet_corrector_localized.jacobian(quadrature[quadraturePoint], grad_global_dirichlet_corrector);

              RangeType value_global_dirichlet_corrector;
              glob_dirichlet_corrector_localized.evaluate(quadrature[quadraturePoint],
                                                          value_global_dirichlet_corrector);

              JacobianRangeType grad_global_neumann_corrector;
              glob_neumann_corrector_localized.jacobian(quadrature[quadraturePoint], grad_global_neumann_corrector);

              RangeType value_global_neumann_corrector;
              glob_neumann_corrector_localized.evaluate(quadrature[quadraturePoint], value_global_neumann_corrector);

              RangeType position = value_previous_solution + value_local_dirichlet_extension +
                                   value_global_dirichlet_corrector - value_global_neumann_corrector;
              JacobianRangeType direction(0.0);
              direction[0] = grad_previous_solution[0] + grad_local_dirichlet_extension[0] +
                             grad_global_dirichlet_corrector[0] - grad_global_neumann_corrector[0];

              RangeType value_F_x(0.0), value_derivative_1_F_x(0.0);
              JacobianRangeType value_derivative_2_F_x(0.0);
              nonlinear_term.evaluate(global_point, position, direction, value_F_x);
              nonlinear_term.position_derivative(global_point, position, direction, value_derivative_1_F_x);
              nonlinear_term.direction_derivative(global_point, position, direction, value_derivative_2_F_x);

              newton_system_matrix[row][col] += weight * (diffusive_flux[0] * grad_func_2[0]);
              newton_system_matrix[row][col] += weight * value_derivative_1_F_x * value_func_1 * value_func_2;
              newton_system_matrix[row][col] += weight * (value_derivative_2_F_x[0] * grad_func_1[0]) * value_func_2;

              if ((row == col) && (first_cycle == true)) {

                RangeType f_x(0.0);
                f.evaluate(global_point, f_x);

                JacobianRangeType diffusive_flux_in_grad_previous_solution_direction(0.0);
                diffusion_op.diffusiveFlux(global_point, direction, diffusive_flux_in_grad_previous_solution_direction);

                newton_step_rhs[col] += weight * (value_func_1 * f_x);
                newton_step_rhs[col] -=
                    weight * (diffusive_flux_in_grad_previous_solution_direction[0] * grad_func_1[0]);
                newton_step_rhs[col] -= weight * (value_func_1 * value_F_x);
              }
            }
          }

          unsigned int row_copy = row;
          row = col;
          col = row_copy;
          switch_row_col += 1;
        }
      }

      VectorType copy_newton_step_rhs(number_of_relevant_coarse_nodes);
      for (size_t col = 0; col != newton_step_rhs.N(); ++col)
        copy_newton_step_rhs[col] = newton_step_rhs[col];

      MatrixOperatorType newton_matrix_op(newton_system_matrix);
      PreconditionerType newton_preconditioner(newton_system_matrix, 100, 0.9);
      Dune::InverseOperatorResult newton_result_data;

      Dune::BiCGSTABSolver<VectorType> newton_solver(newton_matrix_op, newton_preconditioner, tol, num_iterations,
                                                     true);
      newton_solver.apply(newton_step_solution_vector, newton_step_rhs,
                          newton_result_data); // this changes newton_step_rhs!

      for (size_t col = 0; col != newton_step_rhs.N(); ++col)
        newton_step_rhs[col] = copy_newton_step_rhs[col];

      // assemble copys and check if the error is small and if damping is required
      double newton_error = 10000.0;
      if (first_cycle == true) {
        newton_error = 0.0;
        for (size_t col = 0; col != newton_solution_vector.N(); ++col) {
          newton_error += (copy_newton_step_rhs[col] * copy_newton_step_rhs[col]);
        }
        newton_error = sqrt(newton_error);
        previous_newton_error = newton_error;
        std::cout << "Newton error = " << newton_error << std::endl;
      }

      // while-loop not visited in the first Newton iteration cycle
      while (newton_error > previous_newton_error) {

        DiscreteFunction pre_solution("pre_solution", fine_space);
        pre_solution.clear();

        VectorType pre_newton_solution_vector(number_of_relevant_coarse_nodes);

        for (size_t col = 0; col != newton_solution_vector.N(); ++col) {
          pre_newton_solution_vector[col] =
              newton_solution_vector[col] + (damping_parameter * newton_step_solution_vector[col]);
          newton_step_rhs[col] = 0.0;
        }

        for (int internal_id = 0; internal_id < number_of_relevant_coarse_nodes; internal_id += 1) {
          DiscreteFunction aux("auxilliary function", fine_space);
          aux.clear();
          aux += *(msfem_basis_function[internal_id]);
          aux *= pre_newton_solution_vector[internal_id];
          pre_solution += aux;
        }

        // evaluate G(solution) = newton_step_rhs
        for (unsigned int t = 0; t < relevant_constellations.size(); ++t) {
          unsigned int row = get<0>(relevant_constellations[t]);
          unsigned int col = get<1>(relevant_constellations[t]);

          if (row != col) {
            continue;
          }

          int polOrder = 2 * DiscreteFunctionSpace::polynomialOrder + 2;
          for (int it_id = 0; it_id < support_of_ms_basis_func_intersection[row][row].size(); ++it_id) {

            HostEntityPointer it =
                discreteFunctionSpace_.grid().entityPointer(support_of_ms_basis_func_intersection[row][col][it_id]);

            LocalFunction loc_func = msfem_basis_function[row]->localFunction(*it);
            LocalFunction local_dirichlet_extension = dirichlet_extension.localFunction(*it);
            LocalFunction glob_dirichlet_corrector_localized = global_dirichlet_corrector.localFunction(*it);
            LocalFunction glob_neumann_corrector_localized = global_neumann_corrector.localFunction(*it);

            const auto& geometry = (*it).geometry();

            for (const auto& intersection :
                 Dune::Stuff::Common::intersectionRange(discreteFunctionSpace_.gridPart(), (*it))) {
              if (!intersection.boundary())
                continue;
              // boundaryId 1 = Dirichlet face; boundaryId 2 = Neumann face;
              if (intersection.boundary() && (intersection.boundaryId() != 2))
                continue;

              const HostFaceQuadrature faceQuadrature(discreteFunctionSpace_.gridPart(), intersection, polOrder,
                                                      HostFaceQuadrature::INSIDE);
              const int numFaceQuadraturePoints = faceQuadrature.nop();

              static const int faceCodim = 1;
              for (int faceQuadraturePoint = 0; faceQuadraturePoint < numFaceQuadraturePoints; ++faceQuadraturePoint) {
                RangeType func_in_x;
                loc_func.evaluate(faceQuadrature[faceQuadraturePoint], func_in_x);

                const auto local_point_entity = faceQuadrature.point(faceQuadraturePoint);
                const auto global_point = geometry.global(local_point_entity);
                const auto local_point_face = intersection.geometry().local(global_point);

                RangeType neumann_value(0.0);
                neumann_bc.evaluate(global_point, neumann_value);

                const double face_weight = intersection.geometry().integrationElement(local_point_face) *
                                           faceQuadrature.weight(faceQuadraturePoint);

                newton_step_rhs[col] += face_weight * (func_in_x * neumann_value);
              }
            }

            const Fem::CachingQuadrature<GridPart, 0> quadrature(*it, polOrder);
            const int numQuadraturePoints = quadrature.nop();
            for (int quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint) {

              DomainType global_point = geometry.global(quadrature.point(quadraturePoint));

              // weight
              double weight = geometry.integrationElement(quadrature.point(quadraturePoint));
              weight *= quadrature.weight(quadraturePoint);

              // gradients of func1 and func2
              JacobianRangeType grad_func, grad_previous_solution;
              loc_func.jacobian(quadrature[quadraturePoint], grad_func);

              RangeType value_func, value_previous_solution;
              loc_func.evaluate(quadrature[quadraturePoint], value_func);

              pre_solution.localFunction(*it).jacobian(quadrature[quadraturePoint], grad_previous_solution);
              pre_solution.localFunction(*it).evaluate(quadrature[quadraturePoint], value_previous_solution);

              RangeType value_local_dirichlet_extension;
              local_dirichlet_extension.evaluate(quadrature[quadraturePoint], value_local_dirichlet_extension);

              JacobianRangeType grad_local_dirichlet_extension;
              local_dirichlet_extension.jacobian(quadrature[quadraturePoint], grad_local_dirichlet_extension);

              JacobianRangeType grad_global_dirichlet_corrector;
              glob_dirichlet_corrector_localized.jacobian(quadrature[quadraturePoint], grad_global_dirichlet_corrector);

              RangeType value_global_dirichlet_corrector;
              glob_dirichlet_corrector_localized.evaluate(quadrature[quadraturePoint],
                                                          value_global_dirichlet_corrector);

              JacobianRangeType grad_global_neumann_corrector;
              glob_neumann_corrector_localized.jacobian(quadrature[quadraturePoint], grad_global_neumann_corrector);

              RangeType value_global_neumann_corrector;
              glob_neumann_corrector_localized.evaluate(quadrature[quadraturePoint], value_global_neumann_corrector);

              RangeType position = value_previous_solution + value_local_dirichlet_extension +
                                   value_global_dirichlet_corrector - value_global_neumann_corrector;
              JacobianRangeType direction(0.0);
              direction[0] = grad_previous_solution[0] + grad_local_dirichlet_extension[0] +
                             grad_global_dirichlet_corrector[0] - grad_global_neumann_corrector[0];

              RangeType value_F_x(0.0);
              nonlinear_term.evaluate(global_point, position, direction, value_F_x);

              RangeType f_x(0.0);
              f.evaluate(global_point, f_x);

              JacobianRangeType diffusive_flux_in_grad_previous_solution_direction(0.0);
              diffusion_op.diffusiveFlux(global_point, direction, diffusive_flux_in_grad_previous_solution_direction);

              newton_step_rhs[col] += weight * (value_func * f_x);
              newton_step_rhs[col] -= weight * (diffusive_flux_in_grad_previous_solution_direction[0] * grad_func[0]);
              newton_step_rhs[col] -= weight * (value_func * value_F_x);
            }
          }
        }

        newton_error = 0.0;
        for (size_t col = 0; col != newton_solution_vector.N(); ++col) {
          newton_error += (newton_step_rhs[col] * newton_step_rhs[col]);
        }
        newton_error = sqrt(newton_error);

        if (newton_error <= newton_tolerance)
          stop_newton_cycle = true;

        if (newton_error <= previous_newton_error) {
          previous_newton_error = newton_error;
          std::cout << "Newton error = " << newton_error << std::endl;
        } else {
          damping_parameter = damping_parameter / 2.0;
          std::cout << "Newton error = " << newton_error << ". Repeat with damping parameter = " << damping_parameter
                    << std::endl;
        }
      }

      for (size_t col = 0; col != newton_solution_vector.N(); ++col)
        newton_solution_vector[col] += (damping_parameter * newton_step_solution_vector[col]);

      // current solution:
      solution.clear();
      for (int internal_id = 0; internal_id < number_of_relevant_coarse_nodes; internal_id += 1) {
        DiscreteFunction aux("auxilliary function", fine_space);
        aux.clear();
        aux += *(msfem_basis_function[internal_id]);
        aux *= newton_solution_vector[internal_id];
        solution += aux;
      }

      first_cycle = false;
      iteration += 1;

    } // end while ( stop_newton_cycle == false )

    solution += global_dirichlet_corrector;
    solution -= global_neumann_corrector;

    fine_scale_part.clear();
    coarse_scale_part.clear();
    for (int internal_id = 0; internal_id < number_of_relevant_coarse_nodes; internal_id += 1) {
      DiscreteFunction aux("auxilliary function", fine_space);
      aux.clear();
      aux += *(standard_basis_function[internal_id]);
      aux *= newton_solution_vector[internal_id];
      coarse_scale_part += aux;
    }

    fine_scale_part.assign(solution);
    fine_scale_part -= coarse_scale_part;

  } // end if nonlinear

} // solve

} // namespace MsFEM {
} // namespace Multiscale {
} // namespace Dune {
