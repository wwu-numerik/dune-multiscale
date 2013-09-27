// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DiscreteEllipticMsFEMConservativeFluxSolver_HH
#define DiscreteEllipticMsFEMConservativeFluxSolver_HH

#include <config.h>
#include <vector>

#include <dune/multiscale/msfem/msfem_traits.hh>

#include <dune/multiscale/tools/subgrid_io.hh>
#include <dune/multiscale/tools/discretefunctionwriter.hh>
#include <dune/multiscale/tools/misc.hh>
#include <dune/multiscale/tools/misc/outputparameter.hh>
#include <dune/multiscale/hmm/cell_problem_numbering.hh>

#include <dune/subgrid/subgrid.hh>

#include <dune/fem/operator/matrix/spmatrix.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/gridpart/common/gridpart.hh>
#include <dune/fem/operator/linear/petscoperator.hh>
#include <dune/fem/solver/cginverseoperator.hh>
#include <dune/fem/operator/common/stencil.hh>

#include <dune/stuff/common/profiler.hh>
#include <dune/stuff/fem/localmatrix_proxy.hh>
#include <dune/stuff/fem/matrix_object.hh>
#include <dune/stuff/common/math.hh>

namespace Dune {
namespace Multiscale {
namespace MsFEM {

static const bool FLUX_SOLVER_VERBOSE = false;

/** \brief define output parameters for local problems
 *  appends "local_problems" for path
 **/
struct ConFluxProblemDataOutputParameters : public OutputParameters {
public:
  explicit ConFluxProblemDataOutputParameters()
    : OutputParameters(DSC_CONFIG_GET("global.datadir", "data") + "/cf_problems/") {}
};

//! \TODO docme
template <class SubGridDiscreteFunctionType, class DiscreteFunctionType, class DiffusionOperatorType,
          class MacroMicroGridSpecifierType>
class ConservativeFluxOperator : public Operator<typename SubGridDiscreteFunctionType::RangeFieldType,
                                                 typename SubGridDiscreteFunctionType::RangeFieldType,
                                                 SubGridDiscreteFunctionType, SubGridDiscreteFunctionType> {
  typedef ConservativeFluxOperator<SubGridDiscreteFunctionType, DiscreteFunctionType, DiffusionOperatorType,
                                   MacroMicroGridSpecifierType> This;

private:
  typedef SubGridDiscreteFunctionType SubGridDiscreteFunction;
  typedef DiscreteFunctionType DiscreteFunction;
  typedef DiffusionOperatorType DiffusionModel;

  typedef typename SubGridDiscreteFunction::DiscreteFunctionSpaceType SubGridDiscreteFunctionSpace;
  typedef typename DiscreteFunction::DiscreteFunctionSpaceType DiscreteFunctionSpace;

  typedef typename DiscreteFunctionSpace::GridPartType GridPart;
  typedef typename SubGridDiscreteFunctionSpace::GridPartType SubGridPart;

  typedef typename DiscreteFunctionSpace::RangeFieldType RangeFieldType;
  typedef typename DiscreteFunctionSpace::DomainType DomainType;
  typedef typename DiscreteFunctionSpace::RangeType RangeType;
  typedef typename DiscreteFunctionSpace::JacobianRangeType JacobianRangeType;

  static const int dimension = GridPart::GridType::dimension;
  static const int polynomialOrder = DiscreteFunctionSpace::polynomialOrder;

  typedef typename DiscreteFunctionSpace::BasisFunctionSetType BaseFunctionSet;
  typedef typename GridPart::template Codim<0>::EntityType Entity;
  typedef typename Entity::EntityPointer EntityPointer;
  typedef typename SubGridDiscreteFunctionSpace::BasisFunctionSetType SubGridBaseFunctionSet;
  typedef typename SubGridDiscreteFunctionSpace::EntityType SubGridEntity;

public:
  ConservativeFluxOperator(const SubGridDiscreteFunctionSpace& subDiscreteFunctionSpace,
                           const DiscreteFunctionSpace& discreteFunctionSpace, const DiffusionModel& diffusion_op,
                           const MacroMicroGridSpecifierType& specifier)
    : subDiscreteFunctionSpace_(subDiscreteFunctionSpace)
    , discreteFunctionSpace_(discreteFunctionSpace)
    , diffusion_operator_(diffusion_op)
    , specifier_(specifier) {}

private:
  ConservativeFluxOperator(const This&) = delete;

public:
  // dummy operator
  virtual void operator()(const SubGridDiscreteFunction& u, SubGridDiscreteFunction& w) const;

  // assemble stiffness matrix for local problems
  template <class MatrixType>
  void assemble_matrix(const int sub_grid_id, MatrixType& global_matrix) const;

  // the right hand side assembler methods
  void assemble_RHS( // direction 'e_i'
      JacobianRangeType& e_i,
      // solution of the local corrector problem
      const SubGridDiscreteFunction& local_corrector_e_i, const int sub_grid_id,
      // rhs local msfem problem:
      SubGridDiscreteFunction& rhs_flux_problem) const;

  double normRHS(SubGridDiscreteFunction& rhs) const;

private:
  const SubGridDiscreteFunctionSpace& subDiscreteFunctionSpace_;
  const DiscreteFunctionSpace& discreteFunctionSpace_;
  const DiffusionModel& diffusion_operator_;
  const MacroMicroGridSpecifierType& specifier_;
};

//! dummy implementation of "operator()"
//! 'w' = effect of the discrete operator on 'u'
template <class SubGridDiscreteFunctionImp, class DiscreteFunctionImp, class DiffusionImp,
          class MacroMicroGridSpecifierImp>
void
ConservativeFluxOperator<SubGridDiscreteFunctionImp, DiscreteFunctionImp, DiffusionImp, MacroMicroGridSpecifierImp>::
operator()(const SubGridDiscreteFunctionImp& /*u*/, SubGridDiscreteFunctionImp& /*w*/) const {
  DUNE_THROW(Dune::NotImplemented,
             "the ()-operator of the ConservativeFluxOperator class is not yet implemented and still a dummy.");
}

//! assemble system matrix
template <class SubGridDiscreteFunctionImp, class DiscreteFunctionImp, class DiffusionImp,
          class MacroMicroGridSpecifierImp>
template <class MatrixType>
void ConservativeFluxOperator<SubGridDiscreteFunctionImp, DiscreteFunctionImp, DiffusionImp,
                              MacroMicroGridSpecifierImp>::assemble_matrix(const int sub_grid_id,
                                                                           MatrixType& global_matrix) const {
  global_matrix.reserve(DSFe::diagonalAndNeighborStencil(global_matrix));
  global_matrix.clear();

  // local grid basis functions:
  std::vector<RangeType> phi(subDiscreteFunctionSpace_.mapper().maxNumDofs());

  for (const auto& sub_grid_entity : subDiscreteFunctionSpace_) {
    auto host_entity_pointer =
        subDiscreteFunctionSpace_.gridPart().grid().template getHostEntity<0>(sub_grid_entity);

    const auto& coarseGridLeafIndexSet = specifier_.coarseSpace().gridPart().grid().leafIndexSet();

    auto father_of_sub_grid_entity =
        DSG::make_father(coarseGridLeafIndexSet, host_entity_pointer, specifier_.getLevelDifference());
    const int coarse_index = coarseGridLeafIndexSet.index(*father_of_sub_grid_entity);
    assert(sub_grid_entity.partitionType() == InteriorEntity);

    DSFe::LocalMatrixProxy<MatrixType> local_matrix(global_matrix, sub_grid_entity, sub_grid_entity);

    const auto& baseSet = local_matrix.domainBasisFunctionSet();
    const auto numBaseFunctions = baseSet.size();

    for (const auto& intersection : DSC::intersectionRange(discreteFunctionSpace_.gridPart(), *host_entity_pointer)) {
      const auto faceQuadrature = make_quadrature(intersection, discreteFunctionSpace_);
      const auto& faceGeometry = intersection.geometry();

      bool set_zero = false;
      if (coarse_index != sub_grid_id) {
        set_zero = true;
      }

      if ((intersection.neighbor()) && (!set_zero)) {
        const auto outside_it = intersection.outside();

        const auto father_of_neighbor =
            DSG::make_father(coarseGridLeafIndexSet, outside_it, specifier_.getLevelDifference());
        if (DSG::entities_identical(*father_of_sub_grid_entity, *father_of_neighbor)) {
          set_zero = true;
        }
      }

      const auto numQuadraturePoints = faceQuadrature.nop();

      RangeType check_sum(0.0);

      for (size_t quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint) {
        // das liefert immer den gleichen 'local_point' egal fuer welchen Quadraturpunkt. // WTF ?!?!
        const auto local_point = faceGeometry.local(faceQuadrature.point(quadraturePoint));
        const double integrationFactor = faceGeometry.integrationElement(local_point);
        const double quadratureWeight = faceQuadrature.weight(quadraturePoint);

        check_sum += integrationFactor * quadratureWeight;

        baseSet.evaluateAll(faceQuadrature[quadraturePoint], phi);

        for (unsigned int i = 0; i < numBaseFunctions; ++i) {
          for (unsigned int j = 0; j < numBaseFunctions; ++j) {
            if (!set_zero) {
              local_matrix.add(j, i, integrationFactor * quadratureWeight * (phi[i][0] * phi[j][0]));
            } else {
              // stabilization (should be close to zero):
              local_matrix.add(j, i,
                               /*0.0*/ 0.00000001 * integrationFactor * quadratureWeight * (phi[i][0] * phi[j][0]));
            }
          }
        }
      } // done loop over all quadrature points

      if (check_sum != faceGeometry.volume()) {
        DUNE_THROW(Dune::InvalidStateException, "Error in Face Quadrature.");
      }
    }
  }
} // assemble_matrix

template <class SubGridDiscreteFunctionImp, class DiscreteFunctionImp, class DiffusionImp,
          class MacroMicroGridSpecifierImp>
double ConservativeFluxOperator<SubGridDiscreteFunctionImp, DiscreteFunctionImp, DiffusionImp,
                                MacroMicroGridSpecifierImp>::normRHS(SubGridDiscreteFunctionImp& rhs) const {
  double norm = 0.0;
  const auto& discreteFunctionSpace = rhs.space();

  for (const auto& entity : discreteFunctionSpace) {
    const auto quadrature = make_quadrature(entity, discreteFunctionSpace);
    const auto& geo = entity.geometry();
    const auto localRHS = rhs.localFunction(entity);
    // integrate
    for (auto quadraturePoint : DSC::valueRange(quadrature.nop())) {
      const double weight =
          quadrature.weight(quadraturePoint) * geo.integrationElement(quadrature.point(quadraturePoint));

      RangeType value(0.0);
      localRHS.evaluate(quadrature[quadraturePoint], value);

      norm += weight * value * value;
    }
  }

  return norm;
} // end method

// assemble the right hand side of the conservative flux problem
// ----------------------------------------------

template <class SubGridDiscreteFunctionImp, class DiscreteFunctionImp, class DiffusionImp,
          class MacroMicroGridSpecifierImp>
// template< class MatrixType >
    void ConservativeFluxOperator<SubGridDiscreteFunctionImp, DiscreteFunctionImp, DiffusionImp,
                                  MacroMicroGridSpecifierImp>::assemble_RHS( // direction 'e'
        JacobianRangeType& e_i,
        // solution of the local corrector problem
        const SubGridDiscreteFunction& local_corrector_e_i, const int sub_grid_id,
        // rhs flux problem:
        SubGridDiscreteFunction& rhs_flux_problem) const {
  const SubGridDiscreteFunctionSpace& subDiscreteFunctionSpace = rhs_flux_problem.space();

  // set entries to zero:
  rhs_flux_problem.clear();

  // gradient of micro scale base function:
  std::vector<JacobianRangeType> gradient_phi(subDiscreteFunctionSpace.mapper().maxNumDofs());

  for (const auto& local_grid_entity : subDiscreteFunctionSpace) {
    auto host_entity_pointer =
        subDiscreteFunctionSpace.gridPart().grid().template getHostEntity<0>(local_grid_entity);

    const auto& coarseGridLeafIndexSet = specifier_.coarseSpace().gridPart().grid().leafIndexSet();

    auto father_of_sub_grid_entity =
        DSG::make_father(coarseGridLeafIndexSet, host_entity_pointer, specifier_.getLevelDifference());
    const int coarse_index = coarseGridLeafIndexSet.index(*father_of_sub_grid_entity);

    if (coarse_index != sub_grid_id) {
      continue;
    }

    const auto& geometry = local_grid_entity.geometry();
    assert(local_grid_entity.partitionType() == InteriorEntity);

    auto elementOfRHS = rhs_flux_problem.localFunction(local_grid_entity);

    const auto& baseSet = elementOfRHS.basisFunctionSet();
    const auto numBaseFunctions = baseSet.size();

    const auto quadrature = make_quadrature(local_grid_entity, subDiscreteFunctionSpace);
    const auto numQuadraturePoints = quadrature.nop();
    for (size_t quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint) {
      const auto& local_point = quadrature.point(quadraturePoint);

      // remember, we are concerned with: - \int_{U(T)} (A^eps)(x) e · ∇ \phi(x)

      // global point in the subgrid
      const auto global_point = geometry.global(local_point);
      const double weight = quadrature.weight(quadraturePoint) * geometry.integrationElement(local_point);

      // A^eps(x) ( e
      // diffusion operator evaluated in 'x' multiplied with e
      JacobianRangeType diffusion_in_e_i;
      diffusion_operator_.diffusiveFlux(global_point, e_i, diffusion_in_e_i);

      JacobianRangeType total_diffusive_flux;

      JacobianRangeType grad_corrector_e_i;
      const auto localized_corrector_e_i = local_corrector_e_i.localFunction(local_grid_entity);
      localized_corrector_e_i.jacobian(quadrature[quadraturePoint], grad_corrector_e_i);
      diffusion_operator_.diffusiveFlux(global_point, grad_corrector_e_i, total_diffusive_flux);

      total_diffusive_flux[0] += diffusion_in_e_i[0];

      baseSet.jacobianAll(quadrature[quadraturePoint], gradient_phi);

      for (unsigned int i = 0; i < numBaseFunctions; ++i) {
        elementOfRHS[i] -= weight * (total_diffusive_flux[0] * gradient_phi[i][0]);
      }
    }
  } // end for-loop of 'Iterator'
} // assemble_RHS

//! the essential local msfem problem solver class
template <class SubGridDiscreteFunctionType, class HostDiscreteFunctionType, class DiffusionOperatorType,
          class MacroMicroGridSpecifierImp>
class ConservativeFluxProblemSolver {
private:
  typedef MacroMicroGridSpecifierImp MacroMicroGridSpecifierType;

  //! ---------------- typedefs for the HostDiscreteFunctionSpace -----------------------
  typedef typename HostDiscreteFunctionType::DiscreteFunctionSpaceType HostDiscreteFunctionSpaceType;
  typedef typename HostDiscreteFunctionSpaceType::FunctionSpaceType FunctionSpaceType;
  typedef typename HostDiscreteFunctionSpaceType::GridPartType HostGridPartType;
  typedef typename HostDiscreteFunctionSpaceType::GridType HostGridType;
  typedef typename HostDiscreteFunctionSpaceType::EntityType HostEntityType;
  typedef typename HostEntityType::EntityPointer HostEntityPointerType;

  static const int dimension = HostGridType::dimension;

  //! ---------------- typedefs for the SubGridDiscreteFunctionSpace -----------------------
  typedef typename SubGridDiscreteFunctionType::DiscreteFunctionSpaceType SubGridDiscreteFunctionSpaceType;
  typedef typename SubGridDiscreteFunctionSpaceType::RangeFieldType RangeFieldType;
  typedef typename SubGridDiscreteFunctionSpaceType::DomainType DomainType;
  typedef typename SubGridDiscreteFunctionSpaceType::RangeType RangeType;
  typedef typename SubGridDiscreteFunctionSpaceType::JacobianRangeType JacobianRangeType;
  typedef typename SubGridDiscreteFunctionSpaceType::EntityType SubGridEntityType;
  typedef typename SubGridEntityType::EntityPointer SubGridEntityPointerType;
  //!-----------------------------------------------------------------------------------------

  //! ------------------ Matrix Traits for the local Problems ---------------------
  static const int faceCodim = 1;
  //! polynomial order of base functions
  static const int polynomialOrder = SubGridDiscreteFunctionSpaceType::polynomialOrder;

  typedef Dune::Fem::PetscLinearOperator<SubGridDiscreteFunctionType, SubGridDiscreteFunctionType> FluxProbFEMMatrix;
  typedef Dune::Fem::PetscInverseOperator<SubGridDiscreteFunctionType, FluxProbFEMMatrix> InverseFluxProbFEMMatrix;

private:
  const DiffusionOperatorType& diffusion_;
  const HostDiscreteFunctionSpaceType& hostDiscreteFunctionSpace_;

  const MacroMicroGridSpecifierType& specifier_;
  mutable boost::format filename_template_;

public:
  //! constructor - with diffusion operator A^{\epsilon}(x)
  ConservativeFluxProblemSolver(const HostDiscreteFunctionSpaceType& hostDiscreteFunctionSpace,
                                const DiffusionOperatorType& diffusion_operator,
                                const MacroMicroGridSpecifierType& specifier)
    : diffusion_(diffusion_operator)
    , hostDiscreteFunctionSpace_(hostDiscreteFunctionSpace)
    , specifier_(specifier)
    , filename_template_(boost::format("_conservativeFlux_e_%d_sg_%d")) {}

  //! ----------- method: solve the local MsFEM problem ------------------------------------------

  void solve(JacobianRangeType& e_i, // direction 'e_i'
             const SubGridDiscreteFunctionType& local_corrector_e_i, const int sub_grid_id, const int direction_index,
             SubGridDiscreteFunctionType& conservative_flux) const {
    // set solution equal to zero:
    conservative_flux.clear();

    const SubGridDiscreteFunctionSpaceType& localDiscreteFunctionSpace = local_corrector_e_i.space();

    //! the matrix in our linear system of equations
    FluxProbFEMMatrix flux_prob_system_matrix("Conservative Flux Problem System Matrix", localDiscreteFunctionSpace,
                                              localDiscreteFunctionSpace);

    //! define the discrete (elliptic) local MsFEM problem operator
    // ( effect of the discretized differential operator on a certain discrete function )
    // discrete elliptic operator describing the elliptic local msfem problems
    typedef ConservativeFluxOperator<SubGridDiscreteFunctionType, HostDiscreteFunctionType, DiffusionOperatorType,
                                     MacroMicroGridSpecifierType> ConservativeFluxOperatorType;
    ConservativeFluxOperatorType cf_problem_operator(localDiscreteFunctionSpace, hostDiscreteFunctionSpace_, diffusion_,
                                                     specifier_);

    //! right hand side vector of the algebraic local MsFEM problem
    SubGridDiscreteFunctionType rhs("RHS of Conservative Flux Problem", localDiscreteFunctionSpace);
    rhs.clear();

    // assemble the stiffness matrix
    cf_problem_operator.assemble_matrix(sub_grid_id, flux_prob_system_matrix);

    // assemble right hand side of algebraic local msfem problem
    cf_problem_operator.assemble_RHS(e_i, local_corrector_e_i, sub_grid_id, rhs);

    const double norm_rhs = cf_problem_operator.normRHS(rhs);

    if (!(rhs.dofsValid())) {
      DUNE_THROW(Dune::InvalidStateException, "Local Flux Problem RHS invalid.");
    }

    if (norm_rhs < /*1e-06*/ 1e-30) {
      conservative_flux.clear();
      DSC_LOG_INFO << "Local Flux Problem with solution zero." << std::endl;
    } else {
      InverseFluxProbFEMMatrix flux_prob_biCGStab(flux_prob_system_matrix, 1e-8, 1e-8, 20000, FLUX_SOLVER_VERBOSE, "cg",
                                                  DSC_CONFIG_GET("preconditioner_type", std::string("sor")));
      flux_prob_biCGStab(rhs, conservative_flux);
    }

    if (!(conservative_flux.dofsValid())) {
      DUNE_THROW(Dune::InvalidStateException, "Solution of the Local Flux Problem is invalid!");
    }

#ifdef VTK_OUTPUT
    vtk_output(conservative_flux, sub_grid_id, direction_index);
#endif
    file_data_output(conservative_flux, sub_grid_id, direction_index);
  } // solve

  //! ----------- end method: solve local MsFEM problem ------------------------------------------

  // create a hostgrid function from a subgridfunction
  void subgrid_to_hostrid_function(const SubGridDiscreteFunctionType& sub_func,
                                   HostDiscreteFunctionType& host_func) const {
    host_func.clear();
    const auto& subDiscreteFunctionSpace = sub_func.space();
    const auto& subGrid = subDiscreteFunctionSpace.grid();

    for (const auto& sub_entity : subDiscreteFunctionSpace) {
      const auto host_entity_pointer = subGrid.template getHostEntity<0>(sub_entity);
      const auto& host_entity = *host_entity_pointer;

      const auto sub_loc_value = sub_func.localFunction(sub_entity);
      auto host_loc_value = host_func.localFunction(host_entity);

      const auto numBaseFunctions = sub_loc_value.basisFunctionSet().size();
      for (unsigned int i = 0; i < numBaseFunctions; ++i) {
        host_loc_value[i] = sub_loc_value[i];
      }
    }
  } // subgrid_to_hostrid_function

  void vtk_output(const SubGridDiscreteFunctionType& subgrid_disc_func, const int sub_grid_index,
                  const int direction_index) const {
    // --------------- writing data output ---------------------
    // typedefs and initialization
    typedef std::tuple<HostDiscreteFunctionType*> IOTupleType;
    typedef Fem::DataOutput<HostGridType, IOTupleType> DataOutputType;

    // general output parameters
    ConFluxProblemDataOutputParameters outputparam;
    // -------------------------------------------------------

    HostDiscreteFunctionType host_disc_func("Conservative Flux", hostDiscreteFunctionSpace_);
    subgrid_to_hostrid_function(subgrid_disc_func, host_disc_func);

    // create and initialize output class
    IOTupleType conservative_flux_series(&host_disc_func);
    filename_template_ % direction_index % sub_grid_index;
    outputparam.set_prefix(filename_template_.str());
    DataOutputType cf_dataoutput(hostDiscreteFunctionSpace_.gridPart().grid(), conservative_flux_series, outputparam);

    cf_dataoutput.writeData(1.0 /*dummy*/, "conservative-flux");
  } // vtk_output

  void file_data_output(const SubGridDiscreteFunctionType& subgrid_disc_func, const int sub_grid_index,
                        const int direction_index) const {
    const std::string locprob_solution_location =
        std::string("cf_problems/") + (filename_template_ % direction_index % sub_grid_index).str();
    DiscreteFunctionWriter(locprob_solution_location).append(subgrid_disc_func);
  } // file_data_output

  void solve_all(MsFEMTraits::SubGridListType& subgrid_list) const {
    JacobianRangeType e[dimension];

    for (int i = 0; i < dimension; ++i) {
      for (int j = 0; j < dimension; ++j) {
        if (i == j) {
          e[i][0][j] = 1.0;
        } else {
          e[i][0][j] = 0.0;
        }
      }
    }

    // number of coarse grid entities (of codim 0).
    const auto number_of_coarse_grid_entities = specifier_.getNumOfCoarseEntities();

    DSC_PROFILER.startTiming("msfem.conservative_flux_solver.solve_all_subgrids");

    // we want to determine minimum, average and maxiumum time for solving a local msfem problem in the current method
    Dune::Stuff::Common::MinMaxAvg<double> cell_time;

    const auto& coarseSpace = specifier_.coarseSpace();
    const auto& coarseGridLeafIndexSet = coarseSpace.gridPart().grid().leafIndexSet();

    for (const auto& host_entity : coarseSpace) {
      const auto global_index_entity = coarseGridLeafIndexSet.index(host_entity);

      // the sub grid U(T) that belongs to the coarse_grid_entity T
      auto subGridPart = subgrid_list.gridPart(global_index_entity);

      SubGridDiscreteFunctionSpaceType localDiscreteFunctionSpace(subGridPart);

      SubGridDiscreteFunctionType local_problem_solution_e0("Local problem Solution e_0", localDiscreteFunctionSpace);
      local_problem_solution_e0.clear();

      SubGridDiscreteFunctionType local_problem_solution_e1("Local problem Solution e_1", localDiscreteFunctionSpace);
      local_problem_solution_e1.clear();

      // --------- load local solutions -------
      // the file/place, where we saved the solutions of the cell problems
      const std::string local_solution_location =
          (boost::format("local_problems/_localProblemSolutions_%d") %
           coarseSpace.gridPart().grid().globalIdSet().id(host_entity)).str();

      // reader for the cell problem data file:
      DiscreteFunctionReader discrete_function_reader(local_solution_location);
      discrete_function_reader.read(0, local_problem_solution_e0);
      discrete_function_reader.read(1, local_problem_solution_e1);

      SubGridDiscreteFunctionType conservative_flux_e0("Conservative Flux for e_0", localDiscreteFunctionSpace);
      SubGridDiscreteFunctionType conservative_flux_e1("Conservative Flux for e_1", localDiscreteFunctionSpace);

      DSC_LOG_INFO << "Number of the 'conservative flux problem': " << dimension* global_index_entity << " (of "
                   << (dimension * number_of_coarse_grid_entities) - 1 << " problems in total)" << std::endl;

      // take time
      DSC_PROFILER.startTiming("none.local_problem_solution");

      this->solve(e[0], local_problem_solution_e0, global_index_entity, 0, conservative_flux_e0);

      DSC_LOG_INFO << "Number of the 'conservative flux problem': " << (dimension * global_index_entity) + 1 << " (of "
                   << (dimension * number_of_coarse_grid_entities) - 1 << " problems in total)" << std::endl;

      cell_time(DSC_PROFILER.stopTiming("none.local_problem_solution") / 1000.f);
      DSC_PROFILER.resetTiming("none.local_problem_solution");
      DSC_PROFILER.startTiming("none.local_problem_solution");

      this->solve(e[1], local_problem_solution_e1, global_index_entity, 1, conservative_flux_e1);

      cell_time(DSC_PROFILER.stopTiming("none.local_problem_solution") / 1000.f);
      DSC_PROFILER.resetTiming("none.local_problem_solution");
    }
    const auto total_time = DSC_PROFILER.stopTiming("msfem.conservative_flux_solver.solve_all_subgrids") / 1000.f;
    DSC_LOG_INFO << std::endl;
    DSC_LOG_INFO << "In: 'assemble all conservatice fluxes'." << std::endl << std::endl;
    DSC_LOG_INFO << "Conservative Flux determined for " << number_of_coarse_grid_entities << " coarse grid entities."
                 << std::endl;
    DSC_LOG_INFO << dimension* number_of_coarse_grid_entities << " conservative flux problems solved in total."
                 << std::endl;
    DSC_LOG_INFO << "Minimum time for solving a conservative flux problem = " << cell_time.min() << "s." << std::endl;
    DSC_LOG_INFO << "Maximum time for solving a conservative flux problem = " << cell_time.max() << "s." << std::endl;
    DSC_LOG_INFO << "Average time for solving a conservative flux problem = " << cell_time.average() << "s."
                 << std::endl;
    DSC_LOG_INFO << "Total time for computing and saving the conservative flux problems = " << total_time << "s,"
                 << std::endl << std::endl;

  } // solve_all
};  // end class

} // namespace MsFEM {
} // namespace Multiscale {
} // namespace Dune {

#endif // #ifndef DiscreteEllipticMsFEMLocalProblem_HH
