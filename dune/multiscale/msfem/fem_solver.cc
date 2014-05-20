#include <config.h>

#include <dune/common/exceptions.hh>
#include <dune/common/timer.hh>
#include <dune/fem/function/common/function.hh>
#include <dune/multiscale/common/newton_rhs.hh>
#include <dune/multiscale/common/righthandside_assembler.hh>
#include <dune/multiscale/common/traits.hh>
#include <dune/multiscale/fem/elliptic_fem_matrix_assembler.hh>
#include <dune/multiscale/fem/fem_traits.hh>
#include <dune/stuff/common/logging.hh>
#include <dune/stuff/functions/norm.hh>
#include <dune/stuff/common/parameter/configcontainer.hh>
#include <dune/pdelab/finiteelementmap/qkfem.hh>
#include <dune/pdelab/constraints/common/constraints.hh>
#include <dune/pdelab/constraints/conforming.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridoperator/gridoperator.hh>

#include <dune/pdelab/stationary/linearproblem.hh>
#include <dune/pdelab/finiteelementmap/qkfem.hh>
#include <dune/pdelab/backend/seqistlsolverbackend.hh>
#include <dune/pdelab/backend/ovlpistlsolverbackend.hh>

#include <limits>
#include <sstream>
#include <string>

#include "dune/multiscale/common/dirichletconstraints.hh"
#include "fem_solver.hh"

namespace Dune {
namespace Multiscale {

Elliptic_FEM_Solver::Elliptic_FEM_Solver(const CommonTraits::GridFunctionSpaceType &space)
  : space_(space) {}

void Elliptic_FEM_Solver::solve_linear(const CommonTraits::DiffusionType& diffusion_op,
                                       const std::unique_ptr<const CommonTraits::LowerOrderTermType>& lower_order_term,
                                       const CommonTraits::FirstSourceType& f,
                                       CommonTraits::PdelabVectorType& solution) const {
  DSC_LOG_INFO << "Solving linear problem with standard FEM\n";

  // to assemble the computational time
  Dune::Timer timer;

  typedef CommonTraits::FieldType Real;
  typedef CommonTraits::GridFunctionSpaceType GFS;
  typedef typename GFS::template ConstraintsContainer<Real>::Type CC;
  CC constraints_container;
  constraints_container.clear();
  const auto& bc_type = Problem::getModelData()->boundaryInfo();
  Dune::PDELab::constraints(*bc_type, space_, constraints_container);

  FEM::Local_CG_FEM_Operator local_operator(diffusion_op, f, lower_order_term);

  const int magic_number_stencil = 9;
  typedef BackendChooser<CommonTraits::DiscreteFunctionSpaceType>::MatrixBackendType MatrixBackendType;
  MatrixBackendType fem_matrix(magic_number_stencil);

  typedef Dune::PDELab::GridOperator<GFS,GFS,FEM::Local_CG_FEM_Operator,
                                      MatrixBackendType,Real,Real,Real,CC,CC> GridOperatorType;
  GridOperatorType global_operator(space_, constraints_container, space_, constraints_container, local_operator, fem_matrix);

  DSC_LOG_DEBUG << GridOperatorType::Traits::Jacobian(global_operator).patternStatistics() << std::endl;

  Dune::PDELab::interpolate(DS::pdelabAdapted(*Problem::getDirichletData(), space_.gridView()), space_, solution);
//  typedef Dune::PDELab::ISTLBackend_BCGS_AMG_ILU0<GridOperatorType> LinearSolverType;
//  LinearSolverType ls(space_, 5000);

  solution *= 0;
//  typedef Dune::PDELab::ISTLBackend_OVLP_CG_SSORk<CommonTraits::GridFunctionSpaceType, CC> LinearSolverType;
//  LinearSolverType ls(space_, constraints_container, 5000 /*iter*/, 5 /*steps*/, DSC_CONFIG_GET("global.cgsolver_verbose", false));
  typedef Dune::PDELab::ISTLBackend_BCGS_AMG_ILU0<GridOperatorType> LinearSolverType;
  LinearSolverType ls(space_, 5000 /*iter*/, DSC_CONFIG_GET("global.cgsolver_verbose", false));
  typedef Dune::PDELab::StationaryLinearProblemSolver<GridOperatorType,LinearSolverType,CommonTraits::PdelabVectorType> SLP;
  SLP slp(global_operator,ls,solution,1e-10);
  slp.apply();

  DSC_LOG_INFO << "---------------------------------------------------------------------------------" << std::endl;
  DSC_LOG_INFO << "Standard FEM problem solved in " << timer.elapsed() << "s." << std::endl << std::endl
               << std::endl;
}

void
Elliptic_FEM_Solver::solve_nonlinear(const CommonTraits::DiffusionType& /*diffusion_op*/,
                                     const std::unique_ptr<const CommonTraits::LowerOrderTermType>& /*lower_order_term*/,
                                     const CommonTraits::FirstSourceType& /*f*/,
                                     CommonTraits::PdelabVectorType & /*solution*/) const {
  DSC_LOG_INFO << "Solving non-linear problem.\n";
  DSC_LOG_INFO << "Solving nonlinear problem with FEM + Newton-Method.\n";
  DSC_LOG_INFO << "---------------------------------------------------------------------------------" << std::endl;

  DUNE_THROW(NotImplemented, "");

}

void Elliptic_FEM_Solver::apply(const CommonTraits::DiffusionType& diffusion_op,
                                const std::unique_ptr<const CommonTraits::LowerOrderTermType>& lower_order_term,
                                const CommonTraits::FirstSourceType& f, CommonTraits::PdelabVectorType& solution) const {

  if (Problem::getModelData()->linear())
    solve_linear(diffusion_op, lower_order_term, f, solution);
  else
    solve_nonlinear(diffusion_op, lower_order_term, f, solution);
}

} // namespace Multiscale
}
