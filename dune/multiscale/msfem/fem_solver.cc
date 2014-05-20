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

  solution *= 0;

  typedef typename GridOperatorType::Traits::TrialGridFunctionSpace GFS;
  typedef PDELab::istl::ParallelHelper<GFS> PHELPER;
  typedef typename GridOperatorType::Traits::Jacobian M;
  typedef typename M::BaseT MatrixType;
  typedef typename GridOperatorType::Traits::Domain V;
  typedef typename V::BaseT VectorType;
  typedef Dune::Amg::CoarsenCriterion<Dune::Amg::SymmetricCriterion<MatrixType,
  Dune::Amg::FirstDiagonal> > Criterion;
  typedef SeqILU0<MatrixType,VectorType,VectorType,1> Smoother;
  static constexpr int s = 96;
  typedef typename PDELab::istl::CommSelector<s,Dune::MPIHelper::isFake>::type Comm;
  typedef Dune::BlockPreconditioner<VectorType,VectorType,Comm,Smoother> ParSmoother;
  typedef Dune::OverlappingSchwarzOperator<MatrixType,VectorType,VectorType,Comm> Operator;
  typedef Dune::Amg::AMG<Operator,VectorType,ParSmoother,Comm> AMG;
  typedef typename Dune::Amg::SmootherTraits<ParSmoother>::Arguments SmootherArgs;

  Timer watch;
  Comm oocc(space_.gridView().comm());
  M mat(global_operator);
  global_operator.jacobian(solution, mat);
  PHELPER phelper(space_);
  phelper.createIndexSetAndProjectForAMG(mat, oocc);
  Operator oop(mat.base(), oocc);
  Dune::OverlappingSchwarzScalarProduct<VectorType,Comm> sp(oocc);

  SmootherArgs smootherArgs;
  smootherArgs.iterations = 1;
  smootherArgs.relaxationFactor = 1;
  Dune::Amg::Parameters params(15,2000);
  params.setDefaultValuesIsotropic(GFS::Traits::GridViewType::Traits::Grid::dimension);
  params.setDebugLevel(4);
  Criterion criterion(params);

  AMG amg(oop, criterion, smootherArgs, oocc);

  watch.reset();
  Dune::BiCGSTABSolver<VectorType> solver(oop, sp, amg, 1e-8/*reduction*/, 5000/*maxiter*/, 1 /*verb*/);
  Dune::InverseOperatorResult stat;

  auto rhs = solution;
  global_operator.residual(solution, rhs);
  rhs *= -1;
  solver.apply(solution,rhs, stat);


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
