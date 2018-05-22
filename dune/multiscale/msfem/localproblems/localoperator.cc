#include <config.h>
// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "localoperator.hh"

#include <assert.h>
#include <boost/assert.hpp>
#include <dune/common/fmatrix.hh>
#include <dune/common/exceptions.hh>
#include <dune/xt/common/filesystem.hh>
#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/configuration.hh>
#include <dune/xt/common/timings.hh>
#include <dune/multiscale/common/heterogenous.hh>
#include <dune/multiscale/tools/misc.hh>
#include <dune/multiscale/problems/selector.hh>
#include <dune/multiscale/msfem/localproblems/localproblemsolver.hh>
#include <dune/multiscale/msfem/msfem_traits.hh>
#include <dune/multiscale/common/traits.hh>
#include <dune/multiscale/msfem/localproblems/localproblemsolver.hh>
#include <dune/gdt/operators/projections.hh>
#include <dune/gdt/operators/prolongations.hh>
#include <dune/gdt/spaces/constraints.hh>
#include <dune/gdt/functionals/l2.hh>

namespace Dune {
namespace Multiscale {

//! holds functionals,etc for boundary correctors to make conditional usage simpler
template <class LocalNeumann>
class BoundaryValueHelper
{
  typedef GDT::Functionals::L2Face<LocalNeumann, CommonTraits::GdtVectorType, MsFEMTraits::LocalSpaceType>
      NeumannFunctional;
  typedef GDT::Functionals::L2Volume<Problem::LocalDiffusionType,
                                     CommonTraits::GdtVectorType,
                                     MsFEMTraits::LocalSpaceType,
                                     MsFEMTraits::LocalGridViewType,
                                     DirichletProduct>
      DirichletCorrectorFunctionalType;

public:
  BoundaryValueHelper(const DMP::ProblemContainer& problem,
                      const MsFEMTraits::LocalSpaceType& localSpace,
                      const Problem::LocalDiffusionType& local_diffusion_operator,
                      MsFEMTraits::LocalSolutionVectorType& allLocalRHS,
                      std::size_t coarseBaseFunc)
    : localSpace_(localSpace)
    , dirichletExtensionLocal_(localSpace_, "dirichletExtension")
    , local_neumann_(problem.getNeumannData().transfer<MsFEMTraits::LocalEntityType>())
    , neumann_functional_(local_neumann_, allLocalRHS[coarseBaseFunc]->vector(), localSpace_)
    , dl_corrector_functional_(problem.getDiffusion(), dirichletExtensionLocal_, local_diffusion_operator)
    , dirichlet_corrector_(
          local_diffusion_operator, allLocalRHS[++coarseBaseFunc]->vector(), localSpace_, dl_corrector_functional_)
    , problem_(problem)
  {
  }

  void dirichlet_projection(const CommonTraits::SpaceType& coarse_space)
  {
    GDT::SystemAssembler<CommonTraits::SpaceType> coarse_system_assembler(coarse_space);
    const auto& dirichlet_data = problem_.getDirichletData();
    CommonTraits::DiscreteFunctionType dirichletExtensionCoarse(coarse_space, "Dirichlet Extension Coarse");
    GDT::Operators::DirichletProjectionLocalizable<CommonTraits::GridViewType,
                                                   Problem::DirichletDataBase,
                                                   CommonTraits::DiscreteFunctionType>
        coarse_dirichlet_projection_operator(
            coarse_space.grid_view(), problem_.getModelData().boundaryInfo(), dirichlet_data, dirichletExtensionCoarse);
    coarse_system_assembler.add(coarse_dirichlet_projection_operator,
                                new DSG::ApplyOn::BoundaryEntities<CommonTraits::GridViewType>());
    coarse_system_assembler.assemble();
    GDT::Operators::LagrangeProlongation<MsFEMTraits::LocalGridViewType> projection(localSpace_.grid_view());
    projection.apply(dirichletExtensionCoarse, dirichletExtensionLocal_);
  }

  void add_to(GDT::SystemAssembler<MsFEMTraits::LocalSpaceType>& system_assembler)
  {
    typedef DSG::ApplyOn::BoundaryEntities<MsFEMTraits::LocalGridViewType> OnLocalBoundaryEntities;
    system_assembler.add(neumann_functional_, new OnLocalBoundaryEntities());
    system_assembler.add(dirichlet_corrector_, new OnLocalBoundaryEntities());
  }

private:
  const MsFEMTraits::LocalSpaceType& localSpace_;
  MsFEMTraits::LocalGridDiscreteFunctionType dirichletExtensionLocal_;
  LocalNeumann local_neumann_;
  NeumannFunctional neumann_functional_;
  GDT::LocalFunctional::Codim0Integral<DirichletProduct> dl_corrector_functional_;
  DirichletCorrectorFunctionalType dirichlet_corrector_;
  const DMP::ProblemContainer& problem_;
};

LocalProblemOperator::LocalProblemOperator(const DMP::ProblemContainer& problem,
                                           const CommonTraits::SpaceType& coarse_space,
                                           const MsFEMTraits::LocalSpaceType& space)
  : localSpace_(space)
  , local_diffusion_operator_(problem.getDiffusion())
  , coarse_space_(coarse_space)
  , system_matrix_(localSpace_.mapper().size(), localSpace_.mapper().size(), EllipticOperatorType::pattern(localSpace_))
  , system_assembler_(localSpace_)
  , elliptic_operator_(local_diffusion_operator_, system_matrix_, localSpace_)
  , dirichletConstraints_(problem.getModelData().subBoundaryInfo(), localSpace_.mapper().size(), true)
#if HAVE_UMFPACK
  , use_umfpack_(problem.config().get("msfem.local_solver", std::string("umfpack")) == std::string("umfpack"))
#else
  , use_umfpack_(false)
#endif
  , problem_(problem)
{
  system_assembler_.add(elliptic_operator_);
}

void LocalProblemOperator::assemble_all_local_rhs(const MsFEMTraits::CoarseEntityType& coarseEntity,
                                                  MsFEMTraits::LocalSolutionVectorType& allLocalRHS)
{
  BOOST_ASSERT_MSG(allLocalRHS.size() > 0, "You need to preallocate the necessary space outside this function!");

  const bool is_simplex_grid = DSG::is_simplex_grid(coarse_space_);
  if (is_simplex_grid)
    DUNE_THROW(NotImplemented, "special treatment for simplicial grids missing");

  const auto numBoundaryCorrectors = is_simplex_grid ? 1u : 2u;
  const auto numInnerCorrectors = allLocalRHS.size() - numBoundaryCorrectors;

  const auto hasOversamplig = problem_.config().get("msfem.oversampling_layers", 0) > 0;
  typedef BoundaryValueHelper<decltype(problem_.getNeumannData().transfer<MsFEMTraits::LocalEntityType>())> BVHelper;
  std::unique_ptr<BVHelper> bv_helper(nullptr);
  if (hasOversamplig && coarseEntity.hasBoundaryIntersections()) {
    bv_helper = Dune::XT::Common::make_unique<BVHelper>(
        problem_, localSpace_, local_diffusion_operator_, allLocalRHS, numInnerCorrectors);
    bv_helper->dirichlet_projection(coarse_space_);
  }

  typedef GDT::Functionals::L2Volume<Problem::LocalDiffusionType,
                                     CommonTraits::GdtVectorType,
                                     MsFEMTraits::LocalSpaceType,
                                     MsFEMTraits::LocalGridViewType,
                                     CoarseBasisProduct>
      RhsFunctionalType;
  std::vector<std::unique_ptr<RhsFunctionalType>> rhs_functionals(numInnerCorrectors);
  std::size_t coarseBaseFunc = 0;
  const auto coarseBaseFunctionSet = coarse_space_.base_function_set(coarseEntity);
  for (; coarseBaseFunc < numInnerCorrectors; ++coarseBaseFunc) {
    assert(allLocalRHS[coarseBaseFunc]);
    GDT::LocalFunctional::Codim0Integral<CoarseBasisProduct> local_rhs_functional(
        problem_.getDiffusion(), coarseBaseFunctionSet, local_diffusion_operator_, coarseBaseFunc);
    auto& rhs_vector = allLocalRHS[coarseBaseFunc]->vector();
    rhs_functionals[coarseBaseFunc] = Dune::XT::Common::make_unique<RhsFunctionalType>(
        local_diffusion_operator_, rhs_vector, localSpace_, local_rhs_functional);
    system_assembler_.add(*rhs_functionals[coarseBaseFunc]);
  }
  
  if (bv_helper != nullptr)
    bv_helper->add_to(system_assembler_);

  // system_assembler_.assemble();
  // dirichlet-0 for all rhs
  typedef DSG::ApplyOn::BoundaryEntities<MsFEMTraits::LocalGridViewType> OnLocalBoundaryEntities;
  system_assembler_.add(dirichletConstraints_, new OnLocalBoundaryEntities());

  DXTC_TIMINGS.start("msfem.local.rhs.assemble");
  system_assembler_.assemble();
  DXTC_TIMINGS.stop("msfem.local.rhs.assemble");

  dirichletConstraints_.apply(system_matrix_);
  for (auto& rhs : allLocalRHS)
    dirichletConstraints_.apply(rhs->vector());
#if HAVE_UMFPACK
  if (use_umfpack_)
    local_direct_inverse_ = Dune::XT::Common::make_unique<LocalDirectInverseType>(
        system_matrix_.backend(), problem_.config().get("msfem.local_solver_verbose", 0));
#endif
}

void LocalProblemOperator::apply_inverse(const MsFEMTraits::LocalGridDiscreteFunctionType& current_rhs,
                                         MsFEMTraits::LocalGridDiscreteFunctionType& current_solution)
{
  if (!current_rhs.dofs_valid())
    DUNE_THROW(Dune::InvalidStateException, "Local MsFEM Problem RHS invalid.");

  typedef BackendChooser<MsFEMTraits::LocalSpaceType>::InverseOperatorType LocalInverseOperatorType;
  const LocalInverseOperatorType local_inverse(system_matrix_, current_rhs.space().communicator());

  auto options = local_inverse.options(problem_.config().get("msfem.local_solver", "umfpack"));
  options["precision"] = problem_.config().get("msfem.localproblemsolver_precision", 1e-5);
  options["verbose"] = problem_.config().get("msfem.local_solver_verbose", "0");
  {
    InverseOperatorResult stat;
    auto writable_rhs = current_rhs.vector().copy();
#if HAVE_UMFPACK
    if (use_umfpack_)
      local_direct_inverse_->apply(current_solution.vector().backend(), writable_rhs.backend(), stat);
    else
#endif
      local_inverse.apply(current_solution.vector(), writable_rhs, options);
  }
  if (!current_solution.dofs_valid())
    DUNE_THROW(Dune::InvalidStateException, "Current solution of the local msfem problem invalid!");
}

} // namespace Multiscale {
} // namespace Dune {
