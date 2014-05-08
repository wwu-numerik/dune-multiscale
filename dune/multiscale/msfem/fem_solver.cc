#include <config.h>

#include <dune/common/exceptions.hh>
#include <dune/common/timer.hh>
#include <dune/fem/function/common/function.hh>
#include <dune/multiscale/common/righthandside_assembler.hh>
#include <dune/multiscale/common/traits.hh>
#include <dune/multiscale/fem/elliptic_fem_matrix_assembler.hh>
#include <dune/multiscale/fem/fem_traits.hh>
#include <dune/stuff/common/logging.hh>
#include <dune/stuff/common/profiler.hh>
#include <dune/stuff/functions/norm.hh>
#include <dune/stuff/common/parameter/configcontainer.hh>
#include <dune/pdelab/finiteelementmap/qkfem.hh>
#include <dune/pdelab/constraints/common/constraints.hh>
#include <dune/pdelab/constraints/conforming.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridoperator/gridoperator.hh>
#include <dune/pdelab/backend/istlsolverbackend.hh>
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

Elliptic_FEM_Solver::Elliptic_FEM_Solver(const CommonTraits::GdtSpaceType &space)
  : space_(space) {}

void Elliptic_FEM_Solver::apply(const CommonTraits::DiffusionType& diffusion_op,
                                const CommonTraits::SourceType& f,
                                CommonTraits::GdtDiscreteFunctionType &solution) const {
  DSC_LOG_DEBUG << "Solving linear problem with standard FEM\n";
  DSC_PROFILER.startTiming("fem.apply");

  DUNE_THROW(NotImplemented, "");

  DSC_PROFILER.stopTiming("fem.apply");
  DSC_LOG_DEBUG << "Standard FEM problem solved in " << DSC_PROFILER.getTiming("fem.apply") << "ms.\n";
}

} // namespace Multiscale
}
