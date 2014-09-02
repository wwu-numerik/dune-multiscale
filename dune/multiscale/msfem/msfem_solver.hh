// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef Elliptic_MSEM_Solver_HH
#define Elliptic_MSEM_Solver_HH

#include <dune/common/fmatrix.hh>
#include <dune/common/typetraits.hh>
#include <dune/multiscale/common/traits.hh>
#include <dune/multiscale/msfem/msfem_traits.hh>
#include <dune/stuff/discretefunction/projection/heterogenous.hh>
#include <dune/stuff/fem/functions/checks.hh>

#include <dune/multiscale/common/la_backend.hh>

namespace Dune {
namespace Multiscale {

namespace Problem {
struct DiffusionBase;
}
namespace MsFEM {

class LocalGridList;

//! \TODO needs a better name
class Elliptic_MsFEM_Solver {
private:
  typedef CommonTraits::DiscreteFunctionType DiscreteFunctionType;

  static const int faceCodim = 1;

  typedef MsFEMTraits::LocalGridDiscreteFunctionSpaceType LocalGridDiscreteFunctionSpaceType;
  typedef MsFEMTraits::LocalGridDiscreteFunctionType LocalGridDiscreteFunctionType;

  //! identify fine scale part of MsFEM solution (including the projection!)
  void identify_fine_scale_part(LocalGridList& subgrid_list, const DiscreteFunctionType& coarse_msfem_solution,
                                DiscreteFunctionType& fine_scale_part) const;

public:
  /** - ∇ (A(x,∇u)) + b ∇u + c u = f - divG
   then:
   A --> diffusion operator ('DiffusionOperatorType')
   b --> advective part ('AdvectionTermType')
   c --> reaction part ('ReactionTermType')
   f --> 'first' source term, scalar ('SourceTermType')
   homogenous Dirchilet boundary condition!:
   **/
  void apply(const CommonTraits::DiscreteFunctionSpaceType& coarse_space, const Problem::DiffusionBase& diffusion_op,
             DiscreteFunctionType& coarse_scale_part, DiscreteFunctionType& fine_scale_part,
             DiscreteFunctionType& solution) const;
};

} // namespace MsFEM {
} // namespace Multiscale {
} // namespace Dune {

#endif // #ifndef Elliptic_MSEM_Solver_HH
