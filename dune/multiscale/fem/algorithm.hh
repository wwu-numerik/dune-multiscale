// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_FEM_ALGORITHM_HH
#define DUNE_FEM_ALGORITHM_HH

#include <dune/multiscale/common/traits.hh>
#include <dune/multiscale/fem/fem_traits.hh>


namespace Dune {
namespace Multiscale {
namespace FEM {


//! write discrete function to a file + VTK Output
void write_discrete_function(typename CommonTraits::DiscreteFunctionType& discrete_solution );

//! \TODO docme
void solve(typename CommonTraits::DiscreteFunctionType& solution,
           const typename CommonTraits::DiscreteFunctionSpaceType& finerDiscreteFunctionSpace,
           const typename FEMTraits::EllipticOperatorType& discrete_elliptic_op,
           const std::string& filename,
           const Dune::RightHandSideAssembler< typename CommonTraits::DiscreteFunctionType >& rhsassembler);

//! the main FEM computation
void algorithm(typename CommonTraits::GridPointerType& macro_grid_pointer,
               const std::string filename);

//! \TODO docme
void algorithm_hom_fem(typename CommonTraits::GridPointerType& macro_grid_pointer,
                       const std::string filename);

} //namespace FEM {
} //namespace Multiscale {
} //namespace Dune {

#endif // DUNE_FEM_ALGORITHM_HH
