// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_FEM_TYPES_HH
#define DUNE_FEM_TYPES_HH

#include <dune/multiscale/common/la_backend.hh>
#include <dune/multiscale/common/traits.hh>
#include <dune/common/tuples.hh>

#include <dune/multiscale/common/righthandside_assembler.hh>

#include <dune/fem/space/common/adaptmanager.hh>
#include <dune/fem/io/file/dataoutput.hh>

namespace Dune {
namespace Multiscale {
namespace FEM {

//! Type constructions for the FEM problem
struct FEMTraits {
  typedef typename BackendChooser<typename CommonTraits::DiscreteFunctionSpaceType>::InverseOperatorType
  InverseOperatorType;
}; // struct  FEMTraits

} // namespace FEM {
} // namespace Multiscale {
} // namespace Dune {

#endif // DUNE_FEM_TYPES_HH
