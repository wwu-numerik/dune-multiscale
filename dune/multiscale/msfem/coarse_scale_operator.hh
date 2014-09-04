// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_MULTISCALE_MSFEM_COARSESCALE_OPERATOR_HH
#define DUNE_MULTISCALE_MSFEM_COARSESCALE_OPERATOR_HH

#include <ostream>
#include <type_traits>
#include <assert.h>
#include <boost/noncopyable.hpp>

#include <dune/common/fmatrix.hh>
#include <dune/stuff/common/logging.hh>
#include <dune/stuff/fem/functions/checks.hh>
#include <dune/stuff/fem/localmatrix_proxy.hh>
#include <dune/stuff/fem/functions/integrals.hh>
#include <dune/gdt/assembler/system.hh>
#include <dune/multiscale/common/traits.hh>
#include <dune/multiscale/msfem/msfem_traits.hh>
#include <dune/multiscale/msfem/coarse_scale_assembler.hh>

namespace Dune {
namespace Multiscale {
namespace MsFEM {

class MsFEMCodim0Integral;
class MsFemCodim0Matrix;
class LocalSolutionManager;
class LocalGridList;
class CoarseScaleOperator;

class CoarseScaleOperatorTraits {
public:
  typedef CoarseScaleOperator derived_type;
  typedef CommonTraits::LinearOperatorType MatrixType;
  typedef CommonTraits::GdtSpaceType SourceSpaceType;
  typedef CommonTraits::GdtSpaceType RangeSpaceType;
  typedef CommonTraits::GridViewType GridViewType;
}; // class EllipticCGTraits

class CoarseScaleOperator
    : public GDT::Operators::MatrixBased<CoarseScaleOperatorTraits>,
      public GDT::SystemAssembler<CoarseScaleOperatorTraits::RangeSpaceType, CoarseScaleOperatorTraits::GridViewType,
                                  CoarseScaleOperatorTraits::SourceSpaceType> {
  typedef GDT::Operators::EllipticCG<Problem::DiffusionBase, CommonTraits::LinearOperatorType,
                                     CommonTraits::GdtSpaceType> EllipticOperatorType;
  typedef GDT::Operators::MatrixBased<CoarseScaleOperatorTraits> OperatorBaseType;
  typedef MsFEMCodim0Integral LocalOperatorType;
  typedef MsFemCodim0Matrix LocalAssemblerType;
  typedef GDT::SystemAssembler<CoarseScaleOperatorTraits::RangeSpaceType, CoarseScaleOperatorTraits::GridViewType,
                               CoarseScaleOperatorTraits::SourceSpaceType> AssemblerBaseType;
  typedef CommonTraits::DiscreteFunctionType CoarseDiscreteFunction;
  typedef typename CommonTraits::DiscreteFunctionSpaceType CoarseDiscreteFunctionSpace;

public:
  typedef CoarseScaleOperatorTraits Traits;

  typedef typename Traits::MatrixType MatrixType;
  typedef typename Traits::SourceSpaceType SourceSpaceType;
  typedef typename Traits::RangeSpaceType RangeSpaceType;
  typedef typename Traits::GridViewType GridViewType;

  using OperatorBaseType::pattern;

  static Stuff::LA::SparsityPatternDefault pattern(const RangeSpaceType& range_space,
                                                   const SourceSpaceType& source_space, const GridViewType& grid_view);

  CoarseScaleOperator(const SourceSpaceType& coarse_space, LocalGridList& localGridList);

  virtual ~CoarseScaleOperator() {}

  virtual void assemble() DS_OVERRIDE DS_FINAL;

  void apply_inverse(CoarseScaleOperator::CoarseDiscreteFunction& solution);

  MatrixType& system_matrix();
  const MatrixType& system_matrix() const;


private:
  const SourceSpaceType& coarse_space() const;

  MatrixType global_matrix_;
  const LocalOperatorType local_operator_;
  const LocalAssemblerType local_assembler_;
  CommonTraits::DiscreteFunctionType msfem_rhs_;
  CommonTraits::DiscreteFunctionType dirichlet_projection_;
}; // class CoarseScaleOperator

} // namespace MsFEM {
} // namespace Multiscale {
} // namespace Dune {

#endif // #ifndef DUNE_MULTISCALE_MSFEM_COARSESCALE_OPERATOR_HH
