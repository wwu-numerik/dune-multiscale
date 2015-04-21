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
#include <dune/gdt/assembler/system.hh>
#include <dune/multiscale/common/traits.hh>
#include <dune/multiscale/msfem/msfem_traits.hh>
#include <dune/multiscale/msfem/coarse_scale_assembler.hh>

namespace Dune {
namespace Multiscale {

class MsFEMCodim0Integral;
class MsFemCodim0Matrix;
class LocalSolutionManager;
class LocalGridList;
class CoarseScaleOperator;

namespace Problem {
struct ProblemContainer;
}

class CoarseScaleOperatorTraits {
public:
  typedef CoarseScaleOperator derived_type;
  typedef CommonTraits::LinearOperatorType MatrixType;
  typedef CommonTraits::SpaceType SourceSpaceType;
  typedef CommonTraits::SpaceType RangeSpaceType;
  typedef CommonTraits::GridViewType GridViewType;
}; // class EllipticCGTraits

class CoarseScaleOperator
    : public GDT::Operators::MatrixBased<CoarseScaleOperatorTraits>,
      public GDT::SystemAssembler<CoarseScaleOperatorTraits::RangeSpaceType, CommonTraits::InteriorGridViewType,
                                  CoarseScaleOperatorTraits::SourceSpaceType> {
  typedef GDT::Operators::EllipticCG<Problem::DiffusionBase, CommonTraits::LinearOperatorType, CommonTraits::SpaceType>
      EllipticOperatorType;
  typedef GDT::Operators::MatrixBased<CoarseScaleOperatorTraits> OperatorBaseType;
  typedef MsFEMCodim0Integral LocalOperatorType;
  typedef MsFemCodim0Matrix LocalAssemblerType;
  typedef GDT::SystemAssembler<CoarseScaleOperatorTraits::RangeSpaceType, CommonTraits::InteriorGridViewType,
                               CoarseScaleOperatorTraits::SourceSpaceType> AssemblerBaseType;
  typedef CommonTraits::DiscreteFunctionType CoarseDiscreteFunction;
  typedef typename CommonTraits::SpaceType CoarseDiscreteFunctionSpace;

public:
  typedef CoarseScaleOperatorTraits Traits;

  typedef typename Traits::MatrixType MatrixType;
  typedef typename Traits::SourceSpaceType SourceSpaceType;
  typedef typename Traits::RangeSpaceType RangeSpaceType;
  typedef typename Traits::GridViewType GridViewType;

  using OperatorBaseType::pattern;

  static Stuff::LA::SparsityPatternDefault pattern(const RangeSpaceType& range_space,
                                                   const SourceSpaceType& source_space, const GridViewType& grid_view);

  CoarseScaleOperator(DMP::ProblemContainer &problem, const SourceSpaceType& source_space_in, LocalGridList& localGridList);

  virtual ~CoarseScaleOperator() {}

  virtual void assemble() override final;

  void apply_inverse(CoarseScaleOperator::CoarseDiscreteFunction& solution);

private:
  //! used as an alias to test_space()
  const SourceSpaceType& coarse_space() const;

  MatrixType global_matrix_;
  const LocalOperatorType local_operator_;
  const LocalAssemblerType local_assembler_;
  CommonTraits::DiscreteFunctionType msfem_rhs_;
  CommonTraits::DiscreteFunctionType dirichlet_projection_;
  DMP::ProblemContainer& problem_;
}; // class CoarseScaleOperator

} // namespace Multiscale {
} // namespace Dune {

#endif // #ifndef DUNE_MULTISCALE_MSFEM_COARSESCALE_OPERATOR_HH
