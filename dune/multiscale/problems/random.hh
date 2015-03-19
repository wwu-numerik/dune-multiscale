// dune-multiscale
// Copyright Holders: Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_RANDOM
#define DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_RANDOM

#include <dune/multiscale/problems/base.hh>
#include <dune/stuff/grid/boundaryinfo.hh>
#include <memory>
#include <string>

#include "dune/multiscale/common/traits.hh"

#ifdef __GNUC__

#define PURE __attribute__((const))
#define HOT __attribute__((hot))
#define ALWAYS_INLINE __attribute__((always_inline)) inline

#else

#define PURE
#define HOT
#define ALWAYS_INLINE inline

#endif

template< int DIM, typename X, typename R, typename COR >
class Permeability;

namespace Dune {
namespace Multiscale {
namespace Problem {
/** \addtogroup problem_9 Problem::Nine
 * @{ **/
//! ------------ Elliptic Problem 9 -------------------

namespace Random {

/// Class defining exponential correlation function.
/// \tparam DIM   dimension
/// \tparam X     point type
/// \tparam R     scalar type
///
/// \author jan.mohring@itwm.fraunhofer.de
/// \date 2014
class Correlation {
  static constexpr auto DIM = CommonTraits::world_dim;
  typedef DomainType X;
  typedef CommonTraits::DomainFieldType R;
public:
  /// Constructor
  /// \param corrLen   correlation length
  /// \param sigma     standard deviation
  Correlation(R corrLen=0.1, R sigma=1.0)
    : _corrLen(corrLen), _sigma2(sigma*sigma) {}

  Correlation(const Correlation& old) {
    _corrLen = old._corrLen;
    _sigma2  = old._sigma2;
  }

  /// Evaluation
  /// \param d   difference of points to take corretation of
  R operator() (X d) const {
    R sumX2 = 0;
    for(int i=0; i<DIM; ++i) {
       sumX2 += d[i]*d[i];
    }
    return _sigma2 * exp(-sqrt(sumX2)/_corrLen);
  }

private:
  R _corrLen;   //< correlation length
  R _sigma2;    //< standard deviation
};

struct ModelProblemData : public IModelProblemData {
  virtual bool hasExactSolution() const final override { return false; }

  ModelProblemData();

  std::string getMacroGridFile() const final override;
  const BoundaryInfoType& boundaryInfo() const final override;
  const SubBoundaryInfoType& subBoundaryInfo() const final override;
  std::pair<CommonTraits::DomainType, CommonTraits::DomainType> gridCorners() const final override;

  virtual void problem_init(MPIHelper::MPICommunicator global, MPIHelper::MPICommunicator local)  final override;
  virtual void prepare_new_evaluation()  final override;

private:
  Dune::ParameterTree boundary_settings() const;
  std::unique_ptr<DSG::BoundaryInfos::NormalBased<typename View::Intersection>> boundaryInfo_;
  DSG::BoundaryInfos::AllDirichlet<typename SubView::Intersection> subBoundaryInfo_;
};

class Diffusion : public DiffusionBase {
public:
  Diffusion();

  //! currently used in gdt assembler
  virtual void evaluate(const DomainType& x, DiffusionBase::RangeType& y) const final override;
  PURE HOT void diffusiveFlux(const DomainType& x, const Problem::JacobianRangeType& direction,
                              Problem::JacobianRangeType& flux) const final override;

  virtual size_t order() const final override;

  virtual void init(MPIHelper::MPICommunicator global, MPIHelper::MPICommunicator local)  final override;
  virtual void prepare_new_evaluation()  final override;

private:
  typedef Permeability<CommonTraits::world_dim, DomainType, CommonTraits::DomainFieldType, Correlation> PermeabilityType;
  std::unique_ptr<Correlation> correlation_;
#if HAVE_RANDOM_PROBLEM
  std::unique_ptr<PermeabilityType> field_;
#endif
};

class DirichletData : public DirichletDataBase {
public:
  DirichletData() {}
  virtual ~DirichletData() {}

  PURE void evaluate(const DomainType& x, RangeType& y) const final override;
};

class NeumannData : public NeumannDataBase {
public:
  NeumannData() {}

  PURE void evaluate(const DomainType& x, RangeType& y) const final override;
};

MSNULLFUNCTION(ExactSolution)
MSNULLFUNCTION(Source)

} //! @} namespace Random {
}
} // namespace Multiscale {
} // namespace Dune {

#endif // ifndef DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_NINE
