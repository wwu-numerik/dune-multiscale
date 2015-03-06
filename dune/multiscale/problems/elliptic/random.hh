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
template<int DIM, typename X, typename R>
class Correlation {
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
  virtual bool hasExactSolution() const { return false; }

  ModelProblemData();

  std::string getMacroGridFile() const;
  const BoundaryInfoType& boundaryInfo() const;
  const SubBoundaryInfoType& subBoundaryInfo() const;
  std::pair<CommonTraits::DomainType, CommonTraits::DomainType> gridCorners() const;

  virtual void problem_init(MPIHelper::MPICommunicator global, MPIHelper::MPICommunicator local);

private:
  Dune::ParameterTree boundary_settings() const;
  std::unique_ptr<DSG::BoundaryInfos::NormalBased<typename View::Intersection>> boundaryInfo_;
  DSG::BoundaryInfos::AllDirichlet<typename SubView::Intersection> subBoundaryInfo_;
};

class Source : public Dune::Multiscale::CommonTraits::FunctionBaseType {
public:
  Source();

  PURE HOT void evaluate(const DomainType& x, RangeType& y) const;
  virtual size_t order() const;
};

class Diffusion : public DiffusionBase {
public:
  Diffusion();

  //! currently used in gdt assembler
  virtual void evaluate(const DomainType& x, DiffusionBase::RangeType& y) const;
  PURE HOT void diffusiveFlux(const DomainType& x, const Problem::JacobianRangeType& direction,
                              Problem::JacobianRangeType& flux) const;

  virtual size_t order() const;

  virtual void init(MPIHelper::MPICommunicator global, MPIHelper::MPICommunicator local) override;
  virtual void prepare_new_evaluation() override;

private:
  typedef Correlation<CommonTraits::world_dim, DomainType, double> CorrelationType;
  typedef Permeability<CommonTraits::world_dim, DomainType, double, CorrelationType> PermeabilityType;
  std::unique_ptr<CorrelationType> correlation_;
#if HAVE_RANDOM_PROBLEM
  std::unique_ptr<PermeabilityType> field_;
#endif
};

class DirichletData : public DirichletDataBase {
public:
  DirichletData() {}
  virtual ~DirichletData() {}

  PURE void evaluate(const DomainType& x, RangeType& y) const;
  PURE void jacobian(const DomainType& x, JacobianRangeType& y) const;
};

class NeumannData : public NeumannDataBase {
public:
  NeumannData() {}

  PURE void evaluate(const DomainType& x, RangeType& y) const;
};

MSNULLFUNCTION(DirichletBoundaryCondition)
MSNULLFUNCTION(NeumannBoundaryCondition)
MSNULLFUNCTION(ExactSolution)

} //! @} namespace Nine {
}
} // namespace Multiscale {
} // namespace Dune {

#endif // ifndef DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_NINE
