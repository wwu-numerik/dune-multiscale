// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DiscreteEllipticHMMOperator_HH
#define DiscreteEllipticHMMOperator_HH

#include <boost/noncopyable.hpp>
#include <memory>

#include <dune/multiscale/hmm/hmm_traits.hh>
#include <dune/common/fmatrix.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/operator/2order/lagrangematrixsetup.hh>

namespace Dune {
namespace Multiscale {
namespace HMM {

class CellProblemNumberingManager;

//! \TODO docme
class DiscreteEllipticHMMOperator
    : public Operator<typename CommonTraits::DiscreteFunctionType::RangeFieldType,
                      typename CommonTraits::DiscreteFunctionType::RangeFieldType, CommonTraits::DiscreteFunctionType,
                      CommonTraits::DiscreteFunctionType>,
      boost::noncopyable {

private:
  typedef CommonTraits::DiscreteFunctionType DiscreteFunction;
  typedef HMMTraits::PeriodicDiscreteFunctionType PeriodicDiscreteFunction;
  typedef CommonTraits::DiffusionType DiffusionModel;
  typedef typename DiscreteFunction::DiscreteFunctionSpaceType DiscreteFunctionSpace;
  typedef typename PeriodicDiscreteFunction::DiscreteFunctionSpaceType PeriodicDiscreteFunctionSpace;
  typedef typename DiscreteFunctionSpace::RangeFieldType RangeFieldType;
  typedef typename DiscreteFunctionSpace::DomainType DomainType;
  typedef typename DiscreteFunctionSpace::RangeType RangeType;
  typedef typename DiscreteFunctionSpace::JacobianRangeType JacobianRangeType;  
  typedef typename DiscreteFunctionSpace::BasisFunctionSetType BaseFunctionSet;
  typedef typename DiscreteFunctionSpace::EntityType Entity;

  static const int dimension = DiscreteFunctionSpace::GridType::dimension;
  static const int polynomialOrder = DiscreteFunctionSpace::polynomialOrder;

public:
  DiscreteEllipticHMMOperator(const DiscreteFunctionSpace& discreteFunctionSpace,
                              const PeriodicDiscreteFunctionSpace& periodicDiscreteFunctionSpace,
                              const DiffusionModel& diffusion_op, const CellProblemNumberingManager& cp_num_manager)
    : discreteFunctionSpace_(discreteFunctionSpace)
    , periodicDiscreteFunctionSpace_(periodicDiscreteFunctionSpace)
    , diffusion_operator_(diffusion_op)
    , cp_num_manager_(cp_num_manager) {}

private:
  // dummy operator
  virtual void operator()(const DiscreteFunction& u, DiscreteFunction& w) const;
  void boundary_treatment(CommonTraits::FEMMatrix& global_matrix) const;

public:
  void assemble_matrix(CommonTraits::FEMMatrix& global_matrix) const;
  void assemble_jacobian_matrix(DiscreteFunction& old_macro_function, CommonTraits::FEMMatrix& global_matrix) const;

private:
  const DiscreteFunctionSpace& discreteFunctionSpace_;
  const PeriodicDiscreteFunctionSpace& periodicDiscreteFunctionSpace_;
  const DiffusionModel& diffusion_operator_;
  const CellProblemNumberingManager& cp_num_manager_;

  // name of data file, e.g. required if we want to use the saved solutions of the cell problems
  const std::string filename_;
};

} // namespace HMM {
} // namespace Multiscale {
} // namespace Dune {

#endif // #ifndef DiscreteElliptic_HH
