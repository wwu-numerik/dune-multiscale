// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

// das alles klappt momentan nur in 2-D!!! Es laesst sich aber sehr einfach verallgemeinern.
// Implementierung funktioniert unter folgenden Voraussetzungen an den (Diffusions-)Tensor:
// 1. Elliptizitätsbedingung
// 2. A(x,y) = A(y)

#ifndef DUNE_HOMOGENIZER_HH
#define DUNE_HOMOGENIZER_HH

#include <dune/multiscale/fem/elliptic_fem_matrix_assembler.hh>
#include <dune/fem/gridpart/periodicgridpart/periodicgridpart.hh>

// for data output:
#include <dune/fem/io/file/dataoutput.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/io/file/datawriter.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/stuff/fem/functions/analytical.hh>
#include <dune/stuff/common/memory.hh>
#include <dune/multiscale/common/righthandside_assembler.hh>

namespace Dune {
namespace Multiscale {

//! define output traits
struct CellDataOutputParameters : public Fem::DataOutputParameters {
public:
  std::string my_prefix_;

  // path where the data is stored
  std::string path() const { return "data_output_hmm"; }

  void set_prefix(std::string my_prefix) {
    my_prefix_ = my_prefix;
    // std :: cout << "Set prefix. my_prefix_ = " << my_prefix_ << std :: endl;
  }

  // base of file name for data file
  std::string prefix() const {
    if (my_prefix_ == "")
      return "solutions_cell_problems";
    else
      return my_prefix_;
  }

  // format of output:
  int outputformat() const {
    // return 0; // GRAPE (lossless format)
    return 1; // VTK
              // return 2; // VTK vertex data
              // return 3; // gnuplot
  }
};


//! \TODO docme
// (to replace the more general lower order term)
class MassWeight : public Problem::LowerOrderBase {
public:

  MassWeight(double lambda) : lambda_(lambda) {}

  void evaluate(const DomainType& /*x*/, const RangeType& position, const JacobianRangeType& /*direction_gradient*/,
                RangeType& y) const {
    y = lambda_ * position;
  }

  virtual void evaluate(const DomainType& /*x*/, RangeType& /*ret*/) const { DUNE_THROW(Dune::NotImplemented, ""); }
  virtual void evaluate(const DomainType&, const TimeType&, RangeType& ) const { DUNE_THROW(Dune::NotImplemented, ""); }

  virtual void position_derivative(const DomainType&, const RangeType&, const JacobianRangeType&, RangeType& y) const {
    y = RangeType(0);
  }

  virtual void direction_derivative(const DomainType&, const RangeType&, const JacobianRangeType&,
                                    JacobianRangeType& y) const {
    y = JacobianRangeType(0);
  }

  double lambda_;
};

/** since we need to evaluate A( x, \cdot ) to solve cellproblems (in comparison to A( \cdot, \frac{\cdot}{\epsilon} )
 * for the global problem), we must transform the orginal tensor to be able to use a standard FEM operator for solving
 * cell problems (otherwise: calling the method evaluate(i,j,x,y) within the matrixassembler would evaluate A^{\epsilon}
 * instead of A(x,\cdot) )
 **/

class TransformTensor : public Problem::DiffusionBase {
private:
  typedef CommonTraits::FunctionSpaceType FunctionSpaceType;
  typedef typename CommonTraits::DiffusionType TensorType;

  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType RangeType;
  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

  typedef DomainFieldType TimeType;

  static const int dimDomain = DomainType::dimension;

  const TensorType& tensor_;

public:
  // Constructor
  inline explicit TransformTensor(const TensorType& tensor) : tensor_(tensor) {}

  void diffusiveFlux(const DomainType& y, const JacobianRangeType& direction, JacobianRangeType& flux) const;

  inline void evaluate(const int /*i*/, const int /*j*/, const DomainType& /*x*/, const TimeType& /*time*/,
                       RangeType& /*y*/) const;

  inline void evaluate(const DomainType& /*x*/, RangeType& /*y*/) const { DUNE_THROW(Dune::NotImplemented, ""); }

  inline void evaluate(const DomainType& /*x*/, const TimeType& /*time*/, RangeType& /*y*/) const;

  virtual void jacobianDiffusiveFlux(const DomainType& /*x*/, const JacobianRangeType& /*position_gradient*/,
                                     const JacobianRangeType& /*direction_gradient*/,
                                     JacobianRangeType& /*flux*/) const;
};


//! \TODO docme

class Homogenizer {
private:
  typedef CommonTraits::GridType GridType;
  static const int dimension = GridType::dimension;
  typedef typename CommonTraits::DiffusionType TensorType;
  typedef Fem::PeriodicLeafGridPart<GridType> PeriodicGridPartType;
  typedef CommonTraits::FunctionSpaceType FunctionSpaceType;
  typedef Fem::LagrangeDiscreteFunctionSpace<FunctionSpaceType, PeriodicGridPartType, 1>
  PeriodicDiscreteFunctionSpaceType;

  // to avoid confusion:
  typedef PeriodicDiscreteFunctionSpaceType DummySpaceType;
  // (sometimes PeriodicDiscreteFunctionSpaceType is only a dummy)
  typedef Fem::AdaptiveDiscreteFunction<PeriodicDiscreteFunctionSpaceType> PeriodicDiscreteFunctionType;
  // to avoid confusion:
  typedef PeriodicDiscreteFunctionType DummyType;
  // (sometimes PeriodicDiscreteFunctionType is only a dummy)
  typedef MassWeight MassWeightType;
  typedef typename PeriodicDiscreteFunctionSpaceType::JacobianRangeType PeriodicJacobianRangeType;
  typedef typename PeriodicDiscreteFunctionType::LocalFunctionType PeriodicLocalFunctionType;
  typedef typename GridType::Codim<0>::Entity EntityType;

  typedef typename BackendChooser<PeriodicDiscreteFunctionSpaceType>::LinearOperatorType LinearOperatorType;
  typedef typename BackendChooser<PeriodicDiscreteFunctionSpaceType>::InverseOperatorType InverseLinearOperatorType;

  // discrete elliptic operator (corresponds with FEM Matrix)
  typedef Multiscale::FEM::DiscreteEllipticOperator<PeriodicDiscreteFunctionType, TransformTensor>
  EllipticOperatorType;

  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType RangeType;

  static const int spacePolOrd = PeriodicDiscreteFunctionSpaceType::polynomialOrder;

  // dgf file that describes the perforated domain
  const std::string& filename_;

public:
  Homogenizer(const std::string& filename);

  typedef FieldMatrix<RangeType, dimension, dimension> HomTensorType;

private:
  double getEntry(const TransformTensor& tensor,
                  const PeriodicDiscreteFunctionSpaceType& periodicDiscreteFunctionSpace,
                  const PeriodicDiscreteFunctionType& w_i, const PeriodicDiscreteFunctionType& w_j, const int& i,
                  const int& j) const; // end of method

public:
  HomTensorType getHomTensor(const TensorType& tensor) const; // getHomTensor
};  // end of class

} // namespace Multiscale {
} // end namespace
#endif // ifndef DUNE_HOMOGENIZER_HH
